import dnaio
from xopen import xopen
import os
import traceback
import io
from multiprocess import Process, Pipe, Queue
from typing import BinaryIO
from qualtrim_new import quality_trim_index
from ultraplex.modifiers import AdapterCutter, ModificationInfo
from ultraplex.adapters import BackAdapter
import gzip
import glob
import time
import argparse
import shutil
from pathlib import Path
import logging
from math import log10, floor
from collections import defaultdict
from ._version import __version__


def merge_defaultdicts(d1, d2):
    merged = defaultdict(int)
    for key, value in d1.items():
        merged[key] += value
    for key, value in d2.items():
        merged[key] += value
    return merged


# for logging
def setup_logger(name, log_file, level=logging.INFO):
    """To setup as many loggers as you want"""

    formatter = logging.Formatter("%(asctime)-15s %(levelname)-8s %(message)s")
    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger


def round_sig(x, sig=2):
    try:
        to_return = round(x, sig - int(floor(log10(abs(x)))) - 1)
    except:
        to_return = 0
    return to_return


def make_5p_bc_dict(barcodes, min_score, dont_build_reference):
    """
    this function generates a dictionary that matches each possible sequence
    from the read with the best 5' barcode
    """
    if dont_build_reference:
        return {"dont_build": True}
    else:
        first_bc = barcodes[0]
        seq_length = len(first_bc.replace("N", ""))

        # check all the barcodes are the same length (ignoring UMIs)
        for bc in barcodes:
            assert (
                len(bc.replace("N", "")) == seq_length
            ), "Your experimental barcodes are different lengths."

        seqs = make_all_seqs(seq_length)

        # trim sequences to desired length
        # create dictionary
        barcode_dictionary = {}

        for seq in seqs:
            barcode_dictionary[seq] = score_barcode_for_dict(seq, barcodes, min_score)
        return barcode_dictionary


def remove_Ns_from_barcodes(barcodes):
    # returns a dictionary with keys of the N-removed barcode, and values of the original barcode
    barcodes_no_N = {a.replace("N", ""): a for a in barcodes}
    barcodes_no_N["no_match"] = "no_match"
    return barcodes_no_N


def score_barcode_for_dict(seq, barcodes, min_score, Ns_removed=False):
    """
    this function scores a given sequence against all the barcodes. It returns the winner with Ns included.
    """

    if not Ns_removed:
        barcodes = remove_Ns_from_barcodes(barcodes)
    barcodes_no_N = barcodes.keys()

    if seq in barcodes_no_N:  # no need to check all barcodes
        winner = barcodes[seq]  # barcode WITH Ns included
    elif min_score == len(
        barcodes_no_N
    ):  # i.e. no matches allowed, and seq not in barcodes
        winner = "no_match"
    else:  # mismatches allowed so need to check
        scores = {}

        for this_bc in barcodes_no_N:
            # score the barcode against the read, penalty for N in the read
            score = sum(a == b for a, b in zip(this_bc, seq))
            scores[this_bc] = score

        # Find the best score
        best_score = max(scores.values())

        if best_score < min_score:
            winner = "no_match"
        else:
            # check that there is only one barcode with the max score
            filtered = [a for a, b in scores.items() if b == best_score]
            if len(filtered) > 1:
                winner = "no_match"
            else:  # if there is only one
                winner = barcodes[filtered[0]]  # barcode WITH Ns included
    return winner


class ReaderProcess(Process):
    """
    Read chunks of FASTA or FASTQ data (single-end or paired) and send to a worker.

    The reader repeatedly

    - reads a chunk from the file(s)
    - reads a worker index from the Queue
    - sends the chunk to connections[index]

    and finally sends the stop token -1 ("poison pills") to all connections.
    """

    def __init__(self, file, connections, queue, buffer_size, i2):
        # /# Setup the reader process

        """
        queue -- a Queue of worker indices. A worker writes its own index into this
                queue to notify the reader that it is ready to receive more data.
        connections -- a list of Connection objects, one for each worker.
        """
        super().__init__()
        self.file = file
        self.connections = connections
        self.queue = queue
        self.buffer_size = buffer_size
        self.file2 = i2

    def run(self):
        # if self.stdin_fd != -1:
        # 	sys.stdin.close()
        # 	sys.stdin = os.fdopen(self.stdin_fd) #/# not sure why I need this!

        try:
            with xopen(self.file, "rb") as f:
                if self.file2:
                    with xopen(self.file2, "rb") as f2:
                        for chunk_index, (chunk1, chunk2) in enumerate(
                            dnaio.read_paired_chunks(f, f2, self.buffer_size)
                        ):
                            self.send_to_worker(chunk_index, chunk1, chunk2)
                else:
                    for chunk_index, chunk in enumerate(
                        dnaio.read_chunks(f, self.buffer_size)
                    ):
                        self.send_to_worker(chunk_index, chunk)

            # Send poison pills to all workers
            for _ in range(len(self.connections)):
                worker_index = self.queue.get()
                self.connections[worker_index].send(-1)
        except Exception as e:
            # TODO better send this to a common "something went wrong" Queue
            for connection in self.connections:
                connection.send(-2)
                connection.send((e, traceback.format_exc()))

    def send_to_worker(self, chunk_index, chunk, chunk2=None):
        worker_index = self.queue.get()  # get a worker that needs work
        connection = self.connections[
            worker_index
        ]  # find the connection to this worker
        connection.send(chunk_index)  # send the index of this chunk to the worker
        connection.send_bytes(chunk)  # /# send the actual data to this worker
        if chunk2 is not None:
            connection.send_bytes(chunk2)


class InputFiles:
    """
    this is from cutadapt - basically just creates a dnaio object
    """

    def __init__(self, file1: BinaryIO, interleaved: bool = False):
        self.file1 = file1
        self.interleaved = interleaved

    def open(self):
        return dnaio.open(self.file1, interleaved=self.interleaved, mode="r")


def make_3p_bc_dict(bcs, min_score):
    """
    Generates a dictionary that matches a given sequence with
    its best 3' barcode match, or "no_match"
    """
    # First, find positions of barcode bases (not UMIs) relative to 3' end of read

    bcl = len(
        [a for a, b in enumerate(bcs[0]) if b != "N"]
    )  # number of non-N characters in barcode
    all_seqs = make_all_seqs(bcl)
    three_p_match_d = {}
    for seq in all_seqs:
        # now check if there actually is a match
        three_p_match_d[seq] = score_barcode_for_dict(seq, bcs, min_score)
    return three_p_match_d


def make_all_seqs(l):
    """
    Makes all possible sequences, including Ns, of length l
    """
    if l > 8:
        print("Warning - large barcodes detected!")
        print("It may be faster to use option '--dont_build_reference'!")

    nts = ["A", "C", "G", "T", "N"]

    all_seqs = nts

    for i in range(l - 1):
        new_seqs = []
        for seq in all_seqs:
            for nt in nts:
                new_seqs.append(seq + nt)
        all_seqs = new_seqs

    return all_seqs


def careful_rev_c(seq):
    if "no_match".upper() in seq.upper():
        return seq
    else:
        return rev_c(seq)


def rev_c(seq):
    """
    simple function that reverse complements a given sequence
    """
    tab = str.maketrans("ACTGN", "TGACN")
    # first reverse the sequence
    seq = seq[::-1]
    # and then complement
    seq = seq.translate(tab)
    return seq


def three_p_demultiplex(
    read,
    d,
    add_umi,
    linked_bcds,
    linked_bcds_no_N,
    min_score,
    keep_barcodes,
    reverse_complement=False,
):
    """
    read is a dnaio read
    d is the relevant match dictionary
    length is the length of the 3' barcodes (it assumes they're all the same)
    """

    if reverse_complement:  # if it's paired end reverse complement the reverse read
        sequence = rev_c(read.sequence)
        qualities = read.qualities[::-1]  # just reverse - can't complement qualities!
    else:
        sequence = read.sequence
        qualities = read.qualities

    # find the positions of the non-N characters relative to 3' end
    # TODO precalculate these for faster operation - make a dict of positions for each 5' bacordes linked 3' barcodes
    bc_length_reference = len(
        linked_bcds[0]
    )  # can use [0] because non-N guaranteed to be consistent
    non_N_poses = [
        a - bc_length_reference for a, b in enumerate(linked_bcds[0]) if b != "N"
    ]  # all negative

    seq_length = len(sequence)
    positions_to_extract = [x + seq_length for x in non_N_poses]

    umi = ""  # initialise

    if min(positions_to_extract) >= 0:  # ensure long enough
        bc = "".join(sequence[a] for a in positions_to_extract)

        if "dont_build" == d:
            assigned = score_barcode_for_dict(
                bc, linked_bcds_no_N, min_score, Ns_removed=True
            )
        else:
            assigned = d[bc]

        # if reverse_complement:  # the sequence that will be removed
        #     to_remove = read.sequence[0:len(assigned)]
        # else:
        #     # then it's not paired end, so we don't care
        #     to_remove = ""

        if assigned == "no_match":
            to_remove = ""

        else:  # ie there is a match - then we need to trim
            bc_length = len(assigned)

            # Find the positions in the read which contain the UMI
            umi_poses = [a - bc_length for a, b in enumerate(assigned) if b == "N"]
            umi_positions_to_extract = [x + seq_length for x in umi_poses]

            # check it's long enough to contain UMIs, or that there are no UMIs
            long_enough = False  # assume false

            if len(umi_poses) > 0:
                if min(umi_positions_to_extract) >= 0:
                    long_enough = True
            if len(umi_poses) == 0:  # no UMI
                long_enough = True

            # Extract the UMI
            if long_enough:
                umi = "".join(sequence[a] for a in umi_positions_to_extract)
                length = len(assigned)
                if not keep_barcodes:
                    sequence = sequence[0 : (len(read) - length)]
                    qualities = qualities[0 : (len(read.qualities) - length)]

                if reverse_complement:  # the sequence that will be removed
                    to_remove = read.sequence[0 : len(assigned)]
                else:
                    # then it's not paired end, so we don't care
                    to_remove = ""

            else:  # barcode assigned, but read too short to contain full umi
                assigned = "no_match"
                to_remove = ""

    else:  # read was too short for barcode to be extracted
        assigned = "no_match"
        to_remove = ""

    if add_umi:
        if "rbc:" in read.name:
            read.name = read.name + umi
        else:
            read.name = read.name + "rbc:" + umi

    if reverse_complement:
        # spin back round
        read.sequence = rev_c(sequence)
        read.qualities = qualities[::-1]
    else:
        read.sequence = sequence
        read.qualities = qualities

    return read, assigned, umi, to_remove


# TODO link three_p_mismatches to input
def make_dict_of_3p_bc_dicts(linked_bcds, three_p_mismatches, dont_build_reference):
    """
    this function makes a different dictionary for each 5' barcode
    it also checks that they're all the same length
    """

    d = {}
    min_score_d = {}
    for five_p_bc, three_p_bcs in linked_bcds.items():
        if len(three_p_bcs) > 0:  # ie this 5' barcode has 3' barcodes
            # check they're consistent
            check_N_position(three_p_bcs, "3")
            # work out the min score - only do first one because already checked they're consistent
            min_score = (
                len([a for a, b in enumerate(three_p_bcs[0]) if b != "N"])
                - three_p_mismatches
            )
            if dont_build_reference:
                d[five_p_bc] = "dont_build"
            else:
                this_dict = make_3p_bc_dict(three_p_bcs, min_score)
                d[five_p_bc] = this_dict
            min_score_d[five_p_bc] = min_score

    return d, min_score_d


class WorkerProcess(Process):  # /# have to have "Process" here to enable worker.start()
    """
    The worker repeatedly reads chunks of data from the read_pipe, runs the pipeline on it
    and sends the processed chunks to the write_pipe.

    To notify the reader process that it wants data, processes then writes out.
    """

    def __init__(
        self,
        index,
        read_pipe,
        need_work_queue,
        adapter,
        output_directory,
        five_p_bcs,
        three_p_bcs,
        save_name,
        total_demultiplexed,
        total_reads_assigned,
        total_reads_assigned_pass_length_filter,
        read_number,
        total_reads_qtrimmed,
        total_reads_adaptor_trimmed,
        total_reads_5p_no_3p,
        total_reads_5p_no_3p_pass_length_filter,
        ultra_mode,
        min_score_5_p,
        three_p_mismatches,
        linked_bcds,
        three_p_trim_q,
        q5,
        i2,
        adapter2,
        min_trim,
        ignore_no_match,
        final_min_length,
        dont_build_reference,
        keep_barcode,
        three_prime_only,
        TSO_umi_pos,
        TSO_end
    ):
        super().__init__()
        self._id = index  # the worker id
        self._read_pipe = read_pipe  # the pipe the reader reads data from
        self._need_work_queue = (
            need_work_queue  # worker adds its id to this queue when it needs work
        )
        self._end_qc = three_p_trim_q  # quality score to trim qc from 3' end
        self._start_qc = q5  # quality score to trim qc from 5' end
        self._total_demultiplexed = total_demultiplexed  # a queue which keeps track of the total number of reads processed
        self._total_reads_assigned = total_reads_assigned  # a queue which keeps track of the total number of reads assigned to samples
        self._total_reads_assigned_pass_length_filter = total_reads_assigned_pass_length_filter  # a queue which keeps track of the total number of reads assigned to sample files that pass length filter
        self._read_number = read_number  # a queue which has a defaultdict that keeps track of the total number of reads assigned to each sample
        self._total_reads_qtrimmed = total_reads_qtrimmed  # a queue which keeps track of the total number of reads quality trimmed
        self._total_reads_adaptor_trimmed = total_reads_adaptor_trimmed  # a queue which keeps track of the total number of reads adaptor trimmed
        self._total_reads_5p_no_3p = total_reads_5p_no_3p  # a queue which keeps track of how many reads have correct 5p BC but cannot find 3p BC
        self._total_reads_5p_no_3p_pass_length_filter = total_reads_5p_no_3p_pass_length_filter  # a queue which keeps track of how many reads have correct 5p BC but cannot find 3p BC that pass length filter
        self._adapter = adapter  # the 3' adapter
        self._final_min_length = final_min_length  # length of final written reads
        self._three_p_bcs = three_p_bcs  # not sure we need this? A list of all the 3' barcodes. But unecessary because of "linked_bcds"?
        self._save_name = save_name  # the name to save the output fastqs
        self._five_p_barcodes_pos, self._five_p_umi_poses = find_bc_and_umi_pos(
            five_p_bcs
        )
        self._five_p_bc_dict = make_5p_bc_dict(
            five_p_bcs, min_score_5_p, dont_build_reference
        )
        self._min_score_5_p = min_score_5_p  #
        self._linked_bcds = (
            linked_bcds  # which 3' barcodes each 5' bc is linked - a dictionary
        )
        self._three_p_bc_dict_of_dicts, self._min_score_d = make_dict_of_3p_bc_dicts(
            self._linked_bcds, three_p_mismatches, dont_build_reference
        )  # a dict of dicts - each 5' barcodes has
        # a dict of which 3' barcode matches the given sequence, which is also a dict
        self._ultra_mode = ultra_mode
        self._output_directory = output_directory
        self._i2 = i2  # paired end file name, or false
        self._adapter2 = adapter2
        self._min_trim = min_trim
        self._ignore_no_match = ignore_no_match
        self._barcodes_no_N = remove_Ns_from_barcodes(five_p_bcs)
        self._linked_bcds_no_N = {
            key: remove_Ns_from_barcodes(value) for (key, value) in linked_bcds.items()
        }
        self._keep_barcode = keep_barcode
        self._three_prime_only = three_prime_only
        self._TSO_umi_pos = TSO_umi_pos
        self._TSO_end = TSO_end

    def trim_and_cut(self, read, cutter, reads_quality_trimmed, reads_adaptor_trimmed):
        # /# first, trim by quality score
        q_start, q_stop = quality_trim_index(
            read.qualities, self._start_qc, self._end_qc
        )
        prev_length = len(read.sequence)
        read = read[q_start:q_stop]
        if not len(read.sequence) == prev_length:
            # then it was trimmed
            reads_quality_trimmed += 1

        # /# then, trim adapter
        prev_length = len(read.sequence)
        read = cutter(read, ModificationInfo(read))
        if not len(read.sequence) == prev_length:
            # then it was trimmed
            trimmed = True
            reads_adaptor_trimmed += 1
        else:
            trimmed = False

        if prev_length - len(read.sequence) >= self._min_trim:
            min_trimmed = True
        else:
            min_trimmed = False

        return read, trimmed, reads_quality_trimmed, reads_adaptor_trimmed, min_trimmed

    def run(self):
        logger = logging.getLogger("first_logger")
        while (
            True
        ):  # /# once spawned, this keeps running forever, until poison pill recieved
            if self._i2 is False:  # then it's single end
                # Notify reader that we need data
                self._need_work_queue.put(self._id)

                # /# get some data
                chunk_index = self._read_pipe.recv()

                # /# check there's no error
                if chunk_index == -1:  # /# poison pill from Sina
                    # reader is done
                    break
                elif chunk_index == -2:
                    # An exception has occurred in the reader
                    e, tb_str = self._read_pipe.recv()
                    raise e

                # /# otherwise if we have no error, run...
                # /# get some bytes
                data = self._read_pipe.recv_bytes()
                infiles = io.BytesIO(data)

                # Define the cutter
                adapter = [
                    BackAdapter(self._adapter, max_error_rate=0.1, min_overlap=1)
                ]
                cutter = AdapterCutter(adapter, times=1)

                # /# process the reads
                processed_reads = []
                this_buffer_dict = {}
                reads_written = 0
                assigned_reads = 0
                reads_quality_trimmed = 0
                reads_adaptor_trimmed = 0
                five_no_three_reads = 0
                five_no_three_reads_pass_length = 0
                assigned_reads_pass_length = 0

                for read in InputFiles(infiles).open():
                    reads_written += 1
                    umi = ""
                    (
                        read,
                        trimmed,
                        reads_quality_trimmed,
                        reads_adaptor_trimmed,
                        min_trimmed,
                    ) = self.trim_and_cut(
                        read, cutter, reads_quality_trimmed, reads_adaptor_trimmed
                    )

                    # /# demultiplex at the 5' end ###
                    read.name = (
                        read.name.replace(" ", "_").replace("/", "").replace("\\", "")
                    )  # remove bad characters

                    read, five_p_bc, umi, to_remove = five_p_demulti(
                        read,
                        self._five_p_barcodes_pos,
                        self._five_p_umi_poses,
                        self._five_p_bc_dict,
                        add_umi=True,
                        barcodes_no_N=self._barcodes_no_N,
                        min_score=self._min_score_5_p,
                        keep_barcode=self._keep_barcode,
                    )

                    # demultiplex at the 3' end
                    # First, check if this 5' barcode has any 3' barcodes
                    try:
                        linked = self._linked_bcds[five_p_bc]
                    except KeyError:
                        linked = "_none_"

                    if (
                        linked == "_none_"
                    ):  # no 3' barcodes linked to this 3' barcode, this includes "no_match" 5' barcdoes
                        comb_bc = "_5bc_" + five_p_bc

                        if len(read) >= self._final_min_length:
                            try:
                                this_buffer_dict[comb_bc].append(read)
                            except KeyError:
                                this_buffer_dict[comb_bc] = [read]

                        if not five_p_bc == "no_match":
                            assigned_reads += 1
                        if (
                            not five_p_bc == "no_match"
                            and len(read) >= self._final_min_length
                        ):
                            assigned_reads_pass_length += 1

                    elif (
                        min_trimmed
                    ):  # if it is linked to 3' barcodes and has been trimmed
                        read, three_p_bc, umi_3, to_remove_3 = three_p_demultiplex(
                            read,
                            self._three_p_bc_dict_of_dicts[five_p_bc],
                            add_umi=True,
                            linked_bcds=self._linked_bcds[five_p_bc],
                            linked_bcds_no_N=self._linked_bcds_no_N[five_p_bc],
                            min_score=self._min_score_d[five_p_bc],
                            keep_barcodes=self._keep_barcode,
                        )

                        if not three_p_bc == "no_match":
                            assigned_reads += 1
                        if (
                            not three_p_bc == "no_match"
                            and len(read) >= self._final_min_length
                        ):
                            assigned_reads_pass_length += 1

                        if three_p_bc == "no_match":
                            five_no_three_reads += 1
                        if (
                            three_p_bc == "no_match"
                            and len(read) >= self._final_min_length
                        ):
                            five_no_three_reads_pass_length += 1

                        comb_bc = "_5bc_" + five_p_bc + "_3bc_" + three_p_bc

                        if len(read) >= self._final_min_length:
                            try:
                                this_buffer_dict[comb_bc].append(read)
                            except KeyError:
                                this_buffer_dict[comb_bc] = [read]

                    else:  # linked to 3' barcode but wasn't trimmed so can't assign
                        comb_bc = "_5bc_" + five_p_bc + "_3bc_" + "no_match"
                        try:
                            this_buffer_dict[comb_bc].append(read)
                        except KeyError:
                            this_buffer_dict[comb_bc] = [read]

                ## Write out! ##
                read_number = defaultdict(int)
                for demulti_type, reads in this_buffer_dict.items():
                    read_number[demulti_type] = write_tmp_files(
                        output_dir=self._output_directory,
                        save_name=self._save_name,
                        demulti_type=demulti_type,
                        worker_id=self._id,
                        reads=reads,
                        ultra_mode=self._ultra_mode,
                        ignore_no_match=self._ignore_no_match,
                    )

            else:  # if paired end
                # Notify reader that we need data
                self._need_work_queue.put(self._id)

                # /# get some data
                chunk_index = self._read_pipe.recv()

                # /# check there's no error
                if chunk_index == -1:  # /# poison pill from Sina
                    # reader is done
                    break
                elif chunk_index == -2:
                    # An exception has occurred in the reader
                    e, tb_str = self._read_pipe.recv()
                    raise e

                # /# otherwise if we have no error, run...
                # /# get some bytes
                data1 = self._read_pipe.recv_bytes()
                infiles1 = io.BytesIO(data1)
                data2 = self._read_pipe.recv_bytes()
                infiles2 = io.BytesIO(data2)

                # Define the cutter
                adapter1 = [
                    BackAdapter(self._adapter, max_error_rate=0.1, min_overlap=1)
                ]
                adapter2 = [
                    BackAdapter(self._adapter2, max_error_rate=0.1, min_overlap=1)
                ]
                cutter1 = AdapterCutter(adapter1, times=1)
                cutter2 = AdapterCutter(adapter2, times=1)

                # /# process the reads
                processed_reads = []
                this_buffer_dict_1 = {}
                this_buffer_dict_2 = {}
                reads_written = 0
                assigned_reads = 0
                reads_quality_trimmed = 0
                reads_adaptor_trimmed = 0
                five_no_three_reads = 0
                five_no_three_reads_pass_length = 0
                assigned_reads_pass_length = 0

                for read1, read2 in zip(
                    InputFiles(infiles1).open(), InputFiles(infiles2).open()
                ):
                    # First, check that they have the same name
                    # assert read1.name == read2.name, "Paired reads are mismatched"

                    umi = ""
                    reads_written += 1
                    # /# first, trim by quality score
                    (
                        read1,
                        trimmed1,
                        reads_quality_trimmed,
                        reads_adaptor_trimmed,
                        min_trimmed1,
                    ) = self.trim_and_cut(
                        read1, cutter1, reads_quality_trimmed, reads_adaptor_trimmed
                    )

                    (
                        read2,
                        trimmed2,
                        ignore,
                        ignore_also,
                        min_trimmed2,
                    ) = self.trim_and_cut(
                        read2, cutter2, reads_quality_trimmed, reads_adaptor_trimmed
                    )

                    # /# demultiplex at the 5' end ###
                    read1.name = (
                        read1.name.replace(" ", "_").replace("/", "").replace("\\", "")
                    )  # remove bad characters
                    read2.name = (
                        read2.name.replace(" ", "_").replace("/", "").replace("\\", "")
                    )  # remove bad characters

                    # demultiplex based on 5' end
                    read1, five_p_bc, umi, to_remove_5 = five_p_demulti(
                        read1,
                        self._five_p_barcodes_pos,
                        self._five_p_umi_poses,
                        self._five_p_bc_dict,
                        add_umi=False,
                        barcodes_no_N=self._barcodes_no_N,
                        min_score=self._min_score_5_p,
                        keep_barcode=self._keep_barcode,
                    )

                    # add the 5' umi to each
                    for read in [read1, read2]:
                        if not "rbc:" in read.name:
                            read.name = read.name + "rbc:" + umi
                        else:  # if rbc is already there
                            read.name = read.name + umi

                    # demultiplex at the 3' end
                    # First, check if this 5' barcode has any 3' barcodes
                    try:
                        linked = self._linked_bcds[five_p_bc]
                    except:
                        linked = "_none_"

                    if (
                        linked == "_none_"
                    ):  # no 3' barcodes linked to this 5' barcode, this includes "no_match" 5'
                        # barcodes

                        if self._three_prime_only:
                            comb_bc = "_3bc_" + careful_rev_c(five_p_bc)

                            # Need to remove TSO seq and move UMI to header
                            # Have to use 'read2' because this is actually read 1, it's just we've flipped it

                            TSO_umi = read2.sequence[0:self._TSO_umi_pos]
                            read2 = read2[self._TSO_end:]

                            for read in [read1, read2]:
                                if not "rbc:" in read.name:
                                    read.name = read.name + "rbc:" + TSO_umi
                                else:  # if rbc is already there
                                    read.name = read.name + TSO_umi
                        else:
                            comb_bc = "_5bc_" + five_p_bc

                        # remove the 5' barcoded adapter from reverse read
                        read2 = remove_mate_adapter(
                            read=read2,
                            to_remove=to_remove_5,
                            bcd=five_p_bc,
                            trimmed=trimmed2,
                        )

                        if (
                            len(read1) >= self._final_min_length
                            and len(read2) >= self._final_min_length
                        ):
                            try:
                                this_buffer_dict_1[comb_bc].append(read1)
                                this_buffer_dict_2[comb_bc].append(read2)
                            except KeyError:
                                this_buffer_dict_1[comb_bc] = [read1]
                                this_buffer_dict_2[comb_bc] = [read2]

                        if not five_p_bc == "no_match":
                            assigned_reads += 1
                        if (
                            not five_p_bc == "no_match"
                            and len(read1) >= self._final_min_length
                            and len(read2) >= self._final_min_length
                        ):
                            assigned_reads_pass_length += 1

                    else:  # if there's a linked 3' barcode, no need to check if trimmed b/c P.E.
                        read2, three_p_bc, umi_3, to_remove_3 = three_p_demultiplex(
                            read2,
                            self._three_p_bc_dict_of_dicts[five_p_bc],
                            add_umi=False,
                            linked_bcds=self._linked_bcds[five_p_bc],
                            min_score=self._min_score_d[five_p_bc],
                            linked_bcds_no_N=self._linked_bcds_no_N[five_p_bc],
                            reverse_complement=True,
                            keep_barcodes=self._keep_barcode,
                        )

                        # add the second umi
                        for read in [read1, read2]:
                            if not "rbc:" in read.name:
                                read.name = read.name + "rbc:" + umi_3
                            else:  # if rbc is already there
                                read.name = read.name + umi_3

                        # remove the 5' bcded adapter from the reverse read
                        read2 = remove_mate_adapter(
                            read=read2,
                            to_remove=to_remove_5,
                            bcd=five_p_bc,
                            trimmed=trimmed2,
                        )
                        # remove the 3' adapter from the forward read
                        read1 = remove_mate_adapter(
                            read=read1,
                            to_remove=to_remove_3,
                            bcd=three_p_bc,
                            trimmed=trimmed1,
                        )

                        if not three_p_bc == "no_match":
                            assigned_reads += 1
                        if (
                            not three_p_bc == "no_match"
                            and len(read1) >= self._final_min_length
                            and len(read2) >= self._final_min_length
                        ):
                            assigned_reads_pass_length += 1

                        if three_p_bc == "no_match":
                            five_no_three_reads += 1
                        if (
                            three_p_bc == "no_match"
                            and len(read1) >= self._final_min_length
                            and len(read2) >= self._final_min_length
                        ):
                            five_no_three_reads_pass_length += 1

                        comb_bc = "_5bc_" + five_p_bc + "_3bc_" + three_p_bc

                        if (
                            len(read1) >= self._final_min_length
                            and len(read2) >= self._final_min_length
                        ):
                            try:
                                this_buffer_dict_1[comb_bc].append(read1)
                                this_buffer_dict_2[comb_bc].append(read2)
                            except KeyError:
                                this_buffer_dict_1[comb_bc] = [read1]
                                this_buffer_dict_2[comb_bc] = [read2]

                ## Write out! ##
                read_number = defaultdict(int)
                for demulti_type, reads in this_buffer_dict_1.items():

                    # If three_prime_only then we have to label the reads flipped, because the input read files
                    # were also flipped
                    if self._three_prime_only:
                        demulti_type = demulti_type + "_Rev"
                    else:
                        demulti_type = demulti_type + "_Fwd"

                    read_number[demulti_type] = write_tmp_files(
                        output_dir=self._output_directory,
                        save_name=self._save_name,
                        demulti_type=demulti_type,
                        worker_id=self._id,
                        reads=reads,
                        ultra_mode=self._ultra_mode,
                        ignore_no_match=self._ignore_no_match,
                    )

                for demulti_type, reads in this_buffer_dict_2.items():

                    # If three_prime_only then we have to label the reads flipped, because the input read files
                    # were also flipped
                    if self._three_prime_only:
                        demulti_type = demulti_type + "_Fwd"
                    else:
                        demulti_type = demulti_type + "_Rev"

                    read_number[demulti_type] = write_tmp_files(
                        output_dir=self._output_directory,
                        save_name=self._save_name,
                        demulti_type=demulti_type,
                        worker_id=self._id,
                        reads=reads,
                        ultra_mode=self._ultra_mode,
                        ignore_no_match=self._ignore_no_match,
                    )

            # LOG reads processed
            prev_total = self._total_demultiplexed.get()
            new_total = prev_total[0] + reads_written
            if new_total - prev_total[1] >= 1_000_000:
                prog_msg = str(new_total // 1_000_000) + " million reads processed"
                print(prog_msg)
                logger.info(prog_msg)
                last_printed = new_total
            else:
                last_printed = prev_total[1]
            self._total_demultiplexed.put([new_total, last_printed])

            # LOG quality trimming
            prev_total = self._total_reads_qtrimmed.get()
            new_total = prev_total + reads_quality_trimmed
            self._total_reads_qtrimmed.put(new_total)

            # LOG adaptor trimming
            prev_total = self._total_reads_adaptor_trimmed.get()
            new_total = prev_total + reads_adaptor_trimmed
            self._total_reads_adaptor_trimmed.put(new_total)

            # LOG reads assigned
            prev_total = self._total_reads_assigned.get()
            new_total = prev_total + assigned_reads
            self._total_reads_assigned.put(new_total)

            # LOG reads assigned per sample
            prev_dict = self._read_number.get()
            new_dict = merge_defaultdicts(prev_dict, read_number)
            self._read_number.put(new_dict)

            # LOG reads assigned that pass length filter
            prev_total = self._total_reads_assigned_pass_length_filter.get()
            new_total = prev_total + assigned_reads_pass_length
            self._total_reads_assigned_pass_length_filter.put(new_total)

            # LOG 5' no 3' reads
            prev_total = self._total_reads_5p_no_3p.get()
            new_total = prev_total + five_no_three_reads
            self._total_reads_5p_no_3p.put(new_total)

            # LOG 5' no 3' reads that pass length filter
            prev_total = self._total_reads_5p_no_3p_pass_length_filter.get()
            new_total = prev_total + five_no_three_reads_pass_length
            self._total_reads_5p_no_3p_pass_length_filter.put(new_total)


def remove_mate_adapter(read, to_remove, bcd, trimmed):
    """
    This function looks for the barcoded adapter (UMI-specific, not just Ns) in the mate read. It then checks
    that not too much was removed
    The idea is that if the forward read is B1B2B3B4B5AAAGGG..., where B is the barcode, then:
    to_remove = B1B2B3B4B5
    reverse read = CCCTTTb5b4b3b2b1iiiiii where i is the illumina adapter for the reverse read
    the reverse read then gets trimmed to CCCTTTb5b4b3b2b1
    So now we just look for bbbbb, which is the reverse complement of what was removed
    """

    if trimmed:  # then remove the barcode which should, by definition, be there
        read = read[0 : len(read.sequence) - len(bcd)]
    else:
        if len(to_remove) > 0:
            remove_rc = rev_c(to_remove)  # this is now "b5b4b3b2b1"
            adapter = [BackAdapter(remove_rc, max_error_rate=0.1, min_overlap=2)]
            cutter = AdapterCutter(adapter, times=1)
            read = cutter(read, ModificationInfo(read))

    return read


def write_tmp_files(
    output_dir, save_name, demulti_type, worker_id, reads, ultra_mode, ignore_no_match
):
    write_this = True  # assume true
    if "no_match" in demulti_type and ignore_no_match:
        write_this = False

    if ultra_mode and write_this:
        # /# work out this filename
        filename = (
            output_dir
            + "ultraplex_"
            + save_name
            + demulti_type
            + "_tmp_thread_"
            + str(worker_id)
            + ".fastq"
        )

        if os.path.exists(filename):
            append_write = "a"  # append if already exists
        else:
            append_write = "w"  # make a new file if not

        with open(filename, append_write) as file:
            this_out = []
            total_reads = 0
            for counter, read in enumerate(reads):
                # Quality control:
                assert len(read.name.split("rbc:")) <= 2, "Multiple UMIs in header!"
                total_reads += 1
                if counter == 0:
                    umi_l = len(read.name.split("rbc:")[1])
                assert (
                    len(read.name.split("rbc:")[1]) == umi_l
                ), "UMIs are different lengths"
                ## combine into a single list
                this_out.append("@" + read.name)
                this_out.append(read.sequence)
                this_out.append("+")
                this_out.append(read.qualities)

            output = "\n".join(this_out) + "\n"

            file.write(output)
            return total_reads

    elif write_this:
        # /# work out this filename
        filename = (
            output_dir
            + "ultraplex_"
            + save_name
            + demulti_type
            + "_tmp_thread_"
            + str(worker_id)
            + ".fastq.gz"
        )

        if os.path.exists(filename):
            append_write = "ab"  # append if already exists
        else:
            append_write = "wb"  # make a new file if not

        with gzip.open(filename, append_write) as file:
            this_out = []
            total_reads = 0
            for counter, read in enumerate(reads):
                # Quality control:
                assert len(read.name.split("rbc:")) <= 2, "Multiple UMIs in header!"
                total_reads += 1
                if counter == 0:
                    umi_l = len(read.name.split("rbc:")[1])
                assert (
                    len(read.name.split("rbc:")[1]) == umi_l
                ), "UMIs are different lengths"
                ## combine into a single list
                this_out.append("@" + read.name)
                this_out.append(read.sequence)
                this_out.append("+")
                this_out.append(read.qualities)

            output = "\n".join(this_out) + "\n"

            file.write(output.encode())
            return total_reads


def five_p_demulti(
    read,
    five_p_bc_pos,
    five_p_umi_poses,
    five_p_bc_dict,
    add_umi,
    keep_barcode,
    barcodes_no_N=[],
    min_score=0,
):
    """
    this function demultiplexes on the 5' end
    """
    sequence_length = len(read.sequence)
    this_five_p_umi = ""

    # check if "rbc:" already in header
    if "rbc:" in read.name:
        to_add = ""
    else:  # if not
        to_add = "rbc:"

    if sequence_length > max(five_p_bc_pos):
        # find best barcode match
        this_bc_seq = "".join([read.sequence[i] for i in five_p_bc_pos])

        if "dont_build" in five_p_bc_dict:
            winner = score_barcode_for_dict(
                this_bc_seq, barcodes_no_N, min_score, Ns_removed=True
            )
        else:
            winner = five_p_bc_dict[this_bc_seq]

        # store what sequence will be removed
        if sequence_length < len(winner):  # read is too short to contain barcode
            winner = "no_match"
            to_remove = ""
        else:
            to_remove = read.sequence[0 : len(winner)]

        if not winner == "no_match":
            # /# find umi sequence
            this_five_p_umi = "".join(
                [read.sequence[i] for i in five_p_umi_poses[winner]]
            )

            if not keep_barcode:
                # /# update read and read header
                read.sequence = read.sequence[len(winner) :]
                read.qualities = read.qualities[len(winner) :]

            if add_umi:
                # to read header add umi and 5' barcode info
                read.name = read.name.replace(" ", "_") + to_add + this_five_p_umi
        else:  # if no match
            read.name = read.name + to_add

    else:  # if read was too short to assign
        winner = "no_match"
        read.name = read.name + to_add
        to_remove = ""

    return read, winner, this_five_p_umi, to_remove


def find_bc_and_umi_pos(barcodes):
    """
    This function finds the coordinates of the umi and barcodes nucleotides
    """
    bcs_poses = {}
    umi_poses = {}
    for bc in barcodes:
        bcs_poses[bc] = [i for i in range(len(bc)) if bc[i] != "N"]
        umi_poses[bc] = [i for i in range(len(bc)) if bc[i] == "N"]

    # /# we assume that the barcode is always the same
    bc_pos = bcs_poses[barcodes[0]]

    umi_poses["no_match"] = []

    return bc_pos, umi_poses


def start_workers(
    n_workers,
    input_file,
    need_work_queue,
    adapter,
    five_p_bcs,
    three_p_bcs,
    save_name,
    total_demultiplexed,
    total_reads_assigned,
    total_reads_assigned_pass_length_filter,
    read_number,
    total_reads_qtrimmed,
    total_reads_adaptor_trimmed,
    total_reads_5p_no_3p,
    total_reads_5p_no_3p_pass_length_filter,
    min_score_5_p,
    three_p_mismatches,
    linked_bcds,
    three_p_trim_q,
    ultra_mode,
    output_directory,
    final_min_length,
    q5,
    i2,
    adapter2,
    min_trim,
    ignore_no_match,
    dont_build_reference,
    keep_barcode,
    three_prime_only,
    TSO_umi_pos,
    TSO_end
):
    """
    This function starts all the workers
    """
    workers = []
    all_conn_r = []
    all_conn_w = []

    total_demultiplexed.put(
        [0, 0]
    )  # [total written, last time it was printed] - initialise [0,0]
    total_reads_assigned.put(0)
    total_reads_assigned_pass_length_filter.put(0)
    total_reads_qtrimmed.put(0)
    total_reads_adaptor_trimmed.put(0)
    total_reads_5p_no_3p.put(0)
    total_reads_5p_no_3p_pass_length_filter.put(0)
    read_number.put(defaultdict(int))

    for index in range(n_workers):
        # create a pipe to send data to this worker
        conn_r, conn_w = Pipe(duplex=False)
        all_conn_r.append(conn_r)
        all_conn_w.append(conn_w)

        worker = WorkerProcess(
            index=index,
            read_pipe=conn_r,  # this is the "read_pipe"
            need_work_queue=need_work_queue,  # worker tells the reader it needs work
            adapter=adapter,  # the 3' adapter to trim
            output_directory=output_directory,
            five_p_bcs=five_p_bcs,
            three_p_bcs=three_p_bcs,
            save_name=save_name,
            total_demultiplexed=total_demultiplexed,
            total_reads_assigned=total_reads_assigned,
            total_reads_assigned_pass_length_filter=total_reads_assigned_pass_length_filter,
            read_number=read_number,
            total_reads_qtrimmed=total_reads_qtrimmed,
            total_reads_adaptor_trimmed=total_reads_adaptor_trimmed,
            total_reads_5p_no_3p=total_reads_5p_no_3p,
            total_reads_5p_no_3p_pass_length_filter=total_reads_5p_no_3p_pass_length_filter,
            ultra_mode=ultra_mode,
            min_score_5_p=min_score_5_p,
            three_p_mismatches=three_p_mismatches,
            linked_bcds=linked_bcds,
            three_p_trim_q=three_p_trim_q,
            q5=q5,
            i2=i2,
            adapter2=adapter2,
            min_trim=min_trim,
            ignore_no_match=ignore_no_match,
            final_min_length=final_min_length,
            dont_build_reference=dont_build_reference,
            keep_barcode=keep_barcode,
            three_prime_only=three_prime_only,
            TSO_umi_pos=TSO_umi_pos,
            TSO_end=TSO_end
        )

        worker.start()
        workers.append(worker)

    return workers, all_conn_r, all_conn_w


def concatenate_files(
    save_name,
    ultra_mode,
    sbatch_compression,
    output_directory,
    sample_names,
    compression_threads=8,
):
    """
    this function concatenates all the files produced by the
    different workers, then sends an sbatch command to compress
    them all to fastqs.
    """
    # First, file all the unique file names we have, ignoring threads
    all_names = glob.glob(output_directory + "ultraplex_" + save_name + "*")

    all_types = []  # ignoring threads
    for name in all_names:
        this_type = name.split("_tmp_thread_")[0]

        if this_type not in all_types:
            all_types.append(this_type)  # this_type contains directory if applicable

    # now concatenate them
    if ultra_mode:
        for this_type in all_types:
            # find all files with this barcode (or barcode combination)
            filenames = sorted(
                glob.glob(this_type + "*")
            )  # this type already has output directory
            # then concatenate
            command = ""
            for name in filenames:
                command = command + name + " "
            command = "cat " + command + " > " + this_type + ".fastq"
            os.system(command)

            for name in filenames:
                os.remove(name)

            print("Compressing with pigz...")
            c_thread_n = "-p" + str(compression_threads)
            if sbatch_compression:
                os.system(
                    "sbatch -J compression -c "
                    + str(compression_threads)
                    + ' --time 4:00:00 --wrap="pigz '
                    + c_thread_n
                    + " "
                    + this_type
                    + '.fastq"'
                )
            else:
                os.system("pigz " + c_thread_n + " " + this_type + ".fastq")

        # check if compression is complete
        if sbatch_compression:
            finished = False
            print("Compressing....")
            while not finished:
                # assume it's complete
                complete = True
                # now actually check
                for this_type in all_types:
                    filename = glob.glob(this_type + "*")

                    if ".gz" not in filename[0]:
                        complete = False

                if complete:
                    finished = True
                    print("Compression complete!")
                else:
                    time.sleep(1)
    else:  # if not ultra_mode
        for this_type in all_types:
            # find all files with this barcode (or barcode combination)
            filenames = sorted(glob.glob(this_type + "*"))
            # then concatenate
            command = ""
            for name in filenames:
                command = command + name + " "
            command = "cat " + command + " > " + this_type + ".fastq.gz"
            os.system(command)

            for name in filenames:
                os.remove(name)

    if len(sample_names) > 0:  # ie the user has specified at least one sample name
        for this_type in all_types:
            sample_name_list = []
            # check whether any of the keys in the dictionary match this type
            for key, value in sample_names.items():
                if key in this_type and "no_match" not in this_type:
                    # check that if key doesn't contain 3' barcode, then the type does either
                    if "3bc_" in this_type:
                        if "3bc_" in key:
                            sample_name_list.append(key)
                    else:  # if 3' bc is not in the type, the by definition if won't be in the key
                        sample_name_list.append(key)

            if len(sample_name_list) == 1:  # check that only one name matches
                if "Fwd" in this_type:
                    os.system(
                        "mv "
                        + this_type
                        + ".fastq.gz "
                        + output_directory
                        + "ultraplex_"
                        + save_name
                        + "_"
                        + sample_names[sample_name_list[0]]
                        + "_Fwd.fastq.gz"
                    )
                elif "Rev" in this_type:
                    os.system(
                        "mv "
                        + this_type
                        + ".fastq.gz "
                        + output_directory
                        + "ultraplex_"
                        + save_name
                        + "_"
                        + sample_names[sample_name_list[0]]
                        + "_Rev.fastq.gz"
                    )
                else:
                    os.system(
                        "mv "
                        + this_type
                        + ".fastq.gz "
                        + output_directory
                        + "ultraplex_"
                        + save_name
                        + "_"
                        + sample_names[sample_name_list[0]]
                        + ".fastq.gz"
                    )
            if len(sample_name_list) > 1:
                print(
                    "Warning - unable to rename file as filename-sample combination was not unique"
                )


def clean_files(output_directory, save_name):
    files = glob.glob(output_directory + "ultraplex_" + save_name + "*")
    for file in files:
        os.remove(file)


def process_bcs(csv, mismatch_5p, three_prime_only):
    five_p_bcs = []
    three_p_bcs = []
    linked = {}
    sample_names = {}

    counter_5 = 0  # all 5' barcodes must be consistent

    with open(csv, "r") as file:
        for row in file:
            counter_3 = 0  # 3' barcodes only have to be consistent for each 5' barcode
            counter_5 += 1
            # First, find if there's a comma
            line = row.rstrip().replace(" ", "")
            comma_split = line.split(",")

            if len(comma_split) > 1:
                if comma_split[1] == "":
                    just_5p = True
                else:
                    just_5p = False
            elif len(comma_split) == 1:
                just_5p = True
            else:
                just_5p = False

            if just_5p:
                # then there's no 3' barcode (or we only have 3' barcodes but we're pretending they are 5' barcodes)
                five_p_bcs.append(comma_split[0].split(":")[0].upper())
                assert comma_split[0].count(":") <= 1, "multiple colons in 5' barcode"

                if three_prime_only:
                    name_prefix = '3bc_'
                else:
                    name_prefix = '5bc_'

                if len(comma_split[0].split(":")) > 1:
                    if comma_split[0].split(":")[1] != "":
                        sample_names[
                            name_prefix + comma_split[0].split(":")[0].upper()
                        ] = comma_split[0].split(":")[1]

                if counter_5 == 1:
                    fivelength = len(
                        comma_split[0].rstrip().upper().replace("N", "").split(":")[0]
                    )
                else:
                    assert (
                        len(
                            comma_split[0]
                            .rstrip()
                            .upper()
                            .replace("N", "")
                            .split(":")[0]
                        )
                        == fivelength
                    ), "5' barcodes not consistent"
            else:
                # then we have 3 prime bcds as well
                five_p_bcs.append(comma_split[0].split(":")[0].upper())
                assert comma_split[0].count(":") <= 1, "multiple colons in 5' barcode"
                if comma_split[0].count(":") == 1:
                    assert (
                        comma_split[0].count(":")[1] == ""
                    ), "You cannot name the 5' sample if it's linked to 3' barcodes (please remove 5' barcode name)"

                if counter_5 == 1:
                    fivelength = len(
                        comma_split[0].rstrip().upper().replace("N", "").split(":")[0]
                    )
                else:
                    assert (
                        len(
                            comma_split[0]
                            .rstrip()
                            .upper()
                            .replace("N", "")
                            .split(":")[0]
                        )
                        == fivelength
                    ), "5' barcodes not consistent"

                # find the 3' barcodes
                three_ps = [b for a, b in enumerate(comma_split) if a > 0 and b != ""]
                linked[comma_split[0]] = [a.split(":")[0].upper() for a in three_ps]

                for bc_and_sample in three_ps:
                    assert (
                        bc_and_sample.count(":") <= 1
                    ), "multiple colons in 3' barcode"

                    bc = bc_and_sample.split(":")[0].upper()
                    if len(bc_and_sample.split(":")) > 1:
                        if not bc_and_sample.split(":")[1] == "":
                            three_p_sample_name = bc_and_sample.split(":")[1]
                            sample_names[
                                "5bc_" + comma_split[0].split(":")[0] + "_3bc_" + bc
                            ] = three_p_sample_name

                    counter_3 += 1
                    three_p_bcs.append(bc)

                    if counter_3 == 1:
                        threelength = len(bc.replace("N", ""))
                    else:
                        assert (
                            len(bc.replace("N", "")) == threelength
                        ), "3' barcodes not consistent"

    assert len(sample_names) == len(
        set(sample_names.values())
    ), "Sample names are not unique!"

    # remove duplicates
    five_p_bcs = list(dict.fromkeys(five_p_bcs))
    three_p_bcs = list(dict.fromkeys(three_p_bcs))

    match_5p = fivelength - mismatch_5p

    return five_p_bcs, three_p_bcs, linked, match_5p, sample_names


def print_header():
    print("")
    print(
        "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
    )
    print(
        "@@@@@   @@@@   .@@   @@@@@@@          @@         @@@@@@    (@@@@@        @@@    @@@@@@@         @@    @@@    @"
    )
    print(
        "@@@@   @@@@    @@   ,@@@@@@@@@    @@@@@    @@@   @@@@   (   @@@@   @@@   @@    @@@@@@@   @@@@@@@@@@   #   @@@@"
    )
    print(
        "@@@   &@@@    @@    @@@@@@@@@%   @@@@@         @@@@(   @    @@@         @@    @@@@@@@         @@@@@     @@@@@@"
    )
    print(
        "@@    @@@    @@    @@@@@@@@@@   @@@@@    @@    @@@          @@   .@@@@@@@*   @@@@@@@    @@@@@@@@.   @    @@@@@"
    )
    print(
        "@@@       @@@@.        @@@@@   @@@@@&   @@@.   &   &@@@@    @    @@@@@@@@         @          @    @@@@    @@@@"
    )
    print(
        "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
    )
    print("")


def check_enough_space(
    output_directory, input_file, ignore_space_warning, ultra_mode, i2
):
    # First, find the free space on the output directory
    if output_directory == "":
        output_directory = os.getcwd()
    total, used, free = shutil.disk_usage(output_directory)

    if ultra_mode:
        multiplier = 0.098
    else:
        multiplier = 0.98

    if not i2 == False:
        multiplier = multiplier / 2
    # Find the size of the input file
    input_file_size = Path(input_file).stat().st_size
    if ignore_space_warning:
        if not input_file_size < multiplier * free:
            print("WARNING! System may not have enough free space to demultiplex")
            print("(Warning has been ignored)")
    else:
        assert (
            input_file_size < free * multiplier
        ), "Not enough free space. To ignore this warning use option --ignore_space_warning"


def check_N_position(bcds, type):
    # checks that UMI positions result in consistent barcode
    if not len(bcds) == 0:
        for counter, bcd in enumerate(bcds):
            # find positions of non-N
            non_N = [a for a, b in enumerate(bcd) if b != "N"]

            if type == "5":
                # then look for first non N
                ref_pos = min(non_N)
            else:
                # look for last non_n
                ref_pos = len(bcd) - max(
                    non_N
                )  # not just max(non_N) because need to allow for different UMI length at 5' end of 3' bcd

            if counter == 0:
                correct_pos = ref_pos
            else:
                assert ref_pos == correct_pos, "UMI positions not consistent"


def main(buffer_size=int(4 * 1024**2)):  # 4 MB
    print_header()
    start = time.time()

    ## PARSE COMMAND LINE ARGUMENTS ##

    parser = argparse.ArgumentParser(
        description="Ultra-fast demultiplexing of fastq files."
    )
    optional = parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")
    # version
    parser.add_argument("-v", "--version", action="version", version=__version__)
    # input
    required.add_argument(
        "-i",
        "--inputfastq",
        type=str,
        required=True,
        help="fastq file to be demultiplexed",
    )
    # barcodes csv
    required.add_argument(
        "-b",
        "--barcodes",
        type=str,
        required=True,
        help="barcodes for demultiplexing in csv format",
    )
    # output directory
    optional.add_argument(
        "-d",
        "--directory",
        type=str,
        default="",
        nargs="?",
        help="optional output directory",
    )
    # 5' mismatches
    optional.add_argument(
        "-m5",
        "--fiveprimemismatches",
        type=int,
        default=0,
        nargs="?",
        help="number of mismatches allowed for 5prime barcode [DEFAULT 0]",
    )
    # 3' mismatches
    optional.add_argument(
        "-m3",
        "--threeprimemismatches",
        type=int,
        default=0,
        nargs="?",
        help="number of mismatches allowed for 3prime barcode [DEFAULT 0]",
    )
    # 3 prime barcodes only
    optional.add_argument('--three_prime_only',
                          action='store_true',
                          default=False,
                          help='Use this if only using 3 prime barcodes')
    optional.add_argument('--tso_seq',
                          default="",
                          help="Specify this if the read 1 has a TSO. Should be of format NNNNIII or NNNNII where the N "
                               "bases are moved to the UMI and the I bases are ignored")
    # minimum quality score
    optional.add_argument(
        "-q",
        "--phredquality",
        type=int,
        default=30,
        nargs="?",
        help="phred quality score for 3prime end trimming",
    )
    # threads
    optional.add_argument(
        "-t", "--threads", type=int, default=4, nargs="?", help="threads [DEFAULT 4]"
    )
    # adapter sequence
    optional.add_argument(
        "-a",
        "--adapter",
        type=str,
        default="AGATCGGAAGAGCACACGTCTGAA",
        nargs="?",
        help="sequencing adapter to trim [DEFAULT Illumina AGATCGGAAGAGCACACGTCTGAA p7]. Specify 'p3' or "
             "AGATCGGAAGAGCGGTTCAG to trim legacy p3 sequence",
    )
    # name of output file
    optional.add_argument(
        "-o",
        "--outputprefix",
        type=str,
        default="demux",
        nargs="?",
        help="prefix for output sequences [DEFAULT demux]",
    )
    # use sbatch compression in ultra mode
    optional.add_argument(
        "-sb",
        "--sbatchcompression",
        action="store_true",
        default=False,
        help="whether to compress output fastq using SLURM sbatch",
    )
    # ultra mode
    optional.add_argument(
        "-u",
        "--ultra",
        action="store_true",
        default=False,
        help="whether to use ultra mode, which is faster but makes very large temporary files",
    )
    # free space ignore warning
    optional.add_argument(
        "-ig",
        "--ignore_space_warning",
        action="store_true",
        default=False,
        help="whether to ignore warnings that there is not enough free space",
    )
    # minimum final length of read
    optional.add_argument(
        "-l",
        "--final_min_length",
        type=int,
        default=20,
        nargs="?",
        help="minimum length of the final outputted reads",
    )
    # start qc
    optional.add_argument(
        "-q5",
        "--phredquality_5_prime",
        type=int,
        default=0,
        nargs="?",
        help="quality trimming minimum score from 5' end - use with caution!",
    )
    # paired end file
    optional.add_argument(
        "-i2",
        "--input_2",
        type=str,
        default="",
        nargs="?",
        help="Optional second fastq.gz input for paired end data",
    )
    # adaptor for paired end read
    optional.add_argument(
        "-a2",
        "--adapter2",
        type=str,
        default="AGATCGGAAGAGCGTCGTG",
        nargs="?",
        help="sequencing adaptor to trim for the reverse read [Default AGATCGGAAGAGCGTCGTG]",
    )
    optional.add_argument(
        "-mt",
        "--min_trim",
        type=int,
        default=3,
        nargs="?",
        help=(
            "When using single end reads for 3' demultiplexing, this is the minimum adapter trimming amount for a 3'"
            "barcode to be detected. Default = 3"
        ),
    )
    optional.add_argument(
        "-inm",
        "--ignore_no_match",
        action="store_true",
        default=False,
        help="Do not write reads for which there is no match.",
    )
    optional.add_argument(
        "-dbr",
        "--dont_build_reference",
        default=False,
        action="store_true",
        help="Skip the reference building step - for long barcodes this will improve speed.",
    )
    optional.add_argument(
        "-kbc",
        "--keep_barcode",
        default=False,
        action="store_true",
        help="Do not remove barcodes/UMIs from the read (UMIs will still be moved to the "
        "read header)",
    )

    parser._action_groups.append(optional)
    args = parser.parse_args()

    if args.adapter == "p3" or args.adapter == "P3":
        print("p3 adapter specified, therefore setting adapter for read 1 to AGATCGGAAGAGCGGTTCAG")
        args.adapter = 'AGATCGGAAGAGCGGTTCAG'

    output_directory = args.directory

    if not output_directory == "":
        if not output_directory[len(output_directory) - 1] == "/":
            output_directory = output_directory + "/"

    # Make output directory
    if not output_directory == "":
        if not os.path.exists(output_directory):
            os.mkdir(output_directory)

    # It's logging time bbyyy
    logger = setup_logger(
        "first_logger",
        output_directory + "ultraplex_" + str(start) + ".log",
        level=logging.DEBUG,
    )

    print(args)
    logger.info(args)

    file_name = args.inputfastq
    barcodes_csv = args.barcodes
    mismatch_5p = args.fiveprimemismatches
    three_p_mismatches = args.threeprimemismatches
    three_p_trim_q = args.phredquality
    threads = args.threads
    adapter = args.adapter
    adapter2 = args.adapter2
    save_name = args.outputprefix
    sbatch_compression = args.sbatchcompression
    ultra_mode = args.ultra
    ignore_space_warning = args.ignore_space_warning
    final_min_length = args.final_min_length
    q5 = args.phredquality_5_prime
    i2 = args.input_2
    min_trim = args.min_trim
    ignore_no_match = args.ignore_no_match
    dont_build_reference = args.dont_build_reference
    keep_barcode = args.keep_barcode

    if args.tso_seq != "":
        if args.three_prime_only == False:
            assert 0 == 1, 'Can only specify TSO seq when using three prime only!'
        else:

            args.tso_seq = args.tso_seq.replace(" ", "")

            TSO_end = len(args.tso_seq)
            TSO_umi_pos = max([i for i, l in enumerate(args.tso_seq) if l == "N"]) + 1  # this is the last position
    else:
        TSO_end = 0
        TSO_umi_pos = 0

    if i2 == "":
        i2 = False

    if ultra_mode:
        print(
            "Warning - ultra mode selected. This will generate very large temporary files!"
        )

    if not ultra_mode:
        if sbatch_compression:
            print("sbatch_compression can only be used in conjunction with ultra mode")
            print("setting sbatch_compression to false")
            sbatch_compression = False

    # assert output_directory=="" or output_directory[len(output_directory)-1]=="/", "Error! Directory must end with '/'"

    check_enough_space(
        output_directory, file_name, ignore_space_warning, ultra_mode, i2
    )

    # process the barcodes csv
    five_p_bcs, three_p_bcs, linked_bcds, min_score_5_p, sample_names = process_bcs(
        barcodes_csv, mismatch_5p, args.three_prime_only
    )

    # remove files from previous runs
    clean_files(output_directory, save_name)

    # Switch things up if 3 prime only
    if args.three_prime_only:
        assert args.input_2 is not False, 'In --three_prime_only mode you must have paired end sequencing!'
        assert len(three_p_bcs) == 0, 'You appear to have given both 3 and 5 prime barcodes but are using --three_prime_only mode...?'

        # Pretend that read1 is read2 and vice versa. This enables read2 to be passed to the 5' demultiplexing
        file_name = args.input_2
        i2 = args.inputfastq

        # Also need to reverse complement the barcodes that are given. This is because the 3' barcodes are given
        # as they are in read1. However, demultiplexing will use read2, which means it will be the reverse
        # complement
        five_p_bcs = [careful_rev_c(a) for a in five_p_bcs]

        # Also need to switch the adapters
        adapter2 = args.adapter
        adapter = args.adapter2

    # Note that if --three_prime_only then they will have been reverse complemented by the code above, and should
    # now be consistent using the 5' barcode criteria
    check_N_position(
        five_p_bcs, "5"
    )  # check 3' later so that different 5' barcodes can have different types of 3' bcd

    # Make a queue to which workers that need work will add
    # a signal
    need_work_queue = Queue()
    total_demultiplexed = Queue()
    total_reads_assigned = Queue()
    total_reads_assigned_pass_length_filter = Queue()
    total_reads_qtrimmed = Queue()
    total_reads_adaptor_trimmed = Queue()
    total_reads_5p_no_3p = Queue()
    total_reads_5p_no_3p_pass_length_filter = Queue()
    read_number = Queue()

    # make a bunch of workers
    workers, all_conn_r, all_conn_w = start_workers(
        n_workers=threads,
        input_file=file_name,
        need_work_queue=need_work_queue,
        adapter=adapter,
        five_p_bcs=five_p_bcs,
        three_p_bcs=three_p_bcs,
        save_name=save_name,
        total_demultiplexed=total_demultiplexed,
        total_reads_assigned=total_reads_assigned,
        total_reads_assigned_pass_length_filter=total_reads_assigned_pass_length_filter,
        read_number=read_number,
        total_reads_qtrimmed=total_reads_qtrimmed,
        total_reads_adaptor_trimmed=total_reads_adaptor_trimmed,
        total_reads_5p_no_3p=total_reads_5p_no_3p,
        total_reads_5p_no_3p_pass_length_filter=total_reads_5p_no_3p_pass_length_filter,
        min_score_5_p=min_score_5_p,
        three_p_mismatches=three_p_mismatches,
        linked_bcds=linked_bcds,
        three_p_trim_q=three_p_trim_q,
        ultra_mode=ultra_mode,
        output_directory=output_directory,
        q5=q5,
        i2=i2,
        adapter2=adapter2,
        min_trim=min_trim,
        ignore_no_match=ignore_no_match,
        final_min_length=final_min_length,
        dont_build_reference=dont_build_reference,
        keep_barcode=keep_barcode,
        three_prime_only=args.three_prime_only,
        TSO_umi_pos=TSO_umi_pos,
        TSO_end=TSO_end
    )

    print("Demultiplexing...")
    reader_process = ReaderProcess(
        file_name, all_conn_w, need_work_queue, buffer_size, i2
    )
    reader_process.daemon = True
    reader_process.run()

    concatenate_files(
        save_name, ultra_mode, sbatch_compression, output_directory, sample_names
    )

    total_processed_reads = total_demultiplexed.get()[0]
    runtime_seconds = str((time.time() - start) // 1)
    finishing_msg = (
        "Demultiplexing complete! "
        + str(total_processed_reads)
        + " reads processed in "
        + runtime_seconds
        + " seconds"
    )
    print(finishing_msg)
    logger.info(finishing_msg)

    # More stats for logging
    total_qtrim = total_reads_qtrimmed.get()
    total_qtrim_percent = (total_qtrim / total_processed_reads) * 100
    total_adaptortrim = total_reads_adaptor_trimmed.get()
    total_adaptortrim_percent = (total_adaptortrim / total_processed_reads) * 100
    total_5p_no3 = total_reads_5p_no_3p.get()
    total_5p_no3_percent = (total_5p_no3 / total_processed_reads) * 100
    total_ass = total_reads_assigned.get()
    total_ass_percent = (total_ass / total_processed_reads) * 100
    total_ass_passlength = total_reads_assigned_pass_length_filter.get()
    total_ass_passlength_percent = (total_ass_passlength / total_processed_reads) * 100
    qmsg = (
        str(total_qtrim)
        + " ("
        + str(round_sig(total_qtrim_percent, 3))
        + "%) reads quality trimmed"
    )
    logger.info(qmsg)
    amsg = (
        str(total_adaptortrim)
        + " ("
        + str(round_sig(total_adaptortrim_percent, 3))
        + "%) reads adaptor trimmed"
    )
    logger.info(amsg)
    fivemsg = (
        str(total_5p_no3)
        + " ("
        + str(round_sig(total_5p_no3_percent, 3))
        + "%) reads with correct 5' bc but 3' bc not found"
    )
    logger.info(fivemsg)
    total_reads_5p_no_3p_pass_length_filter_log = (
        total_reads_5p_no_3p_pass_length_filter.get()
    )
    total_reads_5p_no_3p_pass_length_filter_percent = (
        total_reads_5p_no_3p_pass_length_filter_log / total_processed_reads
    ) * 100
    fivepassmsg = (
        str(total_reads_5p_no_3p_pass_length_filter_log)
        + " ("
        + str(round_sig(total_reads_5p_no_3p_pass_length_filter_percent, 3))
        + "%) reads with correct 5' bc but 3' bc not found that pass length filter"
    )
    logger.info(fivepassmsg)
    assmsg = (
        str(total_ass)
        + " ("
        + str(round_sig(total_ass_percent, 3))
        + "%) reads correctly assigned to samples"
    )
    logger.info(assmsg)
    asslengthmsg = (
        str(total_ass_passlength)
        + " ("
        + str(round_sig(total_ass_passlength_percent, 3))
        + "%) reads written to sample files (pass length filter)"
    )
    logger.info(asslengthmsg)
    logger.info("Total reads written per sample:")
    reads_samps_Dict = read_number.get()
    for k, v in reads_samps_Dict.items():
        logger.info("{:<8} {:<15}".format(k, v))


if __name__ == "__main__":
    main()
