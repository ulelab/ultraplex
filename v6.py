import dnaio
from xopen import xopen
import sys
import os
import traceback
import io
from multiprocessing import Process, Pipe, Queue
from typing import List, IO, Optional, BinaryIO, TextIO, Any, Tuple
from cutadapt.qualtrim import quality_trim_index
from cutadapt.modifiers import AdapterCutter, ModificationInfo
from cutadapt.adapters import BackAdapter
import pickle
import gzip
import glob
import time
from abc import ABC, abstractmethod
from dnaio import Sequence
import argparse


def make_5p_bc_dict(barcodes, min_score):
	"""
	this function generates a dictionary that matches each possible sequence
	from the read with the best 5' barcode
	"""
	first_bc = barcodes[0]
	seq_length = len(first_bc.replace("N",""))

	# check all the barcodes are the same length (ignoring UMIs)
	for bc in barcodes:
		assert len(bc.replace("N","")) == seq_length, "Your experimental barcodes are different lengths."

	seqs = make_all_seqs(seq_length)

	# trim sequences to desired length
	# create dictionary
	barcode_dictionary = {}

	for seq in seqs:
		barcode_dictionary[seq] = score_barcode_for_dict(seq, barcodes, min_score)

	return barcode_dictionary

def score_barcode_for_dict(seq, barcodes, min_score):
	"""
	this function scores a given sequence against all the barcodes
	it's used for the 5' barcode only
	"""

	# first, remove Ns from barcode list
	barcodes_no_N = []
	for i in range(len(barcodes)):
		this_bc = barcodes[i]
		barcodes_no_N.append(this_bc.replace("N", ""))

	# Now, score the sequence
	score_list = []

	for this_bc in barcodes_no_N:
		# find length of this barcode
		this_bc_l = len(this_bc)

		# score the barcode against the read, penalty for N in the read
		score = sum(a == b for a, b in zip(this_bc, seq))

		# append to score list
		score_list.append(score)

	# Find the best score
	best_score = max(score_list)

	if best_score < min_score:
		winner = "no_match"
	else:
		# check that there is only one barcode with the max score
		winner_indicies = [i for i in range(len(score_list)) if score_list[i] == best_score]

		if len(winner_indicies) > 1:
			winner = "no_match"
		else:  # if there is only one
			winner = barcodes[winner_indicies[0]]

	return (winner)

class ReaderProcess(Process):
	"""
	Read chunks of FASTA or FASTQ data (single-end or paired) and send to a worker.

	The reader repeatedly

	- reads a chunk from the file(s)
	- reads a worker index from the Queue
	- sends the chunk to connections[index]

	and finally sends the stop token -1 ("poison pills") to all connections.
	"""

	def __init__(self, file, connections, queue, buffer_size):
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
		#self.stdin_fd = stdin_fd

	def run(self):
		# if self.stdin_fd != -1:
		# 	sys.stdin.close()
		# 	sys.stdin = os.fdopen(self.stdin_fd) #/# not sure why I need this!

		try:
			with xopen(self.file, 'rb') as f:
				#/# every chunk in the fastq gets a unique chunk_index
				for chunk_index, chunk in enumerate(dnaio.read_chunks(f, self.buffer_size)):
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

	def send_to_worker(self, chunk_index, chunk):
		worker_index = self.queue.get()  # get a worker that needs work
		connection = self.connections[worker_index]  # find the connection to this worker
		connection.send(chunk_index)  # send the index of this chunk to the worker
		connection.send_bytes(chunk)  # /# send the actual data to this worker

class OutputFiles:
	"""
	The attributes are open file-like objects except when demultiplex is True. In that case,
	untrimmed, untrimmed2, out and out2 are file names or templates
	as required by the used demultiplexer ('{name}' etc.).
	Files may also be None.
	"""
	def __init__(
		self,
		out: Optional[BinaryIO] = None,
		out2: Optional[BinaryIO] = None,
		untrimmed: Optional[BinaryIO] = None,
		untrimmed2: Optional[BinaryIO] = None,
		too_short: Optional[BinaryIO] = None,
		too_short2: Optional[BinaryIO] = None,
		too_long: Optional[BinaryIO] = None,
		too_long2: Optional[BinaryIO] = None,
		info: Optional[BinaryIO] = None,
		rest: Optional[BinaryIO] = None,
		wildcard: Optional[BinaryIO] = None,
		demultiplex: bool = False,
		force_fasta: Optional[bool] = None,
	):
		self.out = out
		self.out2 = out2
		self.untrimmed = untrimmed
		self.untrimmed2 = untrimmed2
		self.too_short = too_short
		self.too_short2 = too_short2
		self.too_long = too_long
		self.too_long2 = too_long2
		self.info = info
		self.rest = rest
		self.wildcard = wildcard
		self.demultiplex = demultiplex
		self.force_fasta = force_fasta

	def __iter__(self):
		yield self.out
		yield self.out2
		yield self.untrimmed
		yield self.untrimmed2
		yield self.too_short
		yield self.too_short2
		yield self.too_long
		yield self.too_long2
		yield self.info
		yield self.rest
		yield self.wildcard

class InputFiles:
	"""
	this is from cutadapt - basically just creates a dnaio object
	"""
	def __init__(self, file1: BinaryIO, interleaved: bool = False):
		self.file1 = file1
		self.interleaved = interleaved

	def open(self):
		return dnaio.open(self.file1,
			interleaved=self.interleaved, mode="r")

def make_3p_bc_dict(bcs, min_score):
	"""
	Generates a dictionary that matches a given sequence with
	its best 3' barcode match, or "no_match"
	"""
	# First, find length of all barcodes
	len_d = {}
	for bc in bcs:
		try:
			len_d[len(bc)] +=1
		except:
			len_d[len(bc)] = 1

	if not len(len_d) == 1:
		print("Error, 3' barcodes are not all the same length")
		return -1
	else: #ie if theres only one length
		bcl = len(bcs[0]) # barcode length

		all_seqs = make_all_seqs(bcl)

		three_p_match_d = {}

		for seq in all_seqs:
			# assume no match 
			three_p_match_d[seq] = "no_match"

			# now check if there actually is a match
			correct_bcs = []
			for bc in bcs:
				# Find positions on non-N characters in bc
				nN = [pos for pos, char in enumerate(bc) if char != 'N']

				score = sum(a == b for a, b in zip(bc, seq))
				score += -sum(a == 'N' for a in seq) # penalise N

				if score >= min_score:
					correct_bcs.append(bc)

			if len(correct_bcs) == 1:
				three_p_match_d[seq] = correct_bcs[0]

	return three_p_match_d

def make_all_seqs(l):
	"""
	Makes all possible sequences, including Ns, of length l
	"""
	nts = ['A', "C", "G", "T", "N"]

	all_seqs = nts
	
	for i in range(l-1):
		new_seqs = []
		for seq in all_seqs:
			for nt in nts:
				new_seqs.append(seq + nt)
		all_seqs = new_seqs

	return(all_seqs)


def three_p_demultiplex(read, d, length):
	"""
	read is a dnaio read
	d is the relevant match dictionary
	length is the length of the 3' barcodes (it assumes they're all the same)
	"""
	bc = read.sequence[-length:]

	assigned = d[bc]

	if not assigned == "no_match":
		read.sequence = read.sequence[0:(len(read)-length)]
		read.qualities = read.qualities[0:(len(read.qualities)-length)]
		# add to umi
		umi_poses = [a == 'N' for a in assigned]
		umi = ''.join(bc[a] for a in umi_poses)
		read.name = read.name + umi

	return read, assigned


def make_dict_of_3p_bc_dicts(linked_bcs, min_score):
	"""
	this function makes a different dictionary for each 5' barcode
	it also checks that they're all the same length
	"""
	d = {}
	lengths = {}
	for five_p_bc, three_p_bcs in linked_bcs.items():
		if len(three_p_bcs) > 0: # ie this 5' barcode has 3' barcodes

			this_dict = make_3p_bc_dict(three_p_bcs, min_score)

			d[five_p_bc] = this_dict

			# find the length of each barcode
			for bc in three_p_bcs:
				three_p_length = len(bc)
				lengths[len(bc)] = 1

	# check that barcodes are all the same length, or there are none
	assert len(lengths) <= 1, "Your barcodes are different lengths, this is not allowed."
	return d, three_p_length

class WorkerProcess(Process): #/# have to have "Process" here to enable worker.start()
	"""
	The worker repeatedly reads chunks of data from the read_pipe, runs the pipeline on it
	and sends the processed chunks to the write_pipe.

	To notify the reader process that it wants data, processes then writes out.
	"""
	def __init__(self, index,
				 read_pipe, need_work_queue,
				 adapter, five_p_bcs, three_p_bcs, 
				 save_name, total_demultiplexed,
				 min_score_5_p, min_score_3_p,
				 linked_bcs,
				 start_qc =0, end_qc = 30,
				 min_length = 22):
		super().__init__()
		self._id = index # the worker id
		self._read_pipe = read_pipe # the pipe the reader reads data from
		self._need_work_queue = need_work_queue # worker adds its id to this queue when it needs work
		self._start_qc = start_qc # quality score to trim qc from 5' end
		self._end_qc = end_qc # quality score to trim qc from 3' end
		self._total_demultiplexed = total_demultiplexed # a queue which keeps track of the total number of reads processed
		self._adapter = adapter # the 3' adapter 
		self._min_length = min_length # the minimum length of a read after quality and adapter trimming to include. Remember
									  # that this needs to include the length of the barcodes and umis, so should be quite long eg 22 nt
		self._three_p_bcs = three_p_bcs # not sure we need this? A list of all the 3' barcodes. But unecessary because of "linked_bcs"?
		self._save_name = save_name # the name to save the output fastqs
		self._five_p_barcodes_pos, self._five_p_umi_poses = find_bc_and_umi_pos(five_p_bcs)
		self._five_p_bc_dict = make_5p_bc_dict(five_p_bcs, min_score_5_p)
		self._min_score_5_p = min_score_5_p # 
		self._min_score_3_p = min_score_3_p
		self._linked_bcs = linked_bcs # which 3' barcodes each 5' bc is linked - a dictionary
		self._three_p_bc_dict_of_dicts, self._3p_length = make_dict_of_3p_bc_dicts(self._linked_bcs, 
			min_score = self._min_score_3_p) # a dict of dicts - each 5' barcodes has
					# a dict of which 3' barcode matches the given sequence, which is also a dict
		self._ultra_mode = ultra_mode

	def run(self):
		#try:
		# stats = Statistics()
		while True:  # /# once spawned, this keeps running forever, until poison pill recieved
			# Notify reader that we need data
			self._need_work_queue.put(self._id)

			#/# get some data
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
			#/# get some bytes
			data = self._read_pipe.recv_bytes()
			infiles = io.BytesIO(data)
			
			# Define the cutter
			adapter = [BackAdapter(self._adapter, max_error_rate=0.1)]
			cutter = AdapterCutter(adapter, times=3)

			#/# process the reads
			processed_reads = []
			five_p_bcs = []
			this_buffer_dict = {}
			reads_written = 0

			for read in InputFiles(infiles).open():
				reads_written += 1
				#/# first, trim by quality score
				q_start, q_stop = quality_trim_index(read.qualities, self._start_qc, self._end_qc)
				read = read[q_start:q_stop]
				#/# then, trim adapter
				prev_length = len(read.sequence)
				read = cutter(read, ModificationInfo(read))
				if not len(read.sequence) == prev_length:
					# then it was trimmed
					trimmed = True
				else:
					trimmed = False

				if len(read.sequence) > self._min_length:
					#/# demultiplex at the 5' end ###
					read.name = read.name.replace(" ", "").replace("/", "").replace("\\", "") # remove bad characters

					read, five_p_bc = five_p_demulti(read, five_p_bcs, 
						self._five_p_barcodes_pos,
						self._five_p_umi_poses,
						self._five_p_bc_dict)

					#/# demultiplex at the 3' end
					# First, check if this 5' barcode has any 3' barcodes
					try:
						linked = self._linked_bcs[five_p_bc]
					except:
						linked = "_none_"

					if linked == "_none_": # no 3' barcodes linked to this 3' barcode
						comb_bc = '_5bc_' + five_p_bc

						try:
							this_buffer_dict[comb_bc].append(read)
						except:
							this_buffer_dict[comb_bc] = [read]

					elif trimmed: # if it is linked to 3' barcodes and has been trimmed

						read, three_p_bc = three_p_demultiplex(read, 
							self._three_p_bc_dict_of_dicts[five_p_bc], 
							length = self._3p_length)

						comb_bc = '_5bc_' + five_p_bc + '_3bc_' + three_p_bc

						try:
							this_buffer_dict[comb_bc].append(read)
						except:
							this_buffer_dict[comb_bc] = [read]

			## Write out! ##
			for demulti_type, reads in this_buffer_dict.items():
				if self._ultra_mode:
					#/# work out this filename
					filename = 'mDe_'+self._save_name+demulti_type+'_tmp_thread_'+str(self._id)+'.fastq'

					if os.path.exists(filename):
						append_write = 'a' # append if already exists
					else:
						append_write = 'w' # make a new file if not

					with open(filename, append_write) as file:
						this_out = []
						for counter, read in enumerate(reads):
							
							# Quality control:
							assert len(read.name.split("rbc:")) <= 2, "Multiple UMIs in header!"
							
							if counter == 0:
								umi_l = len(read.name.split("rbc:")[1])
							assert len(read.name.split("rbc:")[1]) == umi_l, "UMIs are different lengths"
							## combine into a single list
							this_out.append("@" + read.name)
							this_out.append(read.sequence)
							this_out.append("+")
							this_out.append(read.qualities)

						output = '\n'.join(this_out) + '\n'
						#print(output)
						file.write(output)
				else:
					#/# work out this filename
					filename = 'mDe_'+self._save_name+demulti_type+'_tmp_thread_'+str(self._id)+'.fastq.gz'

					if os.path.exists(filename):
						append_write = 'ab' # append if already exists
					else:
						append_write = 'wb' # make a new file if not

					with gzip.open(filename, append_write) as file:
						this_out = []
						for counter, read in enumerate(reads):
							
							# Quality control:
							assert len(read.name.split("rbc:")) <= 2, "Multiple UMIs in header!"
							
							if counter == 0:
								umi_l = len(read.name.split("rbc:")[1])
							assert len(read.name.split("rbc:")[1]) == umi_l, "UMIs are different lengths"
							## combine into a single list
							this_out.append("@" + read.name)
							this_out.append(read.sequence)
							this_out.append("+")
							this_out.append(read.qualities)

						output = '\n'.join(this_out) + '\n'
						#print(output)
						file.write(output.encode())

			prev_total = self._total_demultiplexed.get()
			new_total = prev_total[0] + reads_written
			if new_total - prev_total[1] >= 1_000_000:
				print(str(new_total//1_000_000) + ' million reads processed')
				last_printed = new_total
			else:
				last_printed = prev_total[1]
			self._total_demultiplexed.put([new_total, last_printed])


def five_p_demulti(read, five_p_bcs, five_p_bc_pos, five_p_umi_poses,
	five_p_bc_dict):
	"""
	this function demultiplexes on the 5' end
	"""

	#/# find 5' umi
	this_bc_seq = ''.join([read.sequence[i] for i in five_p_bc_pos])
	winner = five_p_bc_dict[this_bc_seq]

	if not winner == "no_match":
		#/# find umi sequence
		this_five_p_umi = ''.join([read.sequence[i] for i in five_p_umi_poses[winner]])

		#/# update read and read header
		read.sequence = read.sequence[len(winner):]
		read.qualities= read.qualities[len(winner):]

		# to read header add umi and 5' barcode info
		read.name = (read.name.replace(" ", "")+"rbc:" + this_five_p_umi)
	else:
		read.name = read.name + 'rbc:'
	


	return read, winner

def find_bc_and_umi_pos(barcodes):
	"""
	This function finds the coordinates of the umi and barcodes nucleotides
	"""
	bcs_poses = {}
	umi_poses = {}
	for bc in barcodes:
		bcs_poses[bc] = [i for i in range(len(barcodes[0])) if barcodes[0][i] != "N"]
		umi_poses[bc] = [i for i in range(len(barcodes[0])) if barcodes[0][i] == "N"]

	#/# we assume that the barcode is always the same
	bc_pos = bcs_poses[barcodes[0]]

	umi_poses["no_match"] = []

	return bc_pos, umi_poses

def start_workers(n_workers, input_file, need_work_queue, adapter,
	five_p_bcs, three_p_bcs,  save_name, total_demultiplexed,
	min_score_5_p, min_score_3_p, linked_bcs, three_p_trim_q,
	ultra_mode):
	"""
	This function starts all the workers
	"""
	workers = []
	all_conn_r = []
	all_conn_w = []

	total_demultiplexed.put([0,0]) # [total written, last time it was printed] - initialise [0,0]

	for index in range(n_workers):
		
		# create a pipe to send data to this worker
		conn_r, conn_w = Pipe(duplex=False)
		all_conn_r.append(conn_r)
		all_conn_w.append(conn_w)

		worker = WorkerProcess(index,
			conn_r, # this is the "read_pipe"
			need_work_queue, # worker tells the reader it needs work
			adapter, # the 3' adapter to trim
			five_p_bcs, 
			three_p_bcs,
			save_name, 
			total_demultiplexed, 
			min_score_5_p,
			min_score_3_p = min_score_3_p,
			linked_bcs = linked_bcs,
			end_qc = three_p_trim_q)

		worker.start()
		workers.append(worker)

	return workers, all_conn_r, all_conn_w

def concatenate_files(save_name, sbatch_compression,
	compression_threads = 4, 
	ultra_mode):
	""" 
	this function concatenates all the files produced by the 
	different workers, then sends an sbatch command to compress
	them all to fastqs
	"""

	# First, file all the unique file names we have, ignoring threads
	all_names = glob.glob("mDe_" + save_name +'*')

	all_types = [] # ignoring threads
	for name in all_names:
		this_type = name.split("_tmp_thread_")[0]
		if this_type not in all_types:
			all_types.append(this_type)

	# now concatenate them
	if ultra_mode:
		for this_type in all_types:
			# find all files with this barcode (or barcode combination)
			filenames = glob.glob(this_type + '*')
			# then concatenate
			command = ''
			for name in filenames:
				command = command + name + ' '
				command = 'cat ' + command + ' > ' + this_type + '.fastq'
				os.system(command)
				print("Compressing with pigz...")
				c_thread_n = '-p' + str(compression_threads)
				if sbatch_compression: 
					os.system('sbatch -J compression --time 4:00:00 --wrap="pigz '+c_thread_n+' '+this_type+'.fastq"')
				else:
					os.system('pigz '+c_thread_n+' '+this_type+'.fastq')

				for name in filenames:
					os.remove(name)

			# check if compression is complete
			finished = False
			print("Compressing....")
			while not finished:
				# assume it's complete
				complete = True
				# now actually check 
				for this_type in all_types:
					filename = glob.glob(this_type + '*')
					
					if '.gz' not in filename[0]:
						complete = False

				if complete:
					finished = True
					print("Compression complete!")
				else:
					time.sleep(1)
	else: # if not ultra_mode
		for this_type in all_types:
			# find all files with this barcode (or barcode combination)
			filenames = glob.glob(this_type + '*')
			# then concatenate
			command = ''
			for name in filenames:
				command = command + name + ' '
				command = 'cat ' + command + ' > ' + this_type + '.fastq.gz'
				os.system(command)

def clean_files(save_name):
	files = glob.glob('mDe_' + save_name +'*')
	for file in files:
		os.remove(file)

def process_bcs(csv, mismatch_5p, mismatch_3p):

	five_p_bcs = []
	three_p_bcs = []
	linked = {}

	with open(csv, 'r') as file:
		for row in file:
			# First, find if theres a comma
			line = row.rstrip()
			comma_split = line.split(',')
			
			if len(comma_split) == 1:
				# then there's no 3' barcode
				five_p_bcs.append(comma_split[0])
				fivelength=len(comma_split[0].replace("N",""))
			elif len(comma_split) == 2:
				# then we have 3 prime bcds
				five_p_bcs.append(comma_split[0])
				# find the 3' barcodes
				three_ps = comma_split[1].split(";")
				linked[comma_split[0]] = three_ps

				for bc in three_ps:
					three_p_bcs.append(bc)
					threelength=len(bc.replace("N",""))

	# remove duplicates
	five_p_bcs = list(dict.fromkeys(five_p_bcs))
	three_p_bcs = list(dict.fromkeys(three_p_bcs))

	match_5p = fivelength - mismatch_5p
	match_3p = threelength - mismatch_3p
	return five_p_bcs, three_p_bcs, linked, match_5p, match_3p

def print_header():
	print("")
	print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
	print("@@@@@   @@@@   .@@   @@@@@@@          @@         @@@@@@    (@@@@@        @@@    @@@@@@@         @@    @@@    @")
	print("@@@@   @@@@    @@   ,@@@@@@@@@    @@@@@    @@@   @@@@   (   @@@@   @@@   @@    @@@@@@@   @@@@@@@@@@   #   @@@@")
	print("@@@   &@@@    @@    @@@@@@@@@%   @@@@@         @@@@(   @    @@@         @@    @@@@@@@         @@@@@     @@@@@@")
	print("@@    @@@    @@    @@@@@@@@@@   @@@@@    @@    @@@          @@   .@@@@@@@*   @@@@@@@    @@@@@@@@.   @    @@@@@")
	print("@@@       @@@@.        @@@@@   @@@@@&   @@@.   &   &@@@@    @    @@@@@@@@         @          @    @@@@    @@@@")
	print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
	print("")
    
def main(buffer_size = int(4*1024**2)): # 4 MB
	print_header()
	start = time.time()

	## PARSE COMMAND LINE ARGUMENTS ##

	parser = argparse.ArgumentParser(description='Ultra-fast demultiplexing of fastq files.')
	optional = parser._action_groups.pop()
	required = parser.add_argument_group('required arguments')
	required.add_argument('-i',"--inputfastq", type=str, required=True,
						help='fastq file to be demultiplexed')
	required.add_argument('-b',"--barcodes", type=str, required=True,
						help='barcodes for demultiplexing in csv format')
	optional.add_argument('-m5',"--fiveprimemismatches", type=int, default=1, nargs='?',
						help='number of mismatches allowed for 5prime barcode [DEFAULT 1]')
	optional.add_argument('-m3',"--threeprimemismatches", type=int, default=0, nargs='?',
						help='number of mismatches allowed for 3prime barcode [DEFAULT 0]')
	optional.add_argument('-q',"--phredquality", type=int, default=30, nargs='?',
						help='phred quality score for 3prime end trimming')
	optional.add_argument('-t',"--threads", type=int, default=8, nargs='?',
						help='threads [DEFAULT 8]')
	optional.add_argument('-a',"--adapter", type=str, default="AGATCGGAAGAGCGGTTCAG", nargs='?',
						help='sequencing adapter to trim [DEFAULT Illumina AGATCGGAAGAGCGGTTCAG]')
	optional.add_argument('-o',"--outputprefix", type=str, default="demux", nargs='?',
						help='prefix for output sequences [DEFAULT demux]')
	optional.add_argument('-sb',"--sbatchcompression", action='store_true', default=False,
						help='whether to compress output fastq using SLURM sbatch')
	optional.add_argument('-u',"--ultra", action='store_true', default=False,
					help='whether to use ultra mode, which is faster but makes very large temporary files')
	parser._action_groups.append(optional)
	args = parser.parse_args()
	print(args)
	file_name = args.inputfastq
	barcodes_csv = args.barcodes
	mismatch_5p = args.fiveprimemismatches
	mismatch_3p = args.threeprimemismatches
	three_p_trim_q = args.phredquality
	threads = args.threads
	adapter = args.adapter
	save_name = args.outputprefix
	sbatch_compression = args.sbatchcompression
	ultra_mode = args.ultra

	if ultra_mode:
		print("Warning - ultra mode selected. This will generate very large temporary files!")

	# process the barcodes csv 
	five_p_bcs, three_p_bcs, linked_bcs, min_score_5_p, min_score_3_p = process_bcs(barcodes_csv, mismatch_5p, mismatch_3p)
	
	# remove files from previous runs
	clean_files(save_name)

	#/# Make a queue to which workers that need work will add
	#/# a signal
	need_work_queue = Queue()
	total_demultiplexed=Queue()

	#/# make a bunch of workers
	workers, all_conn_r, all_conn_w = start_workers(threads, 
		file_name, need_work_queue,
		adapter, five_p_bcs,
		three_p_bcs, save_name, 
		total_demultiplexed,
		min_score_5_p, min_score_3_p, 
		linked_bcs, three_p_trim_q,
		ultra_mode)

	print("Demultiplexing...")
	reader_process = ReaderProcess(file_name, all_conn_w,
			need_work_queue, buffer_size)
	reader_process.daemon = True
	reader_process.run()

	concatenate_files(save_name, sbatch_compression)
	print("Demultiplexing complete! " + str(total_demultiplexed.get()[0])+' reads written in ' +str((time.time()-start)//1) + ' seconds')


if __name__ == "__main__":
	main()
