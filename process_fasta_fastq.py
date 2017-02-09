#!/user/bin/env python

import sys, os
import argparse
from string import maketrans

from Bio import SeqIO

Description="Program to process fasta or fastq file in wide range of aspects Check the help for different types of available options."
usage="""
python {script} --input fasta/fastqFile --read
python {script} --input fasta/fastqFile --stats

""".format(script=sys.argv[0])

parser=argparse.ArgumentParser(description=Description, epilog=usage)
parser.add_argument("-i", "--input", action="store", dest="input", help="Fasta or fastq input file")
parser.add_argument("--stats", action="store_true", default=False, help="Prints the basic statistics of the reads")
parser.add_argument("--interval", action="store", dest="interval", help="gets the sequences in interval specified. Provide start point and then interval. e.g. 3,2")
parser.add_argument("--seqid", action="store", dest="seqid", help="comma separated list of sequence id from input file")
parser.add_argument("-l", "--getlength", action="store_true", dest="getlength", default=False, help="Outputs the sequence ID and sequence length (in tab-delimited)")
parser.add_argument("--filterbylength", action="store_true", dest="filterbylength", help="filter reads by length")
parser.add_argument("--subseq", action="store_true", dest="subseq", default=False, help="get subsequence from sequence reads. Default: gets first 100 bps in every sequence")
parser.add_argument("-x", "--min", action="store", dest="min", type=int, default=1, help="provide minimum position for subsequence. Must provide --subseq")
parser.add_argument("-y", "--max", action="store", dest="max", type=int, help="provide maximum position for subsequence. Must provide --subseq")
parser.add_argument("--leftclip", action="store", dest="leftclip", type=int, help="removes [int] bases from left end or 5' end")
parser.add_argument("--rightclip", action="store", dest="rightclip", type=int, help="removes [int] bases from right end or 3' end")
parser.add_argument("-r", "--reverse", action="store_true", dest="reverse", default=False, help="reverses the sequence, not reverse complement")
parser.add_argument("--reversecomp", action="store_true", dest="reversecomp", default=False, help="reverse complements sequence reads")
parser.add_argument("--random", action="store_true", dest="random", default=False, help="Generate 50 percent of input sequences randomly using seed value, meaning you get same set of random sequences every time.")
parser.add_argument("-o", "--output", action="store", dest="output", help="output filename")




options=parser.parse_args()
options.temp_input = options.input

if (options.leftclip and (options.subseq or options.max)):
	print "You cannot provide --leftclip option along with --subseq or --max/--min"
	exit(1)
if (options.rightclip and (options.subseq or options.max)):
	print "You cannot provide --rightclip option along with --subseq or --max/--min"
	exit(1)


def convert_int_to_ascii_char(scores):
	ascii_val=""

	for each_score in scores:
		ascii_val+=chr(int(each_score)+33)

	return ascii_val



def reverse_sequence():
	#reverses the sequence reads

	fh=open(options.input)
	output="reverse_output" + options.inputfiletype
	out=open(output, "w")

	for record in SeqIO.parse(fh, options.inputfiletype):
		out.write(options.symbol + record.description + "\n" + str(record.seq)[::-1] + "\n")
		if options.inputfiletype == "fastq":
			out.write("+" + "\n" + convert_int_to_ascii_char(record.letter_annotations["phred_quality"][::-1]) + "\n")

	fh.close()
	out.close()
	options.input=output

	return

def reverse_complement_bases(ntseq):
	#reverse the ntseq
	ntseq=ntseq[::-1]
	from_this="actgnACTGN"
	to_this="tgacnTGACN"

	#complement the ntseq and return
	return ntseq.translate(maketrans(from_this, to_this))

def reverse_complement():
		#reverse complements the sequence reads

		fh=open(options.input)
		output="reverse_output" + options.inputfiletype
		out=open(output, "w")

		for record in SeqIO.parse(fh, options.inputfiletype):
			out.write(options.symbol + record.description + "\n" + reverse_complement_bases(str(record.seq)) + "\n")
			if options.inputfiletype == "fastq":
				out.write("+" + "\n" + convert_int_to_ascii_char(record.letter_annotations["phred_quality"][::-1]) + "\n")		#quality scores cannot be reverse complemented, so just reverse them

		fh.close()
		out.close()
		options.input=output

		return


def doleftclip():
	fh=open(options.input)
	output = "clipping." + options.inputfiletype
	out = open(output, "w")

	for record in SeqIO.parse(fh, options.inputfiletype):
		out.write(options.symbol + record.description + "\n" + str(record.seq)[options.leftclip:] + "\n")
		if options.inputfiletype == "fastq":
			out.write("+" + "\n" + convert_int_to_ascii_char(record.letter_annotations["phred_quality"])[options.leftclip:] + "\n")

	fh.close()
	out.close()
	options.input=output
	return

def dorightclip():
	fh=open(options.input)
	output = "clipping." + options.inputfiletype
	out = open(output, "w")

	for record in SeqIO.parse(fh, options.inputfiletype):
		out.write(options.symbol + record.description + "\n" + str(record.seq)[:-options.rightclip] + "\n")
		if options.inputfiletype == "fastq":
			out.write("+" + "\n" + convert_int_to_ascii_char(record.letter_annotations["phred_quality"])[:-options.rightclip] + "\n")

	fh.close()
	out.close()
	options.input=output

	return

def doleftrightclip():
	fh=open(options.input)
	output = "clipping." + options.inputfiletype
	out = open(output, "w")

	for record in SeqIO.parse(fh, options.inputfiletype):
		out.write(options.symbol + record.description + "\n" + str(record.seq)[options.leftclip:-options.rightclip] + "\n")
		if options.inputfiletype == "fastq":
			out.write("+" + "\n" + convert_int_to_ascii_char(record.letter_annotations["phred_quality"])[options.leftclip:-options.rightclip] + "\n")

	fh.close()
	out.close()
	options.input=output

	return




def get_std_dev(scores):
	import math
	total=0
	variance=0
	std_dev=0
	mean=0
	for score in scores:
		total+=score
	mean=total/len(scores)
	#print mean
	total=0
	for score in scores:

		total+=(score-mean)**2
	variance=total/len(scores)
	#print variance
	return(math.sqrt(variance))

def low_score_bases(scores):
	count=0
	for eachscore in scores:
		if eachscore < 20:
			count+=1
	return count

def get_N50(upper_length_list):
	sum_of_all_lengths=sum(upper_length_list)
	total=0
	for contiglen in upper_length_list:

		total+=contiglen
		if float((total * 100)/sum_of_all_lengths) >=50.0:
			return contiglen

	return None

def get_N80(upper_length_list):
	sum_of_all_lengths=sum(upper_length_list)
        total=0
        for contiglen in upper_length_list:
                total+=contiglen
                if float((total * 100)/sum_of_all_lengths) >=80.0:
                        return contiglen

        return None

def get_N90(upper_length_list):
	sum_of_all_lengths=sum(upper_length_list)
        total=0
        for contiglen in upper_length_list:
                total+=contiglen
                if float((total * 100)/sum_of_all_lengths) >=90.0:
                        return contiglen

        return None


def get_stats():
	total_low_score_bases=0
	total_bases=0
	fh=open(options.input)
	upper_length_list=[]
	for record in SeqIO.parse(fh, options.inputfiletype):
		upper_length_list.append(len(str(record.seq)))
		total_bases+=len(str(record.seq))
		if options.inputfiletype=="fastq":
			total_low_score_bases+=low_score_bases(record.letter_annotations["phred_quality"])

	upper_length_list=sorted(upper_length_list)
	print "\nSummary of Statisitcs:\n"
	print "Total number of reads/contigs :", len(upper_length_list)
	print "Minimum read/contig length :", upper_length_list[0]
	print "Maximum read/contig length : ", upper_length_list[-1]
	print "Mean length of total reads/contigs : ", sum(upper_length_list)/len(upper_length_list)
	print "Standard Deviation of read length :", get_std_dev(upper_length_list)
	print "N50 read/contig length :", get_N50(upper_length_list)
	print "N80 read/contig length :", get_N80(upper_length_list)
	print "N90 read/contig length :", get_N90(upper_length_list)
	print "Total bases :", total_bases
	if options.inputfiletype=='fastq':
		print "Total bases with score less than 20 : ", total_low_score_bases, " which is" ,str(total_low_score_bases*100/total_bases)+"% of total bases"


def check_seqid_in_the_list(record):

	for seqid in options.idlist:

		if seqid.lower() == replace_special_chars(record.description.lower()) or seqid.lower() +" " in replace_special_chars(record.description.lower()):
			return True

	return None



def read_fasta_fastq():
	fh=open(options.input)
	for record in SeqIO.parse(fh, options.inputfiletype):
			print options.symbol + record.description
			print str(record.seq)
			if options.inputfiletype == "fastq":
				print "+"
				print convert_int_to_ascii_char(record.letter_annotations["phred_quality"])
	return





def replace_special_chars(sequence_id):
	special_chars="|!$^*(){}[]\/"
	for char in special_chars:
		sequence_id=sequence_id.replace(char, " ")
	return sequence_id


def get_seqid_list():
	options.idlist=[]
	if os.path.exists(options.seqid):
		# filename provided
		#open the file
		fh=open(options.seqid)
		for line in fh:
			line=line.strip()
			options.idlist.append(line)
		fh.close()
	else:
		options.idlist=options.seqid.split(",")


	return options.idlist

def get_seqs_by_id():
	get_seqid_list()
	fh=open(options.input)
	output = "seq_by_id_output." + options.inputfiletype
	out=open(output, "w")
	for record in SeqIO.parse(fh, options.inputfiletype):

		for seqid in options.idlist:
			tmpseqid=seqid
			if seqid.startswith("@") or seqid.startswith(">"):
				seqid=seqid[1:]
			#print seqid, record.id
			if seqid == record.id or seqid.lower() == replace_special_chars(record.description.lower()) or seqid.lower() +" " in replace_special_chars(record.description.lower()):
				#print "printing record"
				#get_reads_with_seqid(record, "fasta")
				if options.inputfiletype=="fastq":
					out.write("@"+record.description + "\n" + str(record.seq) + "\n" + "+" + "\n" + convert_int_to_ascii_char(record.letter_annotations["phred_quality"]) + "\n")
				else:
					options.inputfiletype=="fasta"
					out.write(">"+record.description + "\n" + str(record.seq) + "\n")

				options.idlist.remove(tmpseqid)
				break
	fh.close()
	out.close()
	options.input = output

	return

def get_subseq():
	#get sub sequence from defined position to another defined position
	fh=open(options.input)
	output="subseq." + options.inputfiletype
	out =open(output, "w")

	for record in SeqIO.parse(fh, options.inputfiletype):
		out.write(options.symbol + record.description + "\n")
		out.write( str(record.seq)[options.min - 1:options.max] + "\n")
		if options.inputfiletype == "fastq":
			out.write("+" + "\n" + convert_int_to_ascii_char(record.letter_annotations["phred_quality"][options.min -1:options.max]) + "\n")

	fh.close()
	out.close()
	options.input = output

	return

def filterbylength():

	if not options.min:
			print "Please provide the minimum readlength you want using the option --min"
			print "You can also provide the maximum readlength using option --max. By default there is no max readlength"
			exit(1)

	print "filtering reads by min and max"
	fh=open(options.input)
	output = "filterbylength_output." + options.inputfiletype
	out = open(output, "w")

	for record in SeqIO.parse(fh, options.inputfiletype):
		seqid = record.description.replace("|", "_").split()[0]
		ntseq = str(record.seq)
		if len(ntseq) >= options.min:
			if options.max:
				if len(ntseq) <= options.max:
					out.write(options.symbol + seqid + "\n" + ntseq + "\n")
					if options.inputfiletype == "fastq":
						out.write("+" + convert_int_to_ascii_char(record.letter_annotations["phred_quality"]) + "\n")
				else:
					continue
			else:
				out.write(options.symbol + seqid + "\n" + ntseq + "\n")
				if options.inputfiletype == "fastq":
					out.write("+" + convert_int_to_ascii_char(record.letter_annotations["phred_quality"]) + "\n")
		else:
			continue

	fh.close()
	options.input=output
	out.close()


def seq_by_interval():

	start=int(options.interval.split(",")[0])
	interval=int(options.interval.split(",")[1])
	fh=open(options.input)
	output="seq_by_interval." + options.inputfiletype
	out=open(output, "w")
	seq_counter=0; counter=0
	get_seq = False
	for record in SeqIO.parse(fh, options.inputfiletype):
		seq_counter+=1
		if seq_counter == start:
			get_seq==True; counter=seq_counter
		if seq_counter==counter:
			out.write(options.symbol + record.description + "\n" + str(record.seq) + "\n")
			if options.inputfiletype == "fastq":
				out.write("+" + "\n" + convert_int_to_ascii_char(record.letter_annotations["phred_quality"])+ "\n")

			counter+=interval   # substracting 1 because i added 1 above

		else:
			continue



	fh.close()
	out.close()
	options.input=output

	return


def getlength():
	fh=open(options.input)
	counter=0
	for record in SeqIO.parse(fh, options.inputfiletype):
		seqid = record.description.replace("|", "_").split()[0]
		ntseq = str(record.seq)
		counter+=1
		if options.min and options.max:
			if counter >= options.min and counter <= options.max:
				print seqid + "\t" + str(len(ntseq))
		#elif options.number:
		#	if counter <= options.number:
		#		print seqid + "\t" + str(len(ntseq))

		else:
			print seqid + "\t" + str(len(ntseq))

	return

def generate_random_number_with_seed(seedvalue=3, minvalue=0, maxvalue=10, total_to_generate=10):
	import random

	random.seed(seedvalue)
	random_numbers=[]
	while len(random_numbers) < total_to_generate:
		generated_value=random.randint(minvalue, maxvalue-1)
		if generated_value in random_numbers:
			continue

		random_numbers.append(generated_value)

	return random_numbers

def get_random_sequences_using_random_seed():

	output="generated_random." + options.inputfiletype
	out=open(output, "w")
	with open (options.input, 'rb') as fh:
	    records=list(SeqIO.parse(fh, options.inputfiletype))

	total_sequences=len(records)
	how_many_to_generate=total_sequences/2			#gets 50% of sequences randomly from fasta/fastq file, change this value to get desired percentage.

	random_numbers=generate_random_number_with_seed(maxvalue=total_sequences, total_to_generate=how_many_to_generate)
	
	for value in random_numbers:
		if options.inputfiletype=="fasta":
			out.write(">" + str(records[value].id) + "\n" + str(records[value].seq) + "\n")
		elif options.inputfiletype=="fastq":
			out.write("@" + records[value].id + "\n" + str(records[value].seq) + "\n" + "+" + convert_int_to_ascii_char(records[value].letter_annotations["phred_quality"]) + "\n")

	out.close()
	options.input=output
	return

def Mainprogram():

	if options.input:
		if options.interval and not options.random:
			seq_by_interval()
		elif options.random and not options.interval:
			get_random_sequences_using_random_seed()
		else:
			print "Options interval and random should not go togther. The result will be confusing."
			exit(1)

		if options.seqid:
			get_seqs_by_id()
		if options.getlength == True:
			getlength()
			exit(0)
		if filterbylength == True: #user will need to provide the option --min and --max
			filterbylength()
			exit(0)
		if options.subseq == True:
			get_subseq()

		if options.leftclip and options.rightclip:
			doleftrightclip()
		elif options.leftclip:
			doleftclip()
		elif options.rightclip:
			dorightclip()
		elif options.stats:
			get_stats()
			exit(0)
		elif options.reverse:
			reverse_sequence()
		elif options.reversecomp:
			reverse_complement()


		if options.output:
			#os.rename(options.input, options.output)
			pass
		else:
			read_fasta_fastq()
			#os.remove(options.input)
		options.input=options.temp_input

	else:
		print "No input file provided."
		exit(1)

	exit(0)


if __name__=="__main__":
	if options.input.endswith("fasta") or options.input.endswith("fa"):
		options.inputfiletype="fasta"
		options.symbol=">"
	elif options.input.endswith("fastq") or options.input.endswith("fq"):
		options.inputfiletype="fastq"
		options.symbol="@"
	else:
		print "Could not determine fasta or fastq file type. Fasta file should have .fasta or .fa file extension and fastq file should have .fastq or .fq file extension."
		exit(1)
	Mainprogram()
