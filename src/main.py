import re
import ncbi
import config
from Bio.Align.Applications import MafftCommandline


def mafft(origname=None, destname=None, route=None):
	"""
	DESCRIPTION:
	A function to interact with MAFFT through the system. To align the sequences in a given file.
	:param origname: [string] filename with all sequences in FASTA format.
	:param destname: [string] filename to store the aligned sequences.
	:param destname: [string] route to mafft binary.
	:return: Generates the file with the alignments of the sequences of the given file.
	"""
	# Obtain paths
	origname = config.FASTA_DIR / origname
	destname = config.FASTA_DIR / destname
	# Execute mafft
	mafft_cline = MafftCommandline(route, input=origname)
	stdout, stderr = mafft_cline()
	# Write result into file
	file = open(destname, 'w')
	file.write(stdout)
	file.close()


def filter_complete_genome(header_line):
	"""
	DESCRIPTION:
	Checks whether the header line of the fasta format contains the string "complete genome"
	:param header_line: [string] header line to be checked
	:return: [boolean] true iff the header line contains "complete genome"
	"""
	# no need for regex matching now, but maybe later
	m = re.search('complete genome', header_line)
	if m is not None:
		return True
	return False


def records_to_fasta(records, fasta_dir=None, filename='records.txt', header_filter=None):
	"""
	DESCRIPTION:
	A function to change between a set of records in a list and a normalized fasta file with all the sequences. Be
	notice that all data is downloaded. Not specific fragments.
	:param header_filter: function that that checks condition on the header_line
	:param records: [list] records read.
	:param fasta_dir: [pathlib] route to the folder in which the sequences are to be stored.
	:param filename: [string] name of the file in which the sequences are stores.
	:return: Generates the file with the records.
	"""
	# Check we have records or not
	if records is None:
		return None
	fasta_result = {}
	content = ''
	# Save all the records in fasta format in a string
	for record in records:
		fasta = record.format('fasta')
		[header_line, sequence] = fasta.split('\n', 1)
		# check if header_line passes the specified filter and only then add it to result and content of file
		if header_filter is None or header_filter(header_line):
			if fasta_dir is not None:
				content += fasta
			fasta_result[header_line] = sequence
	# Write the string in a file
	if fasta_dir is not None:
		fasta_dir.mkdir(exist_ok=True)
		f = open(fasta_dir / filename, 'w')
		f.write(content)
		f.close()

	return fasta_result


def main():
	"""
	DESCRIPTION:
	Main method of the program.
	:return: None.
	"""
	result = ncbi.get_all_covid_nucleotide_seqs(cache_dir=config.CACHE_DIR)
	records = result.get('seqrecords')
	records_to_fasta(records, config.FASTA_DIR, 'complete_genomeDNA.txt', filter_complete_genome)
	# I commented this, because I don't have mafft installed (and even if I did, I think I'd get an error, because on
	# Windows the paths are different, maybe use Pathlib to fix that
	# mafft_route = '/usr/bin/mafft'
	# mafft(origname='complete_genomeDNA.txt', destname='complete_gene_align.txt', route=mafft_route)


if __name__ == '__main__':
	main()
