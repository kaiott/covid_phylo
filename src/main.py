import re
import os
import json
import ncbi, config, iqtree, ete
from datetime import datetime
from pathlib import Path
from Bio.Align.Applications import MafftCommandline



def mafft_add(aligned_file, unaligned_file, output_file):
	aligned_file = config.FASTA_DIR / aligned_file
	unaligned_file = config.FASTA_DIR / unaligned_file
	output_file = config.FASTA_DIR / output_file
	command = f'mafft --add {unaligned_file} --reorder {aligned_file} > {output_file}'
	os.system(command)

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
	print('Executing sequences alignment...')
	mafft_cline = MafftCommandline(route, input=origname)
	print('before cline')
	stdout, stderr = mafft_cline()
	print('Alignment completed')
	# Write result into file
	file = open(destname, 'w')
	file.write(stdout)
	file.close()
	file = open(config.FASTA_DIR / 'derroutput', 'w')
	file.write(stderr)
	file.close()
	print('Alignment saved')


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


def get_fasta_info(fasta_dir):
	try:
		with open(fasta_dir / 'alignment_information.txt', 'r') as file:
			json_info = json.load(file)
			return json_info['last_id'], json_info['aligned_ids']
	except FileNotFoundError:
		return None, None


def write_new_records_to_fasta(records, file_id, fasta_dir=None, accounted_for=None, header_filter=None, max_used=None):
	if records is None or file_id is None:
		return None # or better raise some exception

	fasta_result = {}
	content = ''
	if accounted_for is None:
		new_accounted_for = []
	else:
		new_accounted_for = accounted_for[:]

	for record in records:
		if max_used is not None and len(fasta_result) >= max_used:
			break
		if accounted_for is not None and record.id in accounted_for:
			# omitting this record, cause already accounted for
			continue
		fasta = record.format('fasta')
		[header_line, sequence] = fasta.split('\n', 1)
		if header_filter is None or header_filter(header_line):
			if fasta_dir is not None:
				content += fasta
				new_accounted_for.append(record.id)
			fasta_result[header_line] = sequence

	if len(fasta_result)==0:
		return fasta_result

	# Write the string in a file
	if fasta_dir is not None:
		# write fasta
		filename = 'unaligned_complete_'
		if accounted_for is None:
			filename += 'all_'
		else:
			filename += 'new_'
		filename += file_id
		fasta_dir.mkdir(exist_ok=True)
		f = open(fasta_dir / filename, 'w')
		f.write(content)
		f.close()

		# update meta information
		info_dict = {'last_id': file_id,
					 'aligned_ids': new_accounted_for
					}
		with open(fasta_dir / 'alignment_information.txt', 'w') as file:
			file.write(json.dumps(info_dict))

	return fasta_result





def records_to_fasta(records, fasta_dir=None, filename='records.txt', header_filter=None, max_used=None):
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

	return write_new_records_to_fasta(records=records, file_id=filename, fasta_dir=fasta_dir, accounted_for=None, header_filter=header_filter, max_used=max_used)
	# Check we have records or not
	if records is None:
		return None
	fasta_result = {}
	content = ''
	# Save all or maximal max_used records in fasta format in a string
	for record in records:
		if max_used is not None and len(fasta_result) >= max_used:
			break
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





def update_alignment():
	print('retrieving records')
	result = ncbi.get_all_covid_nucleotide_seqs(cache_dir=config.CACHE_DIR)
	records = result.get('seqrecords')
	print('number of records retrieved: ' + str(len(records)))
	timestamp = datetime.fromtimestamp(int(result['request_timestamp'])).strftime("%Y%m%d%H%M%S")
	old_timestamp, accounted_for_genomes = get_fasta_info(config.FASTA_DIR)
	if old_timestamp is None or accounted_for_genomes is None:
		print('we haven\'t made any previous alignments')
		# we haven't made any previous alignments
		result = write_new_records_to_fasta(records, timestamp, fasta_dir=config.FASTA_DIR, accounted_for=accounted_for_genomes, 
											header_filter=filter_complete_genome, max_used=None)
		print('numbers of records used: ' + str(len(result)))
		# align the genomes using mafft
		print('aligning with mafft')
		if os.name == 'posix':
			mafft_route = '/usr/bin/mafft'
			if True or not os.path.exists(config.FASTA_DIR / destname):
				print('Alignment not found in fasta folder')
				print('Starting alignment')
				mafft(origname='unaligned_complete_all_' + timestamp, destname='aligned_complete_' + timestamp, route=mafft_route)
				print('Alignment finished')
			else:
				print('Alignment found in fasta folder')
		else:
			pass


	else:
		print('we have made previous alignments and can thus add our new sequences to that alignment')
		# we have made previous alignments and can thus add our new sequences to that alignment
		result = write_new_records_to_fasta(records, timestamp, fasta_dir=config.FASTA_DIR, accounted_for=accounted_for_genomes, 
											header_filter=filter_complete_genome, max_used=None)
		print('numbers of records used: ' + str(len(result)))
		if (len(result) != 0):
			# new sequences have been retrieved that need to be aligned
			print('Adding new sequences to alignment ' + old_timestamp)
			mafft_add(aligned_file='aligned_complete_' + old_timestamp, unaligned_file='unaligned_complete_new_' + timestamp, output_file='aligned_complete_' + timestamp)
		else:
			print('No new sequences to add: No changes made')


def main():
	"""
	DESCRIPTION:
	Main method of the program.
	:return: None.
	"""
	# retrieve records from database or cache

	update_alignment()
	recent_timestamp, accounted_for = get_fasta_info(config.FASTA_DIR)

	# iqtree visualizing tree	
	print('Select the best alignments')
	origname = 'aligned_complete_' + recent_timestamp
	destname = 'selection.txt'
	n_genomes = 100

	if not Path.exists(config.TREE_DIR / destname.split('.')[0] / destname):
		print('Creating selection file')
	iqtree.align_selector(origname, destname, n_genomes)
	print('Looking for tree inference')
	folder = Path(config.TREE_DIR / destname.split('.')[0])
	if len([x for x in folder.iterdir()]) < 2:
		print('Tree not found in tree folder')
		print('Starting tree inference')
		iqtree.tree_creator(destname)
		print('Inference finished')
	else:
		print('Tree found in tree folder')
	print('Generating tree visualization')
	tree_route = config.TREE_DIR / destname.split('.')[0]
	ete.tree_viewer(tree_route)


if __name__ == '__main__':
	main()
