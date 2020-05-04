import re
import os
import json
import ncbi, config, iqtree, ete
from datetime import datetime
from pathlib import Path
from Bio.Align.Applications import MafftCommandline
from align_tools import SequenceAligner, Filter


def main():
	"""
	DESCRIPTION:
	Main method of the program.
	:return: None.
	"""
	print('retrieving records')
	result = ncbi.get_all_covid_nucleotide_seqs(cache_dir=config.CACHE_DIR)
	records = result.get('seqrecords')
	print('number of records retrieved: ' + str(len(records)))
	timestamp = datetime.fromtimestamp(int(result['request_timestamp'])).strftime("%Y%m%d%H%M%S")

	# update the complete alignment
	my_aligner = SequenceAligner.from_most_recent_existing_alignment('complete', timestamp)
	complete_filter = Filter(['complete genome']).all_filter
	my_aligner.set_records(records)
	my_aligner.add_filter(complete_filter)
	my_aligner.make_alignment()

	# iqtree visualizing tree	
	print('Select the best alignments')
	origname = my_aligner.get_aligned_filename()
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
