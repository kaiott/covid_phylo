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
	timestamp = datetime.fromtimestamp(int(result['request_timestamp'])).strftime("%Y%m%d%H%M%S")

	# all from China that are of the complete genome
	china_aligner = SequenceAligner.from_tag('china', timestamp)
	china_filter = Filter(['CHN', 'complete genome']).all_filter
	china_aligner.set_records(records)
	china_aligner.add_filter(china_filter)
	china_aligner.make_alignment()

	# complete genomes either from Spain or England that have Covid-19-2 in their description
	second_aligner = SequenceAligner.from_tag('test2', timestamp)
	type_filter = Filter(['complete genome', 'Covid-19-2']).all_filter
	country_filter = Filter(['ESP', 'USA']).any_filter
	second_aligner.set_records(records)
	second_aligner.add_filter(type_filter)
	second_aligner.add_filter(country_filter)
	second_aligner.make_alignment()

	# all genomes of protein 'N' neither from US nor China
	third_aligner = SequenceAligner.from_tag('test3', timestamp)
	type_filter = Filter(['N protein or whatever the name was']).all_filter
	country_filter = Filter(['CHN', 'USA']).none_filter
	third_aligner.set_records(records)
	third_aligner.add_filter(type_filter)
	third_aligner.add_filter(country_filter)
	third_aligner.make_alignment()

	# get path of their output files easily
	print(third_aligner.get_aligned_filename())
	# prints filename in fasta folder


if __name__ == '__main__':
	main()