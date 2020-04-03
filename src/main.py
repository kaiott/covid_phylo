import ncbi
import config

if __name__ == '__main__':
	result = ncbi.get_all_covid_nucleotide_seqs(cache_dir=config.CACHE_DIR)
	print(result.keys())
	record = result.get('seqrecords')[0]
	print(record.id)
	print(str(record.seq))