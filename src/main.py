
import ncbi
import config


def records_to_fasta(records, fasta_dir=None, filename='records.txt'):
	if records is None:
		return None
	
	fasta_result = {}
	content = ''
	for record in records:
		fasta = record.format('fasta')
		if fasta_dir is not None:
			content += fasta
		
		[header_line, sequence] = fasta.split('\n', 1)
		print(header_line)
		fasta_result[header_line] = sequence

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
    :return: None. Execution of the instructions.
    """
    result = ncbi.get_all_covid_nucleotide_seqs(cache_dir=config.CACHE_DIR)
    records = result.get('seqrecords')
    records_to_fasta(records, config.FASTA_DIR, 'IUPACAmbiguousDNA.txt')


if __name__ == '__main__':
    main()
