
import ncbi
import config


def main():
    """
    DESCRIPTION:
    Main method of the program.
    :return: None. Execution of the instructions.
    """
    result = ncbi.get_all_covid_nucleotide_seqs(cache_dir=config.CACHE_DIR)
    record = result.get('seqrecords')[0]
    print(record)


if __name__ == '__main__':
    main()
