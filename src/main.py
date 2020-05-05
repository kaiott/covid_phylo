import ncbi, config, iqtree, ete
from align_tools import SequenceAligner, Filter

def get_sequences():
    return ncbi.get_all_covid_nucleotide_seqs(cache_dir=config.CACHE_DIR)


def align_complete(data):
    my_aligner = SequenceAligner.from_tag(tag='complete', data=data)
    complete_filter = Filter(['complete genome']).all_filter
    my_aligner.add_filter(complete_filter)
    my_aligner.make_alignment()

    # for a copy without timestamp
    my_aligner.copy_aligned_file_unstamped()

def align_china(data):
    china_aligner = SequenceAligner.from_tag(tag='china', data=data)
    china_filter = Filter(['CHN', 'complete genome']).all_filter
    china_aligner.add_filter(china_filter)
    china_aligner.make_alignment()

    # for a copy without timestamp
    my_aligner.copy_aligned_file_unstamped()

def align_spain(data):
    aligner = SequenceAligner.from_tag(tag='spain', data=data)
    filt = Filter(['ESP', 'complete genome']).all_filter
    aligner.add_filter(filt)
    aligner.make_alignment()

    # for a copy without timestamp
    aligner.copy_aligned_file_unstamped()

def main():
    """
    DESCRIPTION:
    Main method of the program.
    :return: None.
    """
    print('retrieving records')
    result = get_sequences()
    
    align_complete(data=result)
    align_china(data=result)
    align_spain(data=result)


if __name__ == '__main__':
    main()
