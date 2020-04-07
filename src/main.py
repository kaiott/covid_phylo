import os
import re
from pathlib import Path
import ncbi
import config
from Bio.Align.Applications import MafftCommandline


def mafft(origname=None, destname=None, route=None):
    """
    DESCRIPTION:
    A function to interactuate with MAFFT through the system. To align the sequences in a given file.
    :param origname: [string] filename with all sequences in FASTA format.
    :param destname: [string] filename to store the aligned sequences.
    :param destname: [string] route to mafft binary.
    :return: Generates the file with the alignments of the sequences of the given file.
    """
    # Obtain paths
    fasta_path = Path(Path(__file__).parent.absolute()).parent.absolute() / Path('covid_phylo_data/fasta')
    origname = fasta_path / origname
    destname = fasta_path / destname
    # Execute mafft
    mafft_cline = MafftCommandline(route, input=origname)
    stdout, stderr = mafft_cline()
    # Write result into file
    file = open(destname, 'w')
    file.write(stdout)
    file.close()


def records_to_fasta(records, fasta_dir=None, filename='records.txt'):
    """
    DESCRIPTION:
    A function to change between a set of records in a list and a normalized fasta file with all the sequences. Be
    notice that all data is downloaded. Not specific fragments.
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
        if fasta_dir is not None:
            content += fasta

        [header_line, sequence] = fasta.split('\n', 1)
        fasta_result[header_line] = sequence
    # Write the string in a file
    if fasta_dir is not None:
        fasta_dir.mkdir(exist_ok=True)
        f = open(fasta_dir / filename, 'w')
        f.write(content)
        f.close()

    return fasta_result


def select_fasta(orf=None, origname=None, destname=None):
    """
    DESCRIPTION:
    Given a FASTA file with all the information, generates another FASTA with only the sequences of the selected part of
    the genome.
    :param orf: [string] sequence of letters indicating the ORF.
    :param origname: [string] filename with all sequences in FASTA format.
    :param destname: [string] filename to stores the selected sequences in fasta format. In the same folder.
    :return: Generates the file with the records fo the select region.
    """
    # Take for granted the location of all fasta files
    fasta_path = Path('covid_phylo_data/fasta')
    or_path = Path(Path(__file__).parent.absolute()).parent.absolute()

    # Read from a file, select and write into the other file
    file_r = open(or_path / fasta_path / origname)
    file_w = open(or_path / fasta_path / destname, 'w')
    records = file_r.read()
    records = records.split('>')
    for record in records:
        record = record.split('\n')
        header = re.split(' |, ', record[0])
        orf_a = '(' + orf + ')'
        if orf not in header and orf_a not in header or ';' in record[0]:
            continue
        else:
            header = '\n' + '>' + record[0]
            seq = ''.join(record[1::])
            n = 60
            seq = '\n' + '\n'.join([seq[i:i+n] if n+i < len(seq) else seq[i::] for i in range(0, len(seq), n)])
            file_w.write(header)
            file_w.write(seq)
    file_r.close()
    file_w.close()


def main():
    """
    DESCRIPTION:
    Main method of the program.
    :return: None.
    """
    result = ncbi.get_all_covid_nucleotide_seqs(cache_dir=config.CACHE_DIR)
    records = result.get('seqrecords')
    records_to_fasta(records, config.FASTA_DIR, 'IUPACAmbiguousDNA.txt')
    select_fasta(orf='S', origname='IUPACAmbiguousDNA.txt', destname='S_gene_seq.fasta')
    mafft_route = '/usr/bin/mafft'
    mafft(origname='S_gene_seq.fasta',  destname='S_gene_align.fasta', route=mafft_route)


if __name__ == '__main__':
    main()
