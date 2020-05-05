import config
import json
import numpy as np
import os
from datetime import datetime
from Bio.Align.Applications import MafftCommandline


class SequenceAligner:
    # configuration of the SequenceAligner class, can
    # be changed if needed

    unaligned_pattern = '{tag}_{file_id}_unaligned'
    aligned_pattern = '{tag}_{file_id}_aligned'
    information_pattern = '{tag}_information.txt'

    def __init__(self, tag, file_id, records=None, already_aligned_file_id=None,
                 already_aligned_sequence_ids=None):
        """
        DESCRIPTION:
        Constructor of the SequenceAligner class
        :param tag: [string] desired name for the selection of sequences
        :param file_id: [string] id of this alignment
        :param records: [list] list of records to align
        :param already_aligned_file_id: [string] id of previous alignment we want to work with
        :param already_aligned_sequence_ids: [list] list of ids of sequences already aligned in mentioned previous alignment
        :return: [SequenceAligner] the created object
        """
        if already_aligned_sequence_ids is None:
            already_aligned_sequence_ids = []
        self.filters = []
        self.tag = tag
        self.file_id = file_id
        self.unfiltered_records = records
        self.already_aligned_file_id = already_aligned_file_id
        self.already_aligned_sequence_ids = already_aligned_sequence_ids

    @staticmethod
    def from_tag(tag, data):
        """
        DESCRIPTION:
        Creates a SequenceAligner object based on an already existing alignment or
        in case that doesn't exist creates a new one
        :param tag: [string] name of the alignment selection to build on
        :param data: [dictionary] the data including a timestamp and the records
        :return: [SequenceAligner] the created object
        """
        already_aligned_file_id, already_aligned_sequence_ids = SequenceAligner.get_actual(tag)
        if already_aligned_file_id is None or already_aligned_sequence_ids == []:
            print(f'No previous alignment of {tag} found. Alignment will be done from scratch')
        else:
            print(
                f'Found previous alignment of {tag} with id {already_aligned_file_id} and {len(already_aligned_sequence_ids)} aligned sequences')

        timestamp = datetime.fromtimestamp(int(data.get('request_timestamp'))).strftime("%Y%m%d%H%M%S")
        records = data.get('seqrecords')
        return SequenceAligner(tag, timestamp, records=records, already_aligned_file_id=already_aligned_file_id,
                               already_aligned_sequence_ids=already_aligned_sequence_ids)

    def add_filter(self, description_filter):
        """
        DESCRIPTION:
        Adds a filter to filters
        :param description_filter: [function] filter that will be added
        """
        self.filters.append(description_filter)

    def make_alignment(self, make_copy=False):
        """
        DESCRIPTION:
        Aligns unaligned sequences and writes them to a file. Updates the information file
        that says which sequences have already been aligned
        """
        sequence_ids_written = self._write_filtered_records_to_file()
        if len(sequence_ids_written) == 0:
            print('no unaligned sequences')
            return

        print(f'will align {len(sequence_ids_written)} new sequences')

        if self.already_aligned_file_id is None:
            print('aligning from scratch')
            self._align_from_scratch()
        else:
            print('adding to previous alignment')
            self._align_from_existing()

        if make_copy:
            self._copy_aligned_file_unstamped()

        # at this point assume alignment done successfully
        # update meta information
        sequence_ids_written += self.already_aligned_sequence_ids
        info_dict = {'last_id': self.file_id,
                     'aligned_ids': sequence_ids_written
                     }
        info_filename = SequenceAligner.information_pattern.format(tag=self.tag)
        with open(config.FASTA_DIR / info_filename, 'w') as file:
            file.write(json.dumps(info_dict))

    def _copy_aligned_file_unstamped(self):
        in_filename = config.FASTA_DIR / self.get_aligned_filename()
        out_filename = config.FASTA_DIR / f'{self.tag}_aligned'
        with open(in_filename, 'r') as in_file, open(out_filename, 'w') as out_file:
            content = in_file.read()
            out_file.write(content)

    def get_aligned_filename(self):
        """
        DESCRIPTION:
        Returns the name of the file with the aligned sequences
        :return: [string] the name of the output file
        """
        return SequenceAligner.aligned_pattern.format(tag=self.tag,
                                                      file_id=self.already_aligned_file_id)

    @staticmethod
    def get_actual(tag):
        """
        DESCRIPTION:
        Loads information of a previous alignment of this selection
        :param tag: [string] name of the alignment to check
        :return: [string, list] id of the last alignment done with this tag
        and list of id's of sequences already aligned by that
        """
        try:
            filename = SequenceAligner.information_pattern.format(tag=tag)
            with open(config.FASTA_DIR / filename, 'r') as file:
                json_info = json.load(file)
                return json_info['last_id'], json_info['aligned_ids']
        except FileNotFoundError:
            return None, []

    def get_filtered_records(self):
        """
        DESCRIPTION:
        Filters this object's records according to whether these records fulfill the
        self.filters criteria and whether they have already been aligned
        :return: [list] list of records that pass the filters
        """
        if self.unfiltered_records is None:
            return None

        result_records = []
        for record in self.unfiltered_records:
            if record.id in self.already_aligned_sequence_ids:
                continue

            if all([description_filter(record.description) for description_filter in self.filters]):
                result_records.append(record)

        return result_records

    def set_records(self, records):
        """
        DESCRIPTION:
        Sets this object's filters
        """
        self.unfiltered_records = records

    def _write_filtered_records_to_file(self):
        """
        DESCRIPTION:
        Filters self.unfiltered_records and writes the records that pass into a
        file in the fasta format
        :return: [list] list of the record's ids that were written to the file
        """
        records = self.get_filtered_records()
        if len(records) == 0:
            print('no records to write to file, done nothing')
            return []

        sequence_ids_written = []
        output_file = config.FASTA_DIR / SequenceAligner.unaligned_pattern.format(
            tag=self.tag, file_id=self.file_id)
        with open(output_file, 'w') as file:
            for record in records:
                file.write(record.format('fasta'))
                sequence_ids_written.append(record.id)

        return sequence_ids_written

    def _align_from_scratch(self):
        """
        DESCRIPTION:
        Aligns unaligned sequences in case there was no previous alignment and
        writes the alignment to a file
        """
        origname = config.FASTA_DIR / SequenceAligner.unaligned_pattern.format(
            tag=self.tag, file_id=self.file_id)
        destname = config.FASTA_DIR / SequenceAligner.aligned_pattern.format(
            tag=self.tag, file_id=self.file_id)
        print('Executing sequences alignment...')
        mafft_cline = MafftCommandline(config.MAFFT_DIR, input=origname)
        stdout, stderr = mafft_cline()
        print('Alignment completed')
        # Write result into file
        file = open(destname, 'w')
        file.write(stdout)
        file.close()
        print('Alignment saved')

    def _align_from_existing(self):
        """
        DESCRIPTION:
        Add unaligned sequences to an existing alignment and writes the alignment
        to a file
        """
        unaligned_file = config.FASTA_DIR / SequenceAligner.unaligned_pattern.format(
            tag=self.tag, file_id=self.file_id)
        aligned_file = config.FASTA_DIR / SequenceAligner.aligned_pattern.format(
            tag=self.tag, file_id=self.already_aligned_file_id)
        output_file = config.FASTA_DIR / SequenceAligner.aligned_pattern.format(
            tag=self.tag, file_id=self.file_id)
        print('Executing sequences alignment...')
        command = f'mafft --add {unaligned_file} --reorder {aligned_file} > {output_file}'
        os.system(command)
        print('Alignment completed')


class Filter:
    def __init__(self, key_words):
        """
        DESCRIPTION:
        Constructor of the Filter class
        :param key_words: [list] list of strings the filters use
        :return: [Filter] the created object
        """
        self.key_words = key_words

    def all_filter(self, string_to_check):
        """
        DESCRIPTION:
        Checks whether the a string contains all of the keywords
        :param string_to_check: [string] the string to check
        :return: [boolean] true iff the string to check contains all the keywords
        """
        return all([key_word in string_to_check for key_word in self.key_words])

    def any_filter(self, string_to_check):
        """
        DESCRIPTION:
        Checks whether the a string contains any of the keywords
        :param string_to_check: [string] the string to check
        :return: [boolean] true iff the string to check contains at least one of the keywords
        """
        return any([key_word in string_to_check for key_word in self.key_words])

    def none_filter(self, string_to_check):
        """
        DESCRIPTION:
        Checks whether the a string contains none of the keywords
        :param string_to_check: [string] the string to check
        :return: [boolean] true iff the string to check contains none of the keywords
        """
        return not any([key_word in string_to_check for key_word in self.key_words])


def _get_aligned_content_by_tag(tag):
    file = config.FASTA_DIR / f'{tag}_aligned'
    try:
        with open(file, 'r') as f:
            return f.read()

    except FileNotFoundError:
        print(f'no alignment with tag {tag}')
        return None


def aligned_records_by_tag(tag):
    content = _get_aligned_content_by_tag(tag)
    if content is None:
        return None

    raw_records = content.split('>')[1:]
    records = []
    for raw_record in raw_records:
        header, sequence = raw_record.split('\n', 1)
        records.append({'header': header, 'sequence': sequence.replace('\n', '')})

    return records


def analyse_alignment(aligned_records):
    sequences = [record['sequence'] for record in aligned_records]
    lengths = [len(seq) for seq in sequences]
    if max(lengths) != min(lengths):
        print('sequences don\'t have same length')
        return

    num_gaps = np.zeros(lengths[0], dtype=int)
    num_variation_det = np.zeros(lengths[0], dtype=int)
    num_variation_all = np.zeros(lengths[0], dtype=int)
    for site in range(lengths[0]):
        num_nucleotides_det = {}
        num_nucleotides_undet = {}
        for seq in sequences:
            c = seq[site]
            if c == '-':
                num_gaps[site] += 1
            elif c == 'a' or c == 't' or c == 'g' or c == 'c':
                num_nucleotides_det[c] = True
            else:
                num_nucleotides_undet[c] = True

        num_variation_det[site] = len(num_nucleotides_det)
        num_variation_all[site] = num_variation_det[site] + len(num_nucleotides_undet)

    return num_gaps, num_variation_det, num_variation_all

