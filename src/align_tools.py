import config
import json
import re
import os
from pathlib import Path
from Bio.Align.Applications import MafftCommandline
	
class SequenceAligner:
	# configuration of the SequenceAlinger class, can 
	# be changed if needed
	used_dir = config.FASTA_DIR
	mafft_dir = '/usr/bin/mafft'
	unaligned_pattern = '{selection_specifier}_{file_id}_unaligned'
	aligned_pattern = '{selection_specifier}_{file_id}_aligned'
	information_pattern = '{selection_specifier}_information.txt'

	def __init__(self, selection_specifier, file_id, records = None, already_aligned_file_id = None, already_aligned_sequence_ids = []):
		"""
		DESCRIPTION:
		Constructor of the SequenceAligner class
		:param selection_specifier: [string] desired name for the selection of sequences
		:param file_id: [string] id of this alignment
		:param records: [list] list of records to align
		:param already_aligned_file_id: [string] id of previous alignment we want to work with
		:param already_aligned_sequence_ids: [list] list of ids of sequences already aligned in mentioned previous alignment
		:return: [SequenceAligner] the created object
		"""
		self.filters = []
		self.selection_specifier = selection_specifier
		self.file_id = file_id
		self.unfiltered_records = records
		self.already_aligned_file_id = already_aligned_file_id
		self.already_aligned_sequence_ids = already_aligned_sequence_ids



	def add_filter(self, description_filter):
		"""
		DESCRIPTION:
		Adds a filter to filters
		:param description_filter: [function] filter that will be added
		"""
		self.filters.append(description_filter)


	@staticmethod
	def get_actual(selection_specifier):
		"""
		DESCRIPTION:
		Loads information of a previous alignment of this selection
		:param selection_specifier: [string] name of the alignment to check
		:return: [string, list] id of the last alignment done with this selection_specifier
		and list of id's of sequences already aligned by that
		"""
		try:
			filename = SequenceAligner.information_pattern.format(selection_specifier=selection_specifier)
			with open(SequenceAligner.used_dir / filename, 'r') as file:
				json_info = json.load(file)
				return json_info['last_id'], json_info['aligned_ids']
		except FileNotFoundError:
			return None, []	


	@staticmethod
	def from_most_recent_existing_alignment(selection_specifier, file_id):
		"""
		DESCRIPTION:
		Creates a SequenceAligner object based on an already existing alignment or 
		in case that doesn't exist creates a new one
		:param selection_specifier: [string] name of the alignment selection to build on
		:param file_id: [string] id of this SequenceAligner object
		:return: [SequenceAligner] the created object
		"""
		already_aligned_file_id, already_aligned_sequence_ids = SequenceAligner.get_actual(selection_specifier)
		if already_aligned_file_id is None or already_aligned_sequence_ids == []:
			print(f'No previous alignment of {selection_specifier} found. Alignment will be done from scratch')
		else:
			print(f'Found previous alignment of {selection_specifier} with id {already_aligned_file_id} and {len(already_aligned_sequence_ids)} aligned sequences')
		
		return SequenceAligner(selection_specifier, file_id, already_aligned_file_id=already_aligned_file_id, already_aligned_sequence_ids=already_aligned_sequence_ids)


	def make_alignment(self):
		"""
		DESCRIPTION:
		Aligns unaligned sequences and writes them to a file. Updates the information file
		that says which sequences have already been aligned
		"""
		sequence_ids_written = self.write_filtered_records_to_file()
		if len(sequence_ids_written) == 0:
			print('no unaligned sequences')
			return

		print(f'will align {len(sequence_ids_written)} new sequences')

		if self.already_aligned_file_id is None:
			print('aligning from scratch')
			self.align_from_scratch()
		else:
			print('adding to previous alignment')
			self.align_from_existing()

		# at this point assume alignment done successfully
		# update meta information
		sequence_ids_written += self.already_aligned_sequence_ids
		info_dict = {'last_id': self.file_id,
					 'aligned_ids': sequence_ids_written
					}
		info_filename = SequenceAligner.information_pattern.format(selection_specifier=self.selection_specifier)
		with open(SequenceAligner.used_dir / info_filename, 'w') as file:
			file.write(json.dumps(info_dict))



	def align_from_scratch(self):
		"""
		DESCRIPTION:
		Aligns unaligned sequences in case there was no previous alignment and 
		writes the alignment to a file
		"""
		origname = SequenceAligner.used_dir / SequenceAligner.unaligned_pattern.format(selection_specifier=self.selection_specifier, file_id=self.file_id)
		destname = SequenceAligner.used_dir / SequenceAligner.aligned_pattern.format(selection_specifier=self.selection_specifier, file_id=self.file_id)
		# Execute mafft
		print('Executing sequences alignment...')
		mafft_cline = MafftCommandline(SequenceAligner.mafft_dir, input=origname)
		stdout, stderr = mafft_cline()
		print('Alignment completed')
		# Write result into file
		file = open(destname, 'w')
		file.write(stdout)
		file.close()
		print('Alignment saved')


	def align_from_existing(self):
		"""
		DESCRIPTION:
		Add unaligned sequences to an existing alignment and writes the alignment
		to a file
		"""
		unaligned_file = SequenceAligner.used_dir / SequenceAligner.unaligned_pattern.format(selection_specifier=self.selection_specifier, file_id=self.file_id)
		aligned_file = SequenceAligner.used_dir / SequenceAligner.aligned_pattern.format(selection_specifier=self.selection_specifier, file_id=self.already_aligned_file_id)
		output_file = SequenceAligner.used_dir / SequenceAligner.aligned_pattern.format(selection_specifier=self.selection_specifier, file_id=self.file_id)
		print('Executing sequences alignment...')
		command = f'mafft --add {unaligned_file} --reorder {aligned_file} > {output_file}'
		os.system(command)
		print('Alignment completed')

	def set_records(self, records):
		"""
		DESCRIPTION:
		Sets this object's filters
		"""
		self.unfiltered_records = records


	def write_filtered_records_to_file(self):
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
		output_file = SequenceAligner.used_dir / SequenceAligner.unaligned_pattern.format(selection_specifier=self.selection_specifier, file_id=self.file_id)
		with open(output_file, 'w') as file:
			for record in records:
				file.write(record.format('fasta'))
				sequence_ids_written.append(record.id)

		return sequence_ids_written


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
				# omitting this record, cause already accounted for
				continue

			if all([filt(record.description) for filt in self.filters]):
				result_records.append(record)

		return result_records


	def get_aligned_filename(self):
		"""
		DESCRIPTION:
		Returns the name of the file with the aligned sequences
		:return: [string] the name of the output file
		"""
		return SequenceAligner.aligned_pattern.format(selection_specifier=self.selection_specifier, file_id=self.already_aligned_file_id)


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
