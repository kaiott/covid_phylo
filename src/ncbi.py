import config
from pprint import pprint
import time
import shelve
import io, sys
import requests
from Bio import SeqIO

ENTREZ_COVID_SEARCH_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term=txid2697049[Organism:noexp]&retmode=json&retmax=10000'

ENTREZ_NUCL_DOWNLOAD_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={uids}&retmode=text&rettype={format}'



# update_progress() : Displays or updates a console progress bar
## Accepts a float between 0 and 1. Any int will be converted to a float.
## A value under 0 represents a 'halt'.
## A value at 1 or bigger represents 100%
def update_progress(progress):
    barLength = 30 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\n"
    if progress < 0:
        progress = 0
        status = "Halt...\n"
    if progress >= 1:
        progress = 1
        status = "Done\n"
    block = int(round(barLength*progress))
    text = f'\r\rRetrieving records: [{ "#"*block }{ "-"*(barLength-block) }] {(progress*100):.1f}% {status}'
    sys.stdout.write(text)
    sys.stdout.flush()


def _get_raw_sequence(uid, cache_dir=None, format='gb'):
	if cache_dir is not None:
		shelve_path = cache_dir / config.RAW_SEQUENCE_SHELVE_FNAME
		shelved_raw_seqs = shelve.open(str(shelve_path))
		if uid in shelved_raw_seqs:
			# print('uid: ' + str(uid) + ' already in cache')
			return shelved_raw_seqs[uid]

	#print('uid: ' + str(uid) + ' not in cache, proceding to download')
	nucl_download_url = ENTREZ_NUCL_DOWNLOAD_URL.format(uids=','.join([uid]), format=format)

	response = requests.get(nucl_download_url)
	if response.status_code != 200:
		msg = 'Something went wrong downloading the nucleotide sequences. '
		msg += 'response status: {response.status_code}'
		raise RuntimeError(msg)

	raw_seq = response.text

	if cache_dir is not None:
		shelved_raw_seqs[uid] = raw_seq

	return raw_seq


def get_all_covid_nucleotide_seqs(cache_dir=None):
	if cache_dir is not None:
		cache_dir.mkdir(exist_ok=True)

	response = requests.get(ENTREZ_COVID_SEARCH_URL)
	if response.status_code != 200:
		msg = 'Something went wrong searching for the SARS-CoV-2 nucleotide sequences. '
		msg += 'response status: {response.status_code}'
		raise RuntimeError(msg)

	ncbi_search_result = response.json()['esearchresult']

	n_seqs_found = int(ncbi_search_result['count'])

	print('found ' + str(n_seqs_found) + ' sequences')

	uids = ncbi_search_result['idlist']

	if n_seqs_found > len(uids):
		msg = 'Some sequences were not retrieved, you should implement the search with usehistory'
		raise NotImplementedError(msg)

	seq_records = []
	for uid in uids:
		raw_seq = _get_raw_sequence(uid, cache_dir=cache_dir, format='gb')
		#print('retrieving record ' + str(len(seq_records)+1) + '/' + str(n_seqs_found))
		update_progress((len(seq_records)+1)/n_seqs_found)
		fhand = io.StringIO(raw_seq)
		record = list(SeqIO.parse(fhand, 'gb'))[0]
		seq_records.append(record)  # used to be .append(seq_records) resulting in the list containing itself an exponetial number of times...

	search_result = {'request_timestamp': time.time(),
					 'seqrecords': seq_records
					 }
	return search_result
