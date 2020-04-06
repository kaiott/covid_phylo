
import config
from pprint import pprint
import time
import shelve
import io
import requests
from Bio import SeqIO

ENTREZ_COVID_SEARCH_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term=txid2697049[Organism:noexp]&retmode=json&retmax=10000'

ENTREZ_NUCL_DOWNLOAD_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={uids}&retmode=text&rettype={format}'


def _get_raw_sequence(uid, cache_dir=None, format='gb'):

	if cache_dir is not None:
		shelve_path = cache_dir / config.RAW_SEQUENCE_SHELVE_FNAME
		shelved_raw_seqs = shelve.open(str(shelve_path))
		if uid in shelved_raw_seqs:
			print('uid: ' + str(uid) + ' already in cache')
			return shelved_raw_seqs[uid]

	print('uid: ' + str(uid) + ' not in cache, proceding to download')
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
	uids = ncbi_search_result['idlist']

	if n_seqs_found > len(uids):
		msg = 'Some sequences were not retrieved, you should implement the search with usehistory'
		raise NotImplementedError(msg)

	seq_records = []
	for uid in uids:
		raw_seq = _get_raw_sequence(uid, cache_dir=cache_dir, format='gb')
		fhand = io.StringIO(raw_seq)
		record = list(SeqIO.parse(fhand, 'gb'))[0]
		seq_records.append(record) # used to be .append(seq_records) resulting in the list containing itself an exponetial number of times...

	search_result = {'request_timestamp': time.time(),
					 'seqrecords': seq_records
					}
	return search_result

