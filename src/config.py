
from pathlib import Path
# hola
PROJECT_DIR = Path('.').resolve().parent

BASE_DIR = PROJECT_DIR / 'covid_phylo_data'
BASE_DIR.mkdir(exist_ok=True)

CACHE_DIR = BASE_DIR / 'cache'

FASTA_DIR = BASE_DIR / 'fasta'

RAW_SEQUENCE_SHELVE_FNAME = 'raw_seqs.shelve'
