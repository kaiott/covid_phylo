
from pathlib import Path

PROJECT_DIR = Path('.').resolve().parent

BASE_DIR = PROJECT_DIR / 'covid_phylo_data'
BASE_DIR.mkdir(exist_ok=True)

CACHE_DIR = BASE_DIR / 'cache'

FASTA_DIR = BASE_DIR / 'fasta'

MEDIA_DIR = BASE_DIR / 'media'

TREE_DIR = BASE_DIR / 'tree'

RAW_SEQUENCE_SHELVE_FNAME = 'raw_seqs.shelve'
