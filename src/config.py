
from pathlib import Path

PROJECT_DIR = Path('.').resolve().parent

BASE_DIR = PROJECT_DIR / 'covid_phylo_data'
BASE_DIR.mkdir(exist_ok=True)

CACHE_DIR = BASE_DIR / 'cache'
CACHE_DIR.mkdir(exist_ok=True)

FASTA_DIR = BASE_DIR / 'fasta'
FASTA_DIR.mkdir(exist_ok=True)

MEDIA_DIR = BASE_DIR / 'media'
MEDIA_DIR.mkdir(exist_ok=True)

TREE_DIR = BASE_DIR / 'tree'
TREE_DIR.mkdir(exist_ok=True)

RAW_SEQUENCE_SHELVE_FNAME = 'raw_seqs.shelve'
