import gzip
import pickle
import logging
from pathlib import Path
import pandas as pd
from topological_indices.config import data_dir


def logging_parser(parser):
    parser.add_argument(
        "-d",
        "--debug",
        help="Print lots of debugging statements",
        action="store_const",
        dest="loglevel",
        const=logging.DEBUG,
        default=logging.WARNING,
    )
    parser.add_argument(
        "-v",
        "--verbose",
        help="Be verbose",
        action="store_const",
        dest="loglevel",
        const=logging.INFO,
    )


def get_RNA_sequence(seqid):
    with open(data_dir / "rna" / "sequences" / f"{seqid}.fasta", "r") as f:
        lines = f.readlines()

    return lines[1] if lines[0][0] == ">" else lines[0]


def get_RNA_structures(seqid):
    with open(data_dir / "rna" / "structures" / f"{seqid}.db", "r") as f:
        dbs = f.readlines()

    dbs = [db.strip() for db in dbs]

    return dbs


def get_RNA_trees(seqid):
    with gzip.open(data_dir / "rna" / "trees" / f"{seqid}.gz", "rb") as f:
        Ts = pickle.load(f)

    return Ts


def save_RNA_trees(seqid, Ts):
    with gzip.open(data_dir / "rna" / "trees" / f"{seqid}.tree.gz", "wb") as f:
        pickle.dump(Ts, f)

    return Ts


def save_RNA_structures(seqid, dbs):
    with open(data_dir / "rna" / "structures" / f"{seqid}.db", "w") as f:
        f.writelines("\n".join(dbs))


def load_df(name):
    file = data_dir / Path(name + ".pkl")
    file_compressed = data_dir / Path(name + ".pkl.bz2")

    if file_compressed.exists():
        file = file_compressed

    return pd.read_pickle(file)


def save_df(df, name, compressed=None):
    if compressed:
        df.to_pickle(data_dir / Path(name + ".pkl.bz2"))
    else:
        df.to_pickle(data_dir / Path(name + ".pkl"))
