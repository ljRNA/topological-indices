import multiprocessing
import argparse
import logging
import glob
import pandas as pd
import networkx as nx
from topological_indices import utilities
from topological_indices.topological_indices import calculate_TIs
from topological_indices.config import data_dir


def _calculate_TIs(Ts):
    return calculate_TIs(Ts, ["N", "H2", "ABC", "W", "R"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze random trees.")
    utilities.logging_parser(parser)
    parser.add_argument(
        "-cpu", type=int, default=2, help="Number of CPUs for multiprocessing Pool."
    )

    args = parser.parse_args()
    logging.basicConfig()
    logger = logging.getLogger(__name__)
    logger.setLevel(args.loglevel)

    pool = multiprocessing.Pool(args.cpu)

    rows = []
    for file in glob.glob(str(data_dir / "random_trees" / "*.dat")):
        logger.info("Reading  %s ...", file)
        with open(file, "r") as f:
            lines = f.readlines()

        Ts = [
            nx.from_sparse6_bytes(line.strip().encode()) for i, line in enumerate(lines)
        ]

        logger.info("Calculating TIs ...")
        _rows = list(pool.map(_calculate_TIs, Ts))

        for i in range(len(_rows)):
            _rows[i]["i"] = i

        rows += _rows

    df = pd.DataFrame(rows)
    df = df.sort_values(['N', 'i'])
    
    logger.debug(df)

    utilities.save_df(df, "random_trees/dfRandom-H2-ABC-W-R")
