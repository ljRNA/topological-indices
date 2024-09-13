import os
import argparse
import logging
import networkx as nx
from topological_indices import utilities
from topological_indices.config import data_dir
from topological_indices.tree_generation import random_trees

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze random trees.")
    utilities.logging_parser(parser)
    parser.add_argument(
        "-Ns",
        type=list,
        default=[10, 15, 16, 20, 25, 30, 45, 50, 75, 100, 125, 150, 175, 200, 250, 300],
        help="Tree sizes.",
    )
    parser.add_argument("-M", type=int, default=10000, help="Number of trees for each N.")

    args = parser.parse_args()
    logging.basicConfig()
    logger = logging.getLogger(__name__)
    logger.setLevel(args.loglevel)

    Ns = args.Ns
    M = args.M

    for N in Ns:
        filename = data_dir / "random_trees" / f"random_trees_{N}.dat"
        logger.info("%i \t %s", N, filename)

        if not os.path.exists(filename):
            with open(filename, "wb+") as f:
                logger.info("Generating trees ...")
                Ts = random_trees(N, M)

                for T in Ts:
                    f.write(nx.to_sparse6_bytes(T))
        else:
            logger.info("Trees already exist.")
