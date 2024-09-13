import warnings
import logging
import argparse
from time import time
import numpy as np
import pandas as pd
import scipy.stats as stats
from topological_indices import utilities


def _fit_distribution(data, dist_name, bins="auto"):
    logger = logging.getLogger(__name__)

    y, x = np.histogram(data, bins=bins, density=True)
    x = (x + np.roll(x, -1))[:-1] / 2.0

    distribution = getattr(stats, dist_name)

    try:
        # Ignore warnings from data that can't be fit
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")

            params = distribution.fit(data, method="MLE")

            # Calculate fitted PDF and error with fit in distribution
            pdf = distribution.pdf(x, *params)
            SSE = np.sum(np.power(y - pdf, 2.0))

            pdf = distribution.pdf(data, *params)
            log_likelihood = np.sum(np.log(pdf))

            BIC = len(params) * np.log(len(data)) - 2 * log_likelihood

            return params  # , SSE, BIC

    except Exception as e:
        logger.warning("%s", e)
        return []


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fit analytic distributions to topological indices for all trees."
    )
    utilities.logging_parser(parser)

    args = parser.parse_args()
    logging.basicConfig()
    logger = logging.getLogger(__name__)
    logger.setLevel(args.loglevel)

    dfTrees = utilities.load_df("random_trees/dfRandom-H2-ABC-W-R")
    Ns = sorted(dfTrees["N"].unique())

    rows = []
    inds = ["H2", "ABC", "W", "R"]
    dist_names = ["moyal", "lognorm", "pearson3", "pearson3"]

    for N in Ns:
        logger.info("Fitting N = %i ...", N)

        t = time()

        data = dfTrees[dfTrees["N"] == N]

        for ind, dist_name in zip(inds, dist_names):
            params = _fit_distribution(data[ind].values.T.tolist(), dist_name)

            rows.append(
                {
                    "N": N,
                    "ind": ind,
                    "distribution": dist_name,
                    "parameters": params,
                    "number_of_parameters": len(params),
                }
            )

    dfFit = pd.DataFrame(rows)
    dfFit = dfFit.combine_first(dfFit.parameters.apply(pd.Series).add_prefix('param_'))

    logger.info(dfFit)

    utilities.save_df(dfFit, "random_trees/dfFitRandom-H2-ABC-W-R")
