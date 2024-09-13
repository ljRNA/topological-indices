import argparse
import logging
import numpy as np
import pandas as pd
from topological_indices import utilities, fitting
from topological_indices.topological_indices import dist_names

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fit analytic distributions to topological indices for all trees."
    )
    utilities.logging_parser(parser)

    args = parser.parse_args()
    logging.basicConfig()
    logger = logging.getLogger(__name__)
    logger.setLevel(args.loglevel)
    dfFit = utilities.load_df('random_trees/dfFitRandom-H2-ABC-W-R')

    data = {'ind': [], 'distribution': [], 'param': [], 'alpha': [], 'beta': []}
    for ind, dist_name in dist_names.items():
        df = dfFit[(dfFit['ind'] == ind) & (dfFit['distribution'] == dist_name)]

        for i in range(df['number_of_parameters'].values[0]):
            logger.info("Fitting parameter %s of %s for %s", i, dist_name, ind)
            if ind == 'W' and i == 0:
                alpha, beta = fitting.invfit(df['N'], df[f'param_{i}'])
            elif ind == 'R' and i == 0:
                alpha, beta = fitting.logfit(df['N'], np.abs(df[f'param_{i}']))
            else:
                alpha, beta = fitting.logfit(df['N'], df[f'param_{i}'])

            data['ind'].append(ind)
            data['distribution'].append(dist_name)
            data['param'].append(i)
            data['alpha'].append(alpha)
            data['beta'].append(beta)

    dfFitParameters = pd.DataFrame(data)

    logger.info(dfFitParameters)

    utilities.save_df(dfFitParameters, 'random_trees/dfFitParameters')
