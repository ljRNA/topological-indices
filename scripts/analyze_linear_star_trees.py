import pandas as pd
from topological_indices import utilities, tree_generation, topological_indices

rows = {}
for tree_type in ['linear', 'star']:
    for N in range(1, 1401):
        if N % 50 == 0:
            print(N)

        if tree_type == 'linear':
            T = tree_generation.linear_tree(N)
        elif tree_type == 'star':
            T = tree_generation.star_tree(N)

        res = topological_indices.calculate_TIs(T)
        res['type'] = tree_type

        # Manual correction to get exact result instead of a very small float
        if tree_type == 'star':
            res['logλ2'] = 0

        rows[f'{tree_type}-{N}'] = res

dfLinearStar = pd.DataFrame.from_dict(rows, orient='index').sort_values(['type', 'N'])
utilities.save_df(dfLinearStar, 'all_trees/dfLinearStar')