import networkx as nx
from topological_indices.topological_indices import calculate_TIs


T = nx.Graph(
    [[1, 2], [2, 3], [3, 4], [3, 5], [5, 6], [4, 7], [7, 8], [7, 9], [7, 10], [7, 11]]
)
TIs = calculate_TIs(T)
TIs_correct = {
    "N": 10,
    "MLD": 5,
    "ALD": 2.8363636363636364,
    "NA": 16,
    "λ2": 0.15331604895404052,
    "λN": 6.062608083281499,
    "logλ2": -1.875253808733737,
    "logλN": 1.8021400842748259,
    "M1": 52,
    "M2": 52,
    "W": 156.0,
    "P": 10,
    "ABC": 7.820349451118949,
    "R": 4.744040581781355,
    "SC": 4.507298959743805,
    "H1": 0.7090909090909094,
    "H2": 2.6743801652892607,
    "J": 4.015694083571248,
    "EE": 25.11508677043549,
    "B": 49,
    "I": 2.3685225277282065,
}

for ind, TI in TIs.items():
    assert TI == TIs_correct[ind], f"{ind} not matching"

print("OK.")
