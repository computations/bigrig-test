import numpy

import warnings

with warnings.catch_warnings(action="ignore"):
    import ete3


class VectorCollection:
    def __init__(self, clade_map):
        self._vectors = {}
        self._clade_map = clade_map

    def add_distribution_vector(
        self, node: [ete3.TreeNode], vector: [numpy.ndarray]
    ):
        self._vectors[self.convert_node_to_key(node)] = vector

    def convert_node_to_key(self, node: [ete3.TreeNode]):
        return self._clade_map[frozenset([leaf.name for leaf in node])]

    @property
    def clades(self):
        """The foo property."""
        return sorted(self._vectors)

    def get_vector_by_clade(self, clade: [frozenset]):
        return self._vectors[clade]

    def __getitem__(self, key: [frozenset]):
        return self._vectors[key]
