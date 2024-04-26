import ete3
import os
import pathlib
import numpy


class VectorCollection:

    def __init__(self):
        self._vectors = {}

    def add_distribution_vector(self, node: [ete3.TreeNode],
                                vector: [numpy.ndarray]):
        self._vectors[self.convert_node_to_key(node)] = vector

    @staticmethod
    def convert_node_to_key(node: [ete3.TreeNode]):
        return frozenset([leaf.name for leaf in node])

    @property
    def clades(self):
        """The foo property."""
        return sorted(self._vectors)

    def get_vector_by_clade(self, clade: [frozenset]):
        return self._vectors[clade]

    def __getitem__(self, key: [frozenset]):
        return self._vectors[key]
