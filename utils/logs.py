import json
import os
import numpy
import utils.vector as vector
import utils.graph as graph
import ete3


class BaseLog:

    def __init__(self, logfile):
        self.logfile = logfile
        self._vectors = vector.VectorCollection()
        self._setup()

    @property
    def logfile(self):
        """The logfile property."""
        return self._logfile

    @logfile.setter
    def logfile(self, value):
        self._logfile = os.path.abspath(value)

    @property
    def indices(self):
        """The indices property."""
        return self._indices

    @property
    def vectors(self):
        return self._vectors

    @property
    def taxa(self):
        if getattr(self, "_taxa") is None:
            self._taxa = 0
            for node in self._tree:
                if not node.is_leaf():
                    continue
                self._taxa += 1

        return self._taxa


class LagrangeNGLog(BaseLog):

    def _setup(self):
        self._log = json.load(open(self.logfile))
        self._tree = ete3.Tree(self._log["attributes"]["nodes-tree"], format=1)
        self._regions = self._log["attributes"]["regions"]
        self._add_vectors()

    def _add_vectors(self):
        for node in self._tree.traverse("postorder"):
            if node.is_leaf() or node.name == '':
                continue
            node_results = self._get_node_results(node.name)
            self._vectors.add_distribution_vector(
                node, self._make_vector(node_results["states"]))

    def _make_vector(self, state_list):
        ret = numpy.zeros(2**self._regions)
        for item in state_list:
            ret[item["distribution"]] = item["ratio"]
        return ret

    def _get_node_results(self, number):
        for item in self._log["node-results"]:
            if str(item["number"]) == number:
                return item


class BigrigLog(BaseLog):

    def _setup(self):
        self._log = json.load(open(self.logfile))
        self._tree = ete3.Tree(self._log["tree"] + ";", format=1)
        self._regions = int(self._log["regions"])
        self._make_distributions()
        self._add_vectors()

    def _make_distributions(self):
        self._distributions = {
            k: int(v, 2)
            for k, v in self._log["align"].items()
        }

    def _add_vectors(self):
        for node in self._tree.traverse("postorder"):
            if node.is_leaf() or node.name == '':
                continue
            print(self._distributions)
            self._vectors.add_distribution_vector(
                node, self._make_vector(self._distributions[node.name]))

    def _make_vector(self, distribution):
        ret = numpy.zeros(2**self._regions)
        ret[distribution] = 1.0
        return ret


def DistributionVectorGenerator(log1: BaseLog, log2: BaseLog):
    v1 = log1.vectors
    v2 = log2.vectors
    for c1, c2 in zip(v1.clades, v2.clades):
        if not c1 == c2:
            raise Exception(
                "There was an issue with generating a zipped list of vectors")
        yield (c1, v1[c1], v2[c2])


def compute_node_distance_list(log1: BaseLog, log2: BaseLog):
    distances = {}

    for c, d1, d2 in DistributionVectorGenerator(log1, log2):
        print(d1)
        d = d1 - d2
        distances[c] = graph.Problem(d).normalized_dist()

    return distances


def compute_tree_distance(log1: BaseLog, log2: BaseLog):
    total_dist = 0.0

    for _, d in compute_node_distance_list(log1, log2).items():
        total_dist += d

    return total_dist


def compute_norm_tree_distance(log1: BaseLog, log2: BaseLog):
    if not log1.taxa == log2.taxa:
        raise Exception("Logs do not have the same taxa size")

    return compute_tree_distance(log1, log2) / (log1.taxa - 1)