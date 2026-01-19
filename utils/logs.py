import json
import os
import numpy
import utils.vector as vector
import utils.graph as graph
import re
import gzip

import warnings

with warnings.catch_warnings(action="ignore"):
    import ete3


NANOSECOND_RATIO = 1e-9


class CladeMap:
    def __init__(self):
        self._clades = {}
        self._next_index = 0

    def _add_clade(self, clade):
        self._clades[clade] = self._next_index
        self._next_index += 1

    def __getitem__(self, clade):
        if clade not in self._clades:
            self._add_clade(clade)
        return self._clades[clade]

    @property
    def clades(self):
        return self._clades

    def reverse_lookup(self, id):
        for k, v in self._clades.items():
            if v == id:
                return k


class BaseLog:
    def __init__(self, **kwargs):
        self.logfile = kwargs.get("json_log")
        if kwargs.get("text_log") is not None:
            self.text_log = kwargs.get("text_log")
        self._clade_map = kwargs.get("clade_map")
        if kwargs.get("perf_log") is not None:
            self._perf_log = kwargs.get("perf_log")
        self._setup()

    @property
    def vectors(self):
        if not hasattr(self, "_vectors"):
            self._vectors = vector.VectorCollection(self._clade_map)
        return self._vectors

    @property
    def logfile(self):
        """The logfile property."""
        return self._logfile

    @logfile.setter
    def logfile(self, value):
        if value is not None:
            self._logfile = os.path.abspath(value)
        else:
            self._logfile = None

    @property
    def text_log(self):
        """The text_log property."""
        return self._text_log

    @text_log.setter
    def text_log(self, value):
        self._text_log = os.path.abspath(value)

    @property
    def indices(self):
        """The indices property."""
        return self._indices

    @property
    def taxa(self):
        if not hasattr(self, "_taxa"):
            self._taxa = 0
            for node in self._tree:
                if not node.is_leaf():
                    continue
                self._taxa += 1

        return self._taxa

    def _get_time(self):
        with open(self._perf_log) as perf_log:
            for line in perf_log:
                j = json.loads(line)
                if j["event"] == "duration_time":
                    return float(j["counter-value"]) * NANOSECOND_RATIO

    @property
    def time(self):
        if not hasattr(self, "_time"):
            self._time = self._get_time()
        return self._time


class LagrangeNGLog(BaseLog):
    def _setup(self):
        self._log = json.load(gzip.open(self.logfile, "rb"))
        self._tree = ete3.Tree(
            re.sub(r"\[[^]]*\]", "", self._log["attributes"]["nodes-tree"]),
            format=1,
        )
        self._regions = self._log["attributes"]["regions"]
        self._add_vectors()
        self._time = self._get_execution_time(self.text_log)

    def _add_vectors(self):
        for node in self._tree.traverse("postorder"):
            if node.is_leaf() or node.name == "":
                continue
            node_results = self._get_node_results(node.name)
            self.vectors.add_distribution_vector(
                node, self._make_vector(node_results["states"])
            )

    def _make_vector(self, state_list):
        ret = numpy.zeros(2**self._regions)
        for item in state_list:
            ret[item["distribution"]] = item["ratio"]
        return ret

    def _get_node_results(self, number):
        for item in self._log["node-results"]:
            if str(item["number"]) == number:
                return item

    @staticmethod
    def _get_execution_time(text_log_filename):
        timeline = open(text_log_filename).readlines()[-1]
        return timeline.split()[-1].strip(" s")


class BigrigLog(BaseLog):
    def _setup(self):
        self._log = json.load(gzip.open(self.logfile, "rb"))
        self._tree = ete3.Tree(self._log["tree"] + ";", format=1)
        self._regions = int(self._log["regions"])
        if self._clade_map is not None:
            self._make_distributions()
            self._add_vectors()

    def _make_distributions(self):
        self._distributions = {
            k: int(v, 2) for k, v in self._log["align"].items()
        }

    def _add_vectors(self):
        for node in self._tree.traverse("postorder"):
            if node.is_leaf() or node.name == "":
                continue
            self.vectors.add_distribution_vector(
                node, self._make_vector(self._distributions[node.name])
            )

    def _make_vector(self, distribution):
        ret = numpy.zeros(2**self._regions)
        ret[distribution] = 1.0
        return ret

    def parameters(self):
        p = self._log["periods"][0]
        return {
            "allopatry": p["cladogenesis"]["allopatry"],
            "sympatry": p["cladogenesis"]["sympatry"],
            "copy": p["cladogenesis"]["copy"],
            "jump": p["cladogenesis"]["jump"],
            "dispersion": p["rates"]["dispersion"],
            "extinction": p["rates"]["extinction"],
            "root-range": self._log["align"]["0"],
            "ranges": self.ranges,
        }

    @property
    def time(self):
        """The time property."""
        return self._log["stats"]["time"]

    @property
    def config_time(self):
        """The time property."""
        return self._log["stats"]["config-time"]

    @property
    def execution_time(self):
        """The time property."""
        return self._log["stats"]["execution-time"]

    def time_csv_row(self):
        tmp = self.parameters()
        tmp |= {
            "time": self.time,
            "execution-time": self.execution_time,
            "config-time": self.config_time,
            "taxa": self._log["taxa"],
        }
        return tmp

    @property
    def ranges(self):
        if not hasattr(self, "_ranges"):
            self._ranges = len(self._log["align"]["0"])
        return self._ranges


def DistributionVectorGenerator(log1: BaseLog, log2: BaseLog):
    v1 = log1.vectors
    v2 = log2.vectors
    for c1, c2 in zip(v1.clades, v2.clades):
        if not c1 == c2:
            raise Exception(
                "There was an issue with generating a zipped list of vectors"
            )
        yield (c1, v1[c1], v2[c2])


def compute_node_distance_list(log1: BaseLog, log2: BaseLog):
    distances = {}

    for c, d1, d2 in DistributionVectorGenerator(log1, log2):
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
