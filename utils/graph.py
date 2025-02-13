import numpy
import scipy.optimize


def popcount(x):
    return bin(x).count("1")


class Edge:
    def __init__(self, f, t):
        self._from = f
        self._to = t


class Problem:
    def __init__(self, dist_diff):
        self._edges = []
        self._states = len(dist_diff)
        self._regions = self._states.bit_length() - 1
        self._b_eq = dist_diff[:-1]
        for i in range(self._states):
            for j in range(self._states):
                if popcount(i ^ j) == 1:
                    self._edges.append(Edge(i, j))

    def make_matrix(self):
        matrix = []
        for state in range(self._states - 1):
            row = []
            for e in self._edges:
                if e._from == state:
                    row.append(-1.0)
                elif e._to == state:
                    row.append(1.0)
                else:
                    row.append(0.0)
            matrix.append(row)
        return numpy.array(matrix)

    def dist(self):
        A_eq = self.make_matrix()
        c = numpy.ones(len(self._edges))
        ret = scipy.optimize.linprog(c, A_eq=A_eq, b_eq=self._b_eq)
        return numpy.sum(ret.x)

    def normalized_dist(self):
        return self.dist() / self._regions
