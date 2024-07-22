#!/usr/bin/env python3

import os
import numpy
import yaml
from dataclasses import dataclass

import utils.util

SRC_PATH = os.path.dirname(os.path.abspath(__file__))


@dataclass
class StaticRateDistribution:
    dispersion: float
    extinction: float


@dataclass
class UniformRateDistribution:
    start: float
    end: float


@dataclass
class GammaRateDistribution:
    shape: float
    scale: float


@dataclass
class LogNormalRateDistribution:
    mu: float
    sigma: float


@dataclass
class InverseGaussianRateDistribution:
    mu: float
    l: float


RateDistribution = StaticRateDistribution | UniformRateDistribution \
    | GammaRateDistribution | LogNormalRateDistribution \
    | InverseGaussianRateDistribution


def make_rate_distribution(type, **kwargs):
    match type:
        case "Static":
            return StaticRateDistribution(kwargs['dispersion'],
                                          kwargs['extinction'])
        case "Uniform":
            return UniformRateDistribution(kwargs['a'], kwargs['b'])
        case "Gamma":
            return GammaRateDistribution(kwargs['k'], kwargs['theta'])
        case "LogNormal":
            return LogNormalRateDistribution(kwargs['mu'], kwargs['sigma'])
        case "InverseGaussian":
            return InverseGaussianRateDistribution(
                kwargs['mu'], kwargs['lambda'])


@dataclass
class StaticCladoDistribution:
    allopatry: float
    sympatry: float
    copy: float
    jump: float


@dataclass
class DirichletCladoDistribution:
    allopatry: float
    sympatry: float
    copy: float
    jump: float


CladoDistribution = StaticCladoDistribution | DirichletCladoDistribution


def make_clado_distribution(type, **kwargs):
    match type:
        case "Static":
            return StaticCladoDistribution(
                kwargs['allopatry'],
                kwargs['sympatry'],
                kwargs['copy'],
                kwargs['jump'])
        case "Dirichlet":
            return DirichletCladoDistribution(
                kwargs['allopatry'],
                kwargs['sympatry'],
                kwargs['copy'],
                kwargs['jump'])


class BaseConfig:

    @property
    def tree(self):
        """The treefile property."""
        return self._tree

    @tree.setter
    def tree(self, value):
        self._tree = os.path.abspath(value)

    @property
    def filename(self):
        """The filename property."""
        return self._filename

    @filename.setter
    def filename(self, value):
        self._filename = os.path.abspath(value)

    @property
    def data(self):
        """The data property."""
        return self._data

    @data.setter
    def data(self, value):
        self._data = os.path.abspath(value)

    @property
    def region_count(self):
        """The region_count property."""
        return self._region_count

    @region_count.setter
    def region_count(self, value):
        self._region_count = value

    @property
    def prefix(self):
        """The prefix property."""
        return self._prefix

    @prefix.setter
    def prefix(self, value):
        self._prefix = os.path.abspath(value)

    @property
    def range_count(self):
        """The range_count property."""
        return self._range_count

    @range_count.setter
    def range_count(self, value):
        self._range_count = int(value)


class BigrigConfig(BaseConfig):

    def __init__(self):
        pass

    def write_config(self):
        with open(self.filename, 'w') as outfile:
            outfile.write(
                yaml.dump({
                    "rates": self.rates,
                    "cladogenesis": self.cladogenesis,
                    "root-range": self.root_range,
                    "tree": self.tree,
                    "output-format": self.output_format,
                    "prefix": self.prefix,
                }))

    def roll_params(self):
        self.roll_rate_params()
        self.roll_clado_params()
        self.roll_root_range()

    def roll_rate_params(self):
        match self.rate_distribution:
            case StaticRateDistribution(dis, ext):
                def rate_dist(): return (dis, ext)
            case UniformRateDistribution(a, b):
                def rate_dist(): return (float(f) for f in
                                         numpy.random.uniform(a, b, 2))
            case GammaRateDistribution(k, theta):
                def rate_dist(): return (float(f) for
                                         f in numpy.random.gamma(k, theta, 2))
            case LogNormalRateDistribution(mu, sigma):
                def rate_dist(): return (
                    float(f) for f in numpy.random.lognormal(mu, sigma, 2))
            case InverseGaussianRateDistribution(mu, l):
                def rate_dist(): return (float(f)
                                         for f in numpy.random.wald(mu, l, 2))

        self._dispersion, self._extinction = rate_dist()

    def roll_clado_params(self):
        match self.clado_distribution:
            case StaticCladoDistribution(allopatry, sympatry, copy, jump):
                def clado_dist(): return (allopatry, sympatry, copy, jump)
            case DirichletCladoDistribution(allopatry, sympatry, copy, jump):
                def clado_dist(): 
                    if jump == 0.0:
                        tmp = [
                            float(f) for f in numpy.random.dirichlet(
                                (allopatry,
                                 sympatry,
                                 copy), 1).transpose()]
                        tmp.append(0.0)
                        return tuple(tmp)
                    return (
                        float(f) for f in numpy.random.dirichlet(
                            (allopatry,
                             sympatry,
                             copy, jump), 1).transpose())

        self._allopatry, self._sympatry, self._jump, self._copy = clado_dist()

    def roll_root_range(self):
        self._root_range = ""
        while "1" not in self._root_range:
            self._root_range = "".join(
                numpy.random.choice(["0", "1"], self.range_count))

    @property
    def rates(self):
        if hasattr(self, "_rates"):
            return self._rates
        self._rates = {
            "dispersion": self.dispersion,
            "extinction": self.extinction
        }
        return self._rates

    @property
    def cladogenesis(self):
        if hasattr(self, "_cladogenesis"):
            return self._cladogenesis
        self._cladogenesis = {
            "allopatry": self.allopatry,
            "sympatry": self.sympatry,
            "jump": self.jump,
            "copy": self.copy,
        }
        return self._cladogenesis

    @property
    def rate_distribution(self):
        """The rate_distribution property."""
        return self._rate_distribution

    @rate_distribution.setter
    def rate_distribution(self, value: RateDistribution):
        self._rate_distribution = value

    @property
    def clado_distribution(self):
        """The rate_distribution property."""
        return self._clado_distribution

    @clado_distribution.setter
    def clado_distribution(self, value: CladoDistribution):
        self._clado_distribution = value

    @property
    def dispersion(self):
        return self._dispersion

    @property
    def extinction(self):
        return self._extinction

    @property
    def allopatry(self):
        return self._allopatry

    @property
    def sympatry(self):
        return self._sympatry

    @property
    def jump(self):
        return self._jump

    @property
    def copy(self):
        return self._copy

    @property
    def root_range(self):
        """The root_range property."""
        return self._root_range

    @property
    def output_format(self):
        return "json"


class LagrangeConfig(BaseConfig):

    def __init__(self):
        pass

    def _config_tree_line(self):
        return "treefile = {}".format(self.tree) + "\n"

    def _config_data_line(self):
        return "datafile = {}".format(self.data) + "\n"

    def _config_areanames_line(self):
        area_list = (
            n for _, n in zip(range(self.range_count),
                              utils.util.base26_generator(self.range_count)))
        return "areanames = {}".format(' '.join(area_list)) + "\n"

    def _config_states_line(self):
        return "states" + "\n"

    def _config_workers_line(self):
        return "workers = 1" + "\n"

    def _config_prefix_line(self):
        return "prefix = {}".format(self.prefix) + "\n"

    def write_config(self):
        with open(self.filename, 'w') as outfile:
            outfile.write(self._config_tree_line())
            outfile.write(self._config_data_line())
            outfile.write(self._config_areanames_line())
            outfile.write(self._config_states_line())
            outfile.write(self._config_workers_line())
            outfile.write(self._config_prefix_line())


class BioGeoBEARSConfig(BaseConfig):

    def __init__(self):
        pass

    @property
    def base_script(self):
        return os.path.join(SRC_PATH, "biogeobears.r")

    def _copy_files(self):
        with open(self.filename, 'w') as outfile:
            script = open(self.base_script).read()
            outfile.write(
                script.format(tree=self.tree,
                              data=self.data,
                              results=self.results))

    @property
    def results(self):
        """The results property."""
        return os.path.join(self.prefix, "results.Rdata")

    def write_config(self):
        self._copy_files()
