#!/usr/bin/env python3

import os
import shutil
import subprocess

import numpy
import yaml

import utils.util

SRC_PATH = os.path.dirname(os.path.abspath(__file__))


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

    def load(self, filename):
        with open(filename) as infile:
            y = yaml.load(infile, Loader=yaml.Loader)
            self._dispersion = y['rates']['dispersion']
            self._extinction = y['rates']['extinction']

            self._allopatry = y['cladogenesis']['allopatry']
            self._sympatry = y['cladogenesis']['sympatry']
            self._copy = y['cladogenesis']['copy']
            self._jump = y['cladogenesis']['jump']

            self._prefix = y['prefix']

            self._tree = y['tree']
            self._root_range = y['root-range']

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
        self._dispersion = numpy.random.uniform(0.0, 1.0)
        self._extinction = numpy.random.uniform(0.0, 1.0)

        self._allopatry = 1.0
        self._sympatry = 1.0
        self._jump = 0.0
        self._copy = 1.0

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
