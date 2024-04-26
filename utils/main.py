#!/usr/bin/env python3

import os
import subprocess
import yaml
import numpy
import shutil

BIGRIG_EXE = os.path.abspath("../bigrig/bin/bigrig")
LAGRANGE_EXE = os.path.abspath("../lagrange-ng/bin/lagrange-ng")


class BigrigConfig:

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
                }))

    def _roll_params(self):
        self._dispersion = 1.0
        self._extinction = 2.0

        self._allopatry = 1.0
        self._sympatry = 1.0
        self._jump = 0.0
        self._copy = 1.0

        self._root_range = "".join(numpy.random.choice(["0", "1"], 5))

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
    def tree(self):
        """The tree property."""
        return self._tree

    @tree.setter
    def tree(self, value):
        self._tree = os.path.abspath(value)

    @property
    def prefix_base(self):
        return "results"

    @property
    def output_format(self):
        return "json"

    @property
    def iter(self):
        """The iter property."""
        return self._iter

    @iter.setter
    def iter(self, value):
        self._iter = value

    @property
    def filename(self):
        """The filename property."""
        return self._filename

    @filename.setter
    def filename(self, value):
        self._filename = os.path.abspath(value)


class BigrigIter:

    def __init__(self, iter, config=None):
        self._iter = iter

    def _make_config(self):
        self._config = BigrigConfig()
        self._config.filename = os.path.join(self.run_prefix, "test.conf")
        self._config._tree = "test.nwk"
        self._config._roll_params()

    @property
    def run_prefix(self):
        """The run_prefix property."""
        return self._run_prefix

    @run_prefix.setter
    def run_prefix(self, value):
        self._run_prefix = os.path.abspath(value)
        if not os.path.exists(self._run_prefix):
            os.mkdir(self._run_prefix)

    @property
    def iter_prefix(self):
        """The iter_prefix property."""
        return os.path.join(self.run_prefix, str(self._iter))

    @property
    def prefix(self):
        """The prfix property."""
        return os.path.join(self.iter_prefix, self._config.prefix_base)

    @property
    def treefile(self):
        tree_basename = os.path.basename(self.config.tree)
        return os.path.join(self.iter_prefix, tree_basename)

    @property
    def datafile(self):
        return os.path.abspath(self.prefix + ".phy")

    def _copy_files(self):
        shutil.copy(self.config.tree, self.treefile)

    def run(self):
        self._make_config()
        self._config.write_config()

        cmd = BIGRIG_EXE + " --config {config} --prefix {prefix}".format(
            config=self._config.filename, prefix=self.prefix)
        cmd = cmd.split()
        subprocess.run(cmd)

        self._copy_files()

    @property
    def config(self):
        return self._config


class LagrangeIter:

    def __init__(self, bigrig_iter):
        self._bigrig = bigrig_iter

    @property
    def config_filename(self):
        """The config_filename property."""
        return os.path.join(self._bigrig.iter_prefix, "lagrange.conf")

    def _config_tree_line(self):
        return "treefile = {}".format(self._bigrig.treefile) + "\n"

    def _config_data_line(self):
        return "datafile = {}".format(self._bigrig.datafile) + "\n"

    def _config_areanames_line(self):
        return "areanames = RA RB RC RD RE" + "\n"

    def _config_states_line(self):
        return "states" + "\n"

    def _config_workers_line(self):
        return "workers = 1" + "\n"

    def write_config(self):
        with open(self.config_filename, 'w') as outfile:
            outfile.write(self._config_tree_line())
            outfile.write(self._config_data_line())
            outfile.write(self._config_areanames_line())
            outfile.write(self._config_states_line())
            outfile.write(self._config_workers_line())

    def run(self):
        self.write_config()

        cmd = LAGRANGE_EXE + " {config}".format(config=self.config_filename)
        cmd = cmd.split()
        subprocess.run(cmd)


iter = BigrigIter(0)
iter.run_prefix = "trial"
iter.run()
lagrange = LagrangeIter(iter)
lagrange.run()
