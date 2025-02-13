#!/usr/bin/env python3

import os
import numpy
import yaml
import pathlib
from dataclasses import dataclass

import utils.util

SRC_PATH = os.path.dirname(os.path.abspath(__file__))


@dataclass
class DistributionBase:
    def __call__(self):
        if not hasattr(self, "_val"):
            self._val = self._get()
        return self._val

    def lock(self, p):
        self._val = p


@dataclass
class StaticDistribution(DistributionBase):
    param: float

    def _get(self):
        return self.param


@dataclass
class UniformDistribution(DistributionBase):
    start: float
    end: float

    def _get(self):
        return numpy.random.uniform(self.start, self.end)


@dataclass
class GammaRateDistribution(DistributionBase):
    shape: float
    scale: float

    def _get(self):
        return numpy.random.gamma(self.shape, self.scale)


@dataclass
class LogNormalRateDistribution(DistributionBase):
    mu: float
    sigma: float

    def _get(self):
        return numpy.random.lognormal(self.mu, self.sigma)


@dataclass
class InverseGaussianRateDistribution(DistributionBase):
    mu: float
    lam: float

    def _get(self):
        return numpy.random.wald(self.mu, self.lam)


@dataclass
class BetaDistribution(DistributionBase):
    alpha: float
    beta: float

    def _get(self):
        return numpy.random.beta(self.alpha, self.beta)


@dataclass
class ClampedBetaDistribution(BetaDistribution):
    max: float

    def _get(self):
        while True:
            if (val := numpy.random.beta(self.alpha, self.beta)) <= self.max:
                return val


@dataclass
class RootRangeDistribution(DistributionBase):
    regions: int

    def _get(self):
        while "1" not in (
            tmp := "".join(
                [str(i) for i in numpy.random.randint([2] * self.regions)]
            )
        ):
            pass
        return tmp

    def to_str(self):
        return self._val


Distribution = (
    StaticDistribution
    | UniformDistribution
    | GammaRateDistribution
    | LogNormalRateDistribution
    | InverseGaussianRateDistribution
    | BetaDistribution
    | ClampedBetaDistribution
)


def make_distribution(type, **kwargs) -> Distribution:
    match type:
        case "Static":
            return StaticDistribution(kwargs["value"])
        case "Uniform":
            return UniformDistribution(kwargs["start"], kwargs["end"])
        case "Gamma":
            return GammaRateDistribution(kwargs["k"], kwargs["theta"])
        case "LogNormal":
            return LogNormalRateDistribution(kwargs["mu"], kwargs["sigma"])
        case "InverseGaussian":
            return InverseGaussianRateDistribution(
                kwargs["mu"], kwargs["lambda"]
            )
        case "Beta":
            return BetaDistribution(kwargs["alpha"], kwargs["beta"])
        case "ClampedBeta":
            return ClampedBetaDistribution(
                kwargs["alpha"], kwargs["beta"], kwargs["max"]
            )


@dataclass
class RateParamterSet:
    dispersion: Distribution
    extinction: Distribution

    @property
    def dict(self):
        return {
            "dispersion": self.dispersion(),
            "extinction": self.extinction(),
        }

    def load(self, p):
        self.dispersion.lock(p["dispersion"])
        self.extinction.lock(p["extinction"])

    def __init__(self, model_config):
        rates_config = model_config["rate-parameters"]
        self.dispersion = make_distribution(**rates_config["dispersion"])
        self.extinction = make_distribution(**rates_config["extinction"])


@dataclass()
class CladoParameterSet:
    allopatry: Distribution
    sympatry: Distribution
    copy: Distribution
    jump: Distribution

    @property
    def dict(self):
        return {
            "allopatry": self.allopatry(),
            "sympatry": self.sympatry(),
            "copy": self.copy(),
            "jump": self.jump(),
        }

    def load(self, p):
        self.allopatry.lock(p["allopatry"])
        self.sympatry.lock(p["sympatry"])
        self.copy.lock(p["copy"])
        self.jump.lock(p["jump"])

    def __init__(self, model_config):
        clado_params = model_config["clado-parameters"]
        dist_type = clado_params['type']
        if dist_type == "Dirichlet":
            self.allopatry = GammaRateDistribution(
                shape=clado_params["allopatry"], scale=1.0
            )
            self.sympatry = GammaRateDistribution(
                shape=clado_params["sympatry"], scale=1.0
            )
            self.copy = GammaRateDistribution(
                shape=clado_params["copy"], scale=1.0
            )
            self.jump = GammaRateDistribution(
                shape=clado_params["jump"], scale=1.0
            )
        if dist_type == "Static":
            self.allopatry = StaticDistribution(clado_params["allopatry"])
            self.sympatry = StaticDistribution(clado_params["sympatry"])
            self.copy = StaticDistribution(clado_params["copy"])
            self.jump = StaticDistribution(clado_params["jump"])


@dataclass(init=False)
class BigrigParameterSet:
    rates: RateParamterSet
    clado: CladoParameterSet
    root_range: RootRangeDistribution

    @property
    def dict(self):
        return {
            "rates": self.rates.dict,
            "clado": self.clado.dict,
            "root_range": self.root_range(),
        }

    def load(self, p):
        print(p)
        self.rates.load(p["rates"])
        self.clado.load(p["clado"])
        self.root_range.lock(p["root_range"])

    def __init__(self, model_config):
        self.rates = RateParamterSet(model_config)
        self.clado = CladoParameterSet(model_config)
        self.root_range = RootRangeDistribution(regions=model_config["ranges"])


@dataclass
class BigrigOptions:
    output_format: str = "json"


@dataclass
class BigrigFiles:
    config: pathlib.Path
    tree: pathlib.Path
    prefix: pathlib.Path
    binary: pathlib.Path


@dataclass
class BigrigConfig:
    params: BigrigParameterSet
    files: BigrigFiles
    options: BigrigOptions

    def write_config(self):
        with open(self.files.config, "w") as outfile:
            outfile.write(
                yaml.dump(
                    {
                        "rates": self.params.rates.dict,
                        "cladogenesis": self.params.clado.dict,
                        "root-range": self.params.root_range.to_str(),
                        "tree": self.files.tree,
                        "output-format": self.options.output_format,
                        "prefix": self.files.prefix,
                    }
                )
            )

    def command(self):
        return self.binary + " " + f"--config {self.files.config}"


@dataclass
class LagrangeNGFiles:
    config: pathlib.Path
    tree: pathlib.Path
    data: pathlib.Path
    log: pathlib.Path
    prefix: pathlib.Path
    binary: pathlib.Path


@dataclass
class LagrangeNGParams:
    ranges: int


@dataclass
class LagrangeNGOptions:
    expm_mode: str = "adaptive"
    workers: int = 1
    threads_per_worker: int = 1
    output_type: str = "json"
    opt_method: str = "bfgs"


@dataclass
class LagrangeNGConfig:
    files: LagrangeNGFiles
    params: LagrangeNGParams
    options: LagrangeNGOptions

    def _config_tree_line(self):
        return f"treefile = {self.files.tree}\n"

    def _config_data_line(self):
        return f"datafile = {self.files.data}\n"

    def _config_areanames_line(self):
        area_list = (
            n
            for _, n in zip(
                range(self.params.ranges),
                utils.util.base26_generator(self.params.ranges),
            )
        )
        return "areanames = {}".format(" ".join(area_list)) + "\n"

    def _config_states_line(self):
        return "states\n"

    def _config_workers_line(self):
        return f"workers = {self.options.workers}\n"

    def _config_prefix_line(self):
        return f"prefix = {self.files.prefix}\n"

    def write_config(self):
        with open(self.files.config, "w") as outfile:
            outfile.write(self._config_tree_line())
            outfile.write(self._config_data_line())
            outfile.write(self._config_areanames_line())
            outfile.write(self._config_states_line())
            outfile.write(self._config_workers_line())
            outfile.write(self._config_prefix_line())

    def command(self):
        return f"{self.files.binary} {self.files.config}"


class BioGeoBEARSConfig:
    def __init__(self):
        pass

    @property
    def base_script(self):
        return os.path.join(SRC_PATH, "biogeobears.r")

    def _copy_files(self):
        with open(self.filename, "w") as outfile:
            script = open(self.base_script).read()
            outfile.write(
                script.format(
                    tree=self.tree, data=self.data, results=self.results
                )
            )

    @property
    def results(self):
        """The results property."""
        return os.path.join(self.prefix, "results.Rdata")

    def write_config(self):
        self._copy_files()
