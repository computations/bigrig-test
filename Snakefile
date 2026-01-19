from utils.configs import *
import utils.logs as logs
import csv
import pathlib
import functools
import json
import gzip
import pandas

import warnings


wildcard_constraints:
    tree_iter="[0-9]+",
    model_iter="[0-9]+",


with warnings.catch_warnings(action="ignore"):
    import ete3


workdir: config["prefix"]


BIGRIG_LOG_BASENAME = "bigrig.log"
BIGRIG_JSON_BASENAME = "results.json"

PERF_COMMAND = "perf stat -j --event=duration_time"

for program in config["programs"]:
    cmd = ""
    if program["profile"]:
        cmd += PERF_COMMAND + " "
    cmd += program["binary"]
    program["command"] = cmd

LAGRANGE_PROGRAM_NAMES = [
    program["name"]
    for program in config["programs"]
    if program["type"] == "lagrange-ng"
]

PROGRAM_NAMES = LAGRANGE_PROGRAM_NAMES

node_distance_fields = [
    "clade",
    "software",
    "error",
    "dispersion",
    "extinction",
    "allopatry",
    "sympatry",
    "copy",
    "jump",
    "root-range",
    "tree-iter",
    "bigrig-iter",
    "clade-size",
    "ranges",
    "time",
]

summary_distance_fields = [
    "software",
    "median-error",
    "mean-error",
    "std-error",
    "dispersion",
    "extinction",
    "allopatry",
    "sympatry",
    "copy",
    "jump",
    "root-range",
    "tree-iter",
    "bigrig-iter",
    "ranges",
    "time",
]

bigrig_time_fields = [
    "time",
    "config-time",
    "execution-time",
    "dispersion",
    "extinction",
    "allopatry",
    "sympatry",
    "copy",
    "jump",
    "root-range",
    "ranges",
    "taxa",
]


def expand_iter_list(config_base):
    ret_list = []
    for i in range(len(config_base)):
        iter_config = config_base[i]
        iters = iter_config["iters"] if "iters" in iter_config else 1
        ret_list += [iter_config] * iters
    return ret_list


config["exp_trees"] = expand_iter_list(config["trees"])
config["exp_models"] = expand_iter_list(config["models"])


def make_program_dict():
    programs = {}
    for entry in config["programs"]:
        if entry["name"] in programs:
            raise RuntimeError()
        programs[entry["name"]] = {
            "binary": entry["binary"],
            "profile": entry["profile"],
            "options": entry["options"] if "options" in entry else None,
        }
    return programs

programs = make_program_dict()

rule all:
    input:
        "notebooks/plots.py.ipynb",
        "results_error.html",


rule make_tree:
    output:
        "{tree_iter}/tree.nwk",
    params:
        size=lambda wildcards: config["exp_trees"][int(wildcards["tree_iter"])]["taxa"],
    shell:
        "treegen -u --size  {params.size} > {output}"


rule setup_bigrig_config:
    input:
        tree="{tree_iter}/tree.nwk",
    output:
        config="{tree_iter}/{model_iter}/bigrig/config.yaml",
        model="{tree_iter}/{model_iter}/bigrig/model.json",
    run:
        model = config["exp_models"][int(wildcards.model_iter)]
        bigrig_params = BigrigParameterSet(model)

        with open(output.model, "w") as jsonfile:
            json.dump(bigrig_params.dict, jsonfile)

        bigrig_files = BigrigFiles(
            config=output.config,
            tree=input.tree,
            prefix=os.path.join(os.path.dirname(output.config), "results"),
            binary=programs["bigrig"],
        )

        bigrig_options = BigrigOptions()

        bigrig_config = BigrigConfig(
            params=bigrig_params, files=bigrig_files, options=bigrig_options
        )

        bigrig_config.write_config()


rule run_bigrig:
    input:
        "{tree_iter}/{model_iter}/{program_name}/config.yaml",
    output:
        align="{tree_iter}/{model_iter}/{program_name}/results.phy",
        result="{tree_iter}/{model_iter}/{program_name}/results.json",
    params:
        command=lambda wildcards: [
            program["binary"]
            for program in config["programs"]
            if program["name"] == wildcards.get("program_name")
        ],
    log:
        "{tree_iter}/{model_iter}/{program_name}/bigrig.log",
    shell:
        "{params.command} --config {input} &> {log}"


for result_file in ["results.annotated.nwk", "results.json"]:
  rule:
    input:
      f"{{tree_iter}}/{{model_iter}}/bigrig/{result_file}",
    output:
      f"{{tree_iter}}/{{model_iter}}/bigrig/{result_file}.gz",
    shell:
      "gzip {input}"

rule setup_lagrange_config:
    input:
        tree="{tree_iter}/tree.nwk",
        data="{tree_iter}/{model_iter}/bigrig/results.phy",
    output:
        config="{tree_iter}/{model_iter}/{program_name}/lagrange.conf",
    run:
        lagrange_files = utils.configs.LagrangeNGFiles(
            config=output.config,
            prefix=os.path.join(os.path.dirname(output.config), "analysis"),
            tree=input.tree,
            data=input.data,
            log=os.path.join(os.path.dirname(output.config), "lagrange.log"),
            binary = programs[wildcards.program_name]["binary"]
        )
        lagrange_params = LagrangeNGParams(
          ranges = config["exp_models"][int(wildcards.model_iter)]["ranges"],
        )
        lagrange_options = LagrangeNGOptions()
        lagrange_options.workers = int(programs[wildcards.program_name]["options"]["workers"])
        lagrange_config = LagrangeNGConfig(files = lagrange_files, 
        params = lagrange_params,
        options = lagrange_options
        )
        lagrange_config.write_config()

def get_lagrange_thread_count(wildcards):
  program_name = wildcards.get("program_name")
  if "options" not in programs[program_name]:
    return 1
  options = programs[program_name]["options"]

  workers = int(options["workers"]) if "workers" in options else 1
  threads_per_worker = int(options["workers"]) if "threads_per_worker" in options else 1
  return workers * threads_per_worker

rule run_lagrange:
    input:
        config="{tree_iter}/{model_iter}/{program_name}/lagrange.conf",
    output:
        results="{tree_iter}/{model_iter}/{program_name}/analysis.results.json",
    log:
        "{tree_iter}/{model_iter}/{program_name}/lagrange.log",
    params:
        command=lambda wildcards: [
            program["command"]
            for program in config["programs"]
            if program["name"] == wildcards.get("program_name")
        ],
    threads: get_lagrange_thread_count
    shadow: "full"
    shell:
        "{params.command} {input.config} &> {log}"


for result_file in ["analysis.results.json"]:
  rule:
    input:
      f"{{tree_iter}}/{{model_iter}}/{{program_name}}/{result_file}",
    output:
      f"{{tree_iter}}/{{model_iter}}/{{program_name}}/{result_file}.gz",
    params:
        command=lambda wildcards: [
            program["command"]
            for program in config["programs"]
            if program["name"] == wildcards.get("program_name")
        ],
    shell:
      "gzip {input}"


rule setup_biogeobears_config:
    input:
        tree="{tree_iter}/tree.nwk",
        data="{tree_iter}/{model_iter}/bigrig/results.phy",
    output:
        "{tree_iter}/{model_iter}/biogeobears/script.r",
    run:
        biogebears_config = BioGeoBEARSConfig()
        biogebears_config.filename = output[0]
        biogebears_config.tree = input.tree
        biogebears_config.data = input.data
        biogebears_config.prefix = "biogeobears"
        biogebears_config.write_config()


rule run_biogeobears:
    input:
        tree="{tree_iter}/tree.nwk",
        data="{tree_iter}/{model_iter}/bigrig/results.phy",
        script="{tree_iter}/{model_iter}/biogeobears/script.r",
    output:
        "{tree_iter}/{model_iter}/biogeobears/results.Rdata",
    log:
        "{tree_iter}/{model_iter}/biogeobears/biogeobears.log",
    shell:
        "Rscript {input.script} &> {log}"


rule compute_distances_lagrange:
    input:
        bigrig_json="{tree_iter}/{model_iter}/bigrig/results.json.gz",
        lagrange_json="{tree_iter}/{model_iter}/{software_name}/analysis.results.json.gz",
        lagrange_log="{tree_iter}/{model_iter}/{software_name}/lagrange.log",
    output:
        "{tree_iter}/{model_iter}/{software_name}/distances.csv",
    run:
        clade_map = logs.CladeMap()
        bigrig = logs.BigrigLog(json_log=input.bigrig_json, clade_map=clade_map)
        lagrange = logs.LagrangeNGLog(
            json_log=input.lagrange_json,
            clade_map=clade_map,
            text_log=input.lagrange_log,
        )

        bigrig_params = bigrig.parameters()

        with open(output[0], "w") as outfile:
            writer = csv.DictWriter(outfile, fieldnames=node_distance_fields)
            writer.writeheader()
            distances = logs.compute_node_distance_list(bigrig, lagrange)
            for clade, distance in distances.items():
                clade_list = clade_map.reverse_lookup(clade)
                clade_size = len(clade_list)
                clade_str = "|".join(clade_list)
                row = {
                    "clade": clade_str,
                    "software": wildcards.get("software_name"),
                    "error": distance,
                    "bigrig-iter": wildcards.model_iter,
                    "tree-iter": wildcards.tree_iter,
                    "clade-size": clade_size,
                }
                row = row | bigrig_params
                row = row | {"time": lagrange.time}
                writer.writerow(row)


rule coalece_distances:
    input:
        expand(
            "{{tree_iter}}/{{model_iter}}/{program}/distances.csv",
            program=PROGRAM_NAMES,
        ),
    output:
        "{tree_iter}/{model_iter}/distances.csv",
    run:
        with open(output[0], "w") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=node_distance_fields)
            writer.writeheader()
            for input_file in input:
                reader = csv.DictReader(open(input_file))
                for row in reader:
                    writer.writerow(row)


rule time_bigrig:
    input:
        result_logs=expand(
            ("{tree_iter}/{model_iter}/bigrig/results.json.gz",),
            tree_iter=range(len(config["exp_trees"])),
            model_iter=range(len(config["exp_models"])),
        ),
        bigrig_logs=expand(
            ("{tree_iter}/{model_iter}/bigrig/bigrig.log",),
            tree_iter=range(len(config["exp_trees"])),
            model_iter=range(len(config["exp_models"])),
        ),
    output:
        "bigrig_times.csv",
    run:
        with open(output[0], "w") as csvfile:
            writer = csv.DictWriter(
                csvfile, fieldnames=bigrig_time_fields, extrasaction="ignore"
            )
            writer.writeheader()
            for results_file, log_file in zip(
                input.result_logs, input.bigrig_logs
            ):
                bigrig_log = logs.BigrigLog(
                    json_log=results_file, text_log=log_file
                )
                writer.writerow(bigrig_log.time_csv_row())


rule combine_node_distances:
    input:
        distances=expand(
            "{tree_iter}/{model_iter}/distances.csv",
            tree_iter=range(len(config["exp_trees"])),
            model_iter=range(len(config["exp_models"])),
        ),
    output:
        "node_distances.csv.gz",
    run:
        with gzip.open(output[0], "wt") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=node_distance_fields)
            writer.writeheader()
            for distances in input.distances:
                with open(distances) as infile:
                    reader = csv.DictReader(infile)
                    for row in reader:
                        writer.writerow(row)

rule make_error_plots:
    input:
        node_distances="node_distances.csv.gz",
    log:
        notebook="notebooks/plots.py.ipynb",
    notebook:
        "notebooks/plots.py.ipynb"

rule make_bigrig_plots:
    input:
        node_distances="bigrig_times.csv",
    log:
        notebook="notebooks/bigrig.r.ipynb",
    output:
      "figs/bigrig.times.boxplot.svg",
      "figs/bigrig.times.linreg.svg",
    notebook:
        "notebooks/bigrig.r.ipynb"

rule make_error_plot_html:
    input:
        notebook="notebooks/plots.py.ipynb",
    output:
        "results_error.html",
    shell:
        "jupyter nbconvert --to html --output-dir . --output {output} {input}"

rule make_bigrig_plot_html:
    input:
        notebook="notebooks/bigrig.r.ipynb",
    output:
        "results_bigrig.html",
    shell:
        "jupyter nbconvert --to html --output-dir . --output {output} {input}"

