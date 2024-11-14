from utils.configs import *
import utils.logs as logs
import csv
import pathlib
import functools
import json

import warnings

with warnings.catch_warnings(action="ignore"):
    import ete3


workdir: config["prefix"]


config["tree_count"] = len(config["treeSizes"])
config["range_count"] = len(config["ranges"])
config["rate_params_count"] = len(config["rate-parameters"])
config["clado_params_count"] = len(config["clado-parameters"])

BIGRIG_LOG_BASENAME = "bigrig.log"
BIGRIG_JSON_BASENAME = "results.json"

PERF_COMMAND = "perf stat -j --event=duration_time"

for program in config["programs"]:
    cmd = ""
    if program["profile"]:
        cmd += PERF_COMMAND + " "
    cmd += program["path"]
    program["command"] = cmd

LAGRANGE_PROGRAM_NAMES = [
    program["name"]
    for program in config["programs"]
    if program["type"] == "lagrange-ng"
]

PROGRAM_NAMES = LAGRANGE_PROGRAM_NAMES

distance_fields = [
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

bigrig_time_fields = [
    "time",
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


rule all:
    input:
        "distances.csv",


rule make_tree:
    output:
        "{tree_iter}/tree.nwk",
    params:
        size=lambda wildcards: config["treeSizes"][int(wildcards["tree_iter"])],
    shell:
        "treegen -u --size  {params.size} > {output}"


rule setup_bigrig_config:
    input:
        "{tree_iter}/tree.nwk",
    output:
        "{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/bigrig/config.yaml",
    run:
        bigrig_config = BigrigConfig()
        bigrig_config.filename = output[0]
        bigrig_config.tree = input[0]
        bigrig_config.prefix = os.path.join(os.path.dirname(output[0]), "results")
        bigrig_config.range_count = config["ranges"][int(wildcards.range_iter)]
        bigrig_config.rate_distribution = utils.configs.make_rate_distribution(
            **config["rate-parameters"][int(wildcards.rate_param_iter)]
        )
        bigrig_config.clado_distribution = utils.configs.make_clado_distribution(
            **config["clado-parameters"][int(wildcards.clado_param_iter)]
        )
        bigrig_config.roll_params()
        bigrig_config.write_config()


rule run_bigrig:
    input:
        "{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/bigrig/config.yaml",
    output:
        align="{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/bigrig/results.phy",
        result="{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/bigrig/results.json",
    log:
        "{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/bigrig/bigrig.log",
    shell:
        "~/wrk/hits/bigrig/bin/bigrig --config {input} &> {log}"


rule setup_lagrange_config:
    input:
        tree="{tree_iter}/tree.nwk",
        data="{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/bigrig/results.phy",
    output:
        "{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/{program_name}/lagrange.conf",
    run:
        lagrange_config = LagrangeConfig()
        lagrange_config.filename = output[0]
        lagrange_config.prefix = os.path.join(os.path.dirname(output[0]), "analysis")
        lagrange_config.tree = input.tree
        lagrange_config.data = input.data
        lagrange_config.range_count = config["ranges"][int(wildcards.range_iter)]
        lagrange_config.write_config()


rule run_lagrange:
    input:
        config="{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/{program_name}/lagrange.conf",
    output:
        results="{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/{program_name}/analysis.results.json",
    log:
        "{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/{program_name}/lagrange.log",
    params:
        command=lambda wildcards: [
            program["command"]
            for program in config["programs"]
            if program["name"] == wildcards.get("program_name")
        ],
    shell:
        "{params.command} {input.config} &> {log}"


rule setup_biogeobears_config:
    input:
        tree="{tree_iter}/tree.nwk",
        data="{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/bigrig/results.phy",
    output:
        "{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/biogeobears/script.r",
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
        data="{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/bigrig/results.phy",
        script="{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/biogeobears/script.r",
    output:
        "{tree_iter}/{rate_param_iter}/biogeobears/results.Rdata",
    log:
        "{tree_iter}/{rate_param_iter}/biogeobears/biogeobears.log",
    shell:
        "Rscript {input.script} &> {log}"


rule compute_distances_lagrange:
    input:
        bigrig_json="{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/bigrig/results.json",
        lagrange_json="{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/{software_name}/analysis.results.json",
        lagrange_log="{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/{software_name}/lagrange.log",
    output:
        "{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/{software_name}/distances.csv",
    run:
        clade_map = logs.CladeMap()
        bigrig = logs.BigrigLog(json_log = input.bigrig_json, clade_map = clade_map)
        lagrange = logs.LagrangeNGLog(
            json_log = input.lagrange_json, clade_map = clade_map, text_log =
            input.lagrange_log
        )

        bigrig_params = bigrig.parameters()

        with open(output[0], "w") as outfile:
            writer = csv.DictWriter(outfile, fieldnames=distance_fields)
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
                    "bigrig-iter": wildcards.rate_param_iter,
                    "tree-iter": wildcards.tree_iter,
                    "clade-size": clade_size,
                }
                row = row | bigrig_params
                row = row | {"time": lagrange.time}
                writer.writerow(row)


rule coalece_distances:
    input:
        expand(
            "{{tree_iter}}/{{rate_param_iter}}_{{clado_param_iter}}_{{range_iter}}/{repeat}/{program}/distances.csv",
            program=PROGRAM_NAMES,
            repeat=range(config["repeats"]),
        ),
    output:
        "{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/c_distances.csv",
    run:
        with open(output[0], "w") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=distance_fields)
            writer.writeheader()
            for input_file in input:
                reader = csv.DictReader(open(input_file))
                for row in reader:
                    writer.writerow(row)


rule time_bigrig:
    input:
        result_logs=expand(
            (
                "{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/bigrig/results.json",
            ),
            tree_iter=range(config["tree_count"]),
            rate_param_iter=range(config["rate_params_count"]),
            clado_param_iter=range(config["clado_params_count"]),
            range_iter=range(config["range_count"]),
            repeat=range(config["repeats"]),
        ),
        bigrig_logs=expand(
            (
                "{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/bigrig/bigrig.log",
            ),
            tree_iter=range(config["tree_count"]),
            rate_param_iter=range(config["rate_params_count"]),
            clado_param_iter=range(config["clado_params_count"]),
            range_iter=range(config["range_count"]),
            repeat=range(config["repeats"]),
        ),
        perf_logs=expand(
            (
                "{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/bigrig/perf.json",
            ),
            tree_iter=range(config["tree_count"]),
            rate_param_iter=range(config["rate_params_count"]),
            clado_param_iter=range(config["clado_params_count"]),
            range_iter=range(config["range_count"]),
            repeat=range(config["repeats"]),
        ),
    output:
        "bigrig_times.csv",
    run:
        with open(output[0], "w") as csvfile:
            writer = csv.DictWriter(
                csvfile, fieldnames=bigrig_time_fields, extrasaction="ignore"
            )
            writer.writeheader()
            for results_file, log_file, perf_file in zip(
                input.result_logs, input.bigrig_logs, input.perf_logs
            ):
                bigrig_log = logs.BigrigLog(
                    json_log=results_file, text_log=log_file, perf_log=perf_file
                )

                writer.writerow(bigrig_log.time_csv_row())


rule combine_distances:
    input:
        distances=expand(
            "{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/c_distances.csv",
            tree_iter=range(config["tree_count"]),
            rate_param_iter=range(config["rate_params_count"]),
            clado_param_iter=range(config["clado_params_count"]),
            range_iter=range(config["range_count"]),
        ),
    output:
        "distances.csv",
    run:
        with open(output[0], "w") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=distance_fields)
            writer.writeheader()
            for distances in input.distances:
                with open(distances) as infile:
                    reader = csv.DictReader(infile)
                    for row in reader:
                        writer.writerow(row)
