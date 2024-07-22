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
    "regions",
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
    "regions",
    "taxa"
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


rule make_bigrig_dataset:
    input:
        "{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/bigrig/config.yaml",
    output:
        "{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/bigrig/results.phy",
        "{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/bigrig/results.json",
    log:
        "{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/bigrig/bigrig.log",
    shell:
        "~/wrk/hits/bigrig/bin/bigrig --config {input} &> {log}"


rule setup_lagrange_config:
    input:
        tree="{tree_iter}/tree.nwk",
        data="{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/bigrig/results.phy",
    output:
        "{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/lagrange/lagrange.conf",
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
        config="{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/lagrange/lagrange.conf",
    output:
        "{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/lagrange/analysis.results.json",
    log:
        "{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/lagrange/lagrange.log",
    shell:
        "~/wrk/hits/lagrange-ng/bin/lagrange-ng {input.config} &> {log}"


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
        bigrig="{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/bigrig/results.json",
        lagrange="{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/lagrange/analysis.results.json",
    output:
        "{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/lagrange/distances.csv",
    run:
        clade_map = logs.CladeMap()
        bigrig = logs.BigrigLog(input.bigrig, clade_map)
        lagrange = logs.LagrangeNGLog(input.lagrange, clade_map)

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
                    "software": "lagrange-ng",
                    "error": distance,
                    "bigrig-iter": wildcards.rate_param_iter,
                    "tree-iter": wildcards.tree_iter,
                    "clade-size": clade_size,
                }
                row = row | bigrig_params
                writer.writerow(row)


rule coalece_distances:
    input:
        expand(
            "{{tree_iter}}/{{rate_param_iter}}_{{clado_param_iter}}_{{range_iter}}/{repeat}/{program}/distances.csv",
            program=config["programs"],
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
        logs=expand(
            "{tree_iter}/{rate_param_iter}_{clado_param_iter}_{range_iter}/{repeat}/bigrig/results.json",
            tree_iter=range(config["tree_count"]),
            rate_param_iter=range(config["rate_params_count"]),
            clado_param_iter=range(config["clado_params_count"]),
            range_iter=range(config["range_count"]),
            repeat=range(config["repeats"])
        )
    output:
        "bigrig_times.csv"
    run:
        with open(output[0], "w") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=bigrig_time_fields,
                extrasaction='ignore')
            writer.writeheader()
            for input_file in input.logs:
                js = json.load(open(input_file))

                js['time'] = js['stats']['time']
                js |= js['periods'][0]['cladogenesis']
                js |= js['periods'][0]['rates']

                writer.writerow(js)


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
