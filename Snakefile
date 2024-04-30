from utils.configs import *
import utils.logs as logs
import csv
import ete3
import pathlib


workdir: config["prefix"]
config['tree_count'] = len(config["treeSizes"])
config['range_count'] = len(config['ranges'])

distance_fields = ['clade', 'software', 'error', 'dispersion', 'extinction', 'allopatry', 'sympatry', 'copy', 'jump', 'root-range', 'tree-iter', 'bigrig-iter', 'clade-size', 'ranges']

rule all:
    input:
        "distances.csv",


rule make_tree:
    output:
        "{tree_iter}/tree.nwk",
    params:
      size = lambda wildcards: config['treeSizes'][int(wildcards['tree_iter'])]
    shell:
        "treegen -u --size  {params.size} > {output}"


rule setup_bigrig_config:
    input:
        "{tree_iter}/tree.nwk",
    output:
        "{tree_iter}/{bigrig_iter}_{range_iter}/bigrig/config.yaml",
    run:
        bigrig_config = BigrigConfig()
        bigrig_config.filename = output[0]
        bigrig_config.tree = input[0]
        bigrig_config.prefix = os.path.join(os.path.dirname(output[0]), "results")
        bigrig_config.range_count = config['ranges'][int(wildcards.range_iter)]
        bigrig_config.roll_params()
        bigrig_config.write_config()


rule make_bigrig_dataset:
    input:
        "{tree_iter}/{bigrig_iter}_{range_iter}/bigrig/config.yaml",
    output:
        "{tree_iter}/{bigrig_iter}_{range_iter}/bigrig/results.phy",
        "{tree_iter}/{bigrig_iter}_{range_iter}/bigrig/results.json",
    log:
        "{tree_iter}/{bigrig_iter}_{range_iter}/bigrig/bigrig.log",
    shell:
        "~/wrk/hits/bigrig/bin/bigrig --config {input} &> {log}"


rule setup_lagrange_config:
    input:
        tree="{tree_iter}/tree.nwk",
        data="{tree_iter}/{bigrig_iter}_{range_iter}/bigrig/results.phy",
    output:
        "{tree_iter}/{bigrig_iter}_{range_iter}/lagrange/lagrange.conf",
    run:
        lagrange_config = LagrangeConfig()
        lagrange_config.filename = output[0]
        lagrange_config.prefix = os.path.join(os.path.dirname(output[0]), "analysis")
        lagrange_config.tree = input.tree
        lagrange_config.data = input.data
        lagrange_config.range_count = config['ranges'][int(wildcards.range_iter)]
        lagrange_config.write_config()


rule run_lagrange:
    input:
        config="{tree_iter}/{bigrig_iter}_{range_iter}/lagrange/lagrange.conf",
    output:
        "{tree_iter}/{bigrig_iter}_{range_iter}/lagrange/analysis.results.json",
    log:
        "{tree_iter}/{bigrig_iter}_{range_iter}/lagrange/lagrange.log",
    shell:
        "~/wrk/hits/lagrange-ng/bin/lagrange-ng {input.config} &> {log}"
 

rule setup_biogeobears_config:
    input:
        tree="{tree_iter}/tree.nwk",
        data="{tree_iter}/{bigrig_iter}_{range_iter}/bigrig/results.phy",
    output:
        "{tree_iter}/{bigrig_iter}_{range_iter}/biogeobears/script.r",
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
        data="{tree_iter}/{bigrig_iter}_{range_iter}/bigrig/results.phy",
        script="{tree_iter}/{bigrig_iter}_{range_iter}/biogeobears/script.r",
    output:
        "{tree_iter}/{bigrig_iter}/biogeobears/results.Rdata",
    log:
        "{tree_iter}/{bigrig_iter}/biogeobears/biogeobears.log",
    shell:
        "Rscript {input.script} &> {log}"


rule compute_distances:
    input:
        bigrig="{tree_iter}/{bigrig_iter}_{range_iter}/bigrig/results.json",
        lagrange="{tree_iter}/{bigrig_iter}_{range_iter}/lagrange/analysis.results.json",
    output:
        "{tree_iter}/{bigrig_iter}_{range_iter}/distances.csv",
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
                "bigrig-iter": wildcards.bigrig_iter,
                "tree-iter": wildcards.tree_iter,
                "clade-size": clade_size,
                }
              row = row | bigrig_params
              writer.writerow(row)


rule combine_distances:
    input:
        distances=expand(
            "{tree_iter}/{bigrig_iter}_{range_iter}/distances.csv",
            tree_iter=range(config['tree_count']),
            bigrig_iter=range(config['parameters']),
            range_iter = range(config['range_count']),
        ),
        configs=expand(
            "{tree_iter}/{bigrig_iter}_{range_iter}/bigrig/config.yaml",
            tree_iter=range(config['tree_count']),
            bigrig_iter=range(config['parameters']),
            range_iter = range(config['range_count']),
        ),
    output:
        "distances.csv",
    run:
        with open(output[0], 'w') as csvfile:
          writer = csv.DictWriter(csvfile, fieldnames=distance_fields)
          writer.writeheader()
          for distances, config_filename in zip(input.distances, input.configs):
            with open(distances) as infile:
              reader = csv.DictReader(infile)
              for row in reader:
                writer.writerow(row)

