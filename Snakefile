from utils.configs import *
import utils.logs as logs
import csv
import ete3
import pathlib


workdir: config["prefix"]
config['tree_count'] = len(config["treeSizes"])


rule all:
    input:
        "distances.csv",


rule make_tree:
    output:
        "{tree_iter}/tree.nwk",
    params:
      size = lambda wildcards: config['treeSizes'][int(wildcards['tree_iter'])]
    shell:
        "treegen --size {params.size} > {output}"


rule setup_bigrig_config:
    input:
        "{tree_iter}/tree.nwk",
    output:
        "{tree_iter}/{bigrig_iter}/bigrig/config.yaml",
    run:
        config = BigrigConfig()
        config.filename = output[0]
        config.tree = input[0]
        config.prefix = os.path.join(os.path.dirname(output[0]), "results")
        config.roll_params()
        config.write_config()


rule make_bigrig_dataset:
    input:
        "{tree_iter}/{bigrig_iter}/bigrig/config.yaml",
    output:
        "{tree_iter}/{bigrig_iter}/bigrig/results.phy",
        "{tree_iter}/{bigrig_iter}/bigrig/results.json",
    shell:
        "~/wrk/hits/bigrig/bin/bigrig --config {input}"


rule setup_lagrange_config:
    input:
        tree="{tree_iter}/tree.nwk",
        data="{tree_iter}/{bigrig_iter}/bigrig/results.phy",
    output:
        "{tree_iter}/{bigrig_iter}/lagrange/lagrange.conf",
    run:
        config = LagrangeConfig()
        config.filename = output[0]
        config.prefix = os.path.join(os.path.dirname(output[0]), "analysis")
        config.tree = input.tree
        config.data = input.data
        config.write_config()


rule run_lagrange:
    input:
        tree="{tree_iter}/tree.nwk",
        data="{tree_iter}/{bigrig_iter}/bigrig/results.phy",
        config="{tree_iter}/{bigrig_iter}/lagrange/lagrange.conf",
    output:
        "{tree_iter}/{bigrig_iter}/lagrange/analysis.results.json",
    shell:
        "~/wrk/hits/lagrange-ng/bin/lagrange-ng {input.config}"


rule setup_biogeobears_config:
    input:
        tree="{tree_iter}/tree.nwk",
        data="{tree_iter}/{bigrig_iter}/bigrig/results.phy",
    output:
        "{tree_iter}/{bigrig_iter}/biogeobears/script.r",
    run:
        config = BioGeoBEARSConfig()
        config.filename = output[0]
        config.tree = input.tree
        config.data = input.data
        config.prefix = "biogeobears"
        config.write_config()


rule run_biogeobears:
    input:
        tree="{tree_iter}/tree.nwk",
        data="{tree_iter}/{bigrig_iter}/bigrig/results.phy",
        script="{tree_iter}/{bigrig_iter}/biogeobears/script.r",
    output:
        "{tree_iter}/{bigrig_iter}/biogeobears/results.Rdata",
    shell:
        "Rscript {input.script}"


rule compute_distances:
    input:
        bigrig="{tree_iter}/{bigrig_iter}/bigrig/results.json",
        lagrange="{tree_iter}/{bigrig_iter}/lagrange/analysis.results.json",
    output:
        "{tree_iter}/{bigrig_iter}/distances.csv",
    run:
        clade_map = logs.CladeMap()
        bigrig = logs.BigrigLog(input.bigrig, clade_map)
        lagrange = logs.LagrangeNGLog(input.lagrange, clade_map)

        bigrig_params = bigrig.parameters()

        with open(output[0], "w") as outfile:
            fields = ['clade', 'software', 'error', 'dispersion', 'extinction',
            'allopatry', 'sympatry', 'copy', 'jump', 'root-range', 'tree-iter',
            'bigrig-iter']
            writer = csv.DictWriter(outfile, fieldnames=fields)
            writer.writeheader()
            distances = logs.compute_node_distance_list(bigrig, lagrange)
            for clade, distance in distances.items():
                clade_str = "|".join(clade_map.reverse_lookup(clade))
                row = {
                  "clade": clade_str, 
                  "software": "lagrange-ng", 
                  "error": distance,
                  "bigrig-iter": wildcards.bigrig_iter,
                  "tree-iter": wildcards.tree_iter,
                  }
                row = row | bigrig_params
                writer.writerow( row)


rule combine_distances:
    input:
        distances=expand(
            "{tree_iter}/{bigrig_iter}/distances.csv",
            tree_iter=range(config['tree_count']),
            bigrig_iter=range(config['parameters']),
        ),
        configs=expand(
            "{tree_iter}/{bigrig_iter}/bigrig/config.yaml",
            tree_iter=range(config['tree_count']),
            bigrig_iter=range(config['parameters']),
        ),
    output:
        "distances.csv",
    run:
        with open(output[0], 'w') as csvfile:
          fields = ['clade', 'software', 'error', 'dispersion', 'extinction',
          'allopatry', 'sympatry', 'copy', 'jump', 'root-range', 'tree-iter',
          'bigrig-iter']
          writer = csv.DictWriter(csvfile, fieldnames=fields)
          writer.writeheader()
          for distances, config_filename in zip(input.distances, input.configs):
            with open(distances) as infile:
              reader = csv.DictReader(infile)
              for row in reader:
                writer.writerow(row)

