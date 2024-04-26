from utils.configs import *
import utils.logs as logs
import csv
import ete3

workdir: config["prefix"]

rule all:
  input:
    "distances.csv"

rule make_tree:
  output:
    "tree.nwk"
  shell:
    "treegen --size {config[node_count]} > {output}"

rule setup_bigrig_config:
  input:
    "tree.nwk"
  output:
    "bigrig/config.yaml"
  run:
    config = BigrigConfig()
    config.filename = output[0]
    config.tree = input[0]
    config.roll_params()
    config.write_config()

rule make_bigrig_dataset:
  input:
    "bigrig/config.yaml"
  output:
    "bigrig/results.phy", "bigrig/results.json"
  shell:
    "~/wrk/hits/bigrig/bin/bigrig --config {input}"

rule setup_lagrange_config:
  input:
    tree = "tree.nwk", data = "bigrig/results.phy"
  output:
    "lagrange/lagrange.conf"
  run:
    config = LagrangeConfig()
    config.filename = output[0]
    config.prefix = "lagrange/analysis"
    config.tree = input.tree
    config.data = input.data
    config.write_config()

rule run_lagrange:
  input:
    tree = "tree.nwk", data = "bigrig/results.phy", config = "lagrange/lagrange.conf"
  output:
    "lagrange/analysis.results.json"
  shell:
    "~/wrk/hits/lagrange-ng/bin/lagrange-ng {input.config}"

rule setup_biogeobears_config:
  input:
    tree = "tree.nwk", data = "bigrig/results.phy"
  output:
    "biogeobears/script.r"
  run:
    config = BioGeoBEARSConfig()
    config.filename = output[0]
    config.tree = input.tree
    config.data = input.data
    config.prefix = "biogeobears"
    config.write_config()

rule run_biogeobears:
  input:
    tree = "tree.nwk", data = "bigrig/results.phy", script = "biogeobears/script.r"
  output:
    "biogeobears/results.Rdata"
  shell:
    "Rscript {input.script}"

rule compute_distances:
  input:
    bigrig = "bigrig/results.json", lagrange="lagrange/analysis.results.json"
  output:
    "distances.csv"
  run:
    clade_map = logs.CladeMap()
    bigrig = logs.BigrigLog(input.bigrig, clade_map)
    lagrange = logs.LagrangeNGLog(input.lagrange, clade_map)
    with open(output[0] ,'w') as outfile:
      writer = csv.DictWriter(outfile, ["clade", "software", "error"])
      writer.writeheader()
      distances = logs.compute_node_distance_list(bigrig, lagrange)
      for clade, distance in distances.items():
        clade_str = '|'.join(clade_map.reverse_lookup(clade))
        writer.writerow({'clade': clade_str, 'software': 'lagrange', 'error': distance})
