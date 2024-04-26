#!/usr/bin/env Rscript

library(ape)
library(optimx)
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)

treefile = "{tree}"
datafile = "{data}"
resultsfile = "{results}"

BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = treefile
BioGeoBEARS_run_object$geogfn = datafile
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE

BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

res = bears_optim_run(BioGeoBEARS_run_object)
save(res, file=resultsfile)
