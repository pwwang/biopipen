
library(rlang)
library(RUVcorr)

args <- {{envs.ruvcorr_args | r: todot="-"}}
if (!is.null(seed)) { set.seed(seed) }

args$k <- args$k %||% 10
args$size.alpha <- args$size.alpha %||% 2
args$corr.strength <- args$corr.strength %||% 3
args$g <- args$g %||% NULL
args$Sigma.eps <- args$Sigma.eps %||% 1
args$nc <- args$nc %||% (ngenes %/% 4)
args$ne <- args$ne %||% (ngenes %/% 4)
args$intercept <- args$intercept %||% TRUE
args$check <- args$check %||% TRUE
args$n = ngenes
args$m = nsamples

log_info("Running simulation ...")
sim <- do_call(simulateGEdata, args)
attributes(sim) <- c(attributes(sim), c(simulation_tool = "RUVcorr"))
genes <- paste0("Gene", 1:ngenes)
samples <- paste0("Sample", 1:nsamples)

colnames(sim$Truth) <- genes
rownames(sim$Truth) <- samples
colnames(sim$Y) <- genes
rownames(sim$Y) <- samples
colnames(sim$Noise) <- genes
rownames(sim$Noise) <- samples
colnames(sim$Sigma) <- genes
rownames(sim$Sigma) <- genes

log_info("Saving results ...")
saveRDS(sim, file.path(outdir, "sim.rds"))
saveRDS(sim$Truth, file.path(outdir, "Truth.rds"))
