# Run from the Quarto project with build 2.0.9004 installed.
stopifnot(as.character(packageVersion("skyscapeR"))=="2.0.9004")
src <- function(file) readRDS(system.file("extdata", file, package="skyscapeR", mustWork=TRUE))
dest <- "scripts/illustrations/input"
grid <- src("milky-way-layers.rds")
stopifnot(identical(dim(grid),c(1440L,481L)))
writeBin(as.vector(grid),file.path(dest,"classes.bin"))
fit <- src("ggantija-random-500.rds")
stopifnot(fit$metadata$nsims==500,identical(fit$data$empirical$x,fit$data$null.hyp$x))
write.csv(data.frame(fit$data$empirical, fit$data$null.hyp[-1]),file.path(dest,"fit.csv"),row.names=FALSE)
write.csv(fit$metadata$local.pval,file.path(dest,"clusters.csv"),row.names=FALSE)
