src <- function(file) readRDS(system.file("extdata", file, package="skyscapeR", mustWork=TRUE))
dest <- "scripts/illustrations/input"
fit <- src("ggantija-random-500.rds")
stopifnot(fit$metadata$nsims==500,identical(fit$data$empirical$x,fit$data$null.hyp$x))
write.csv(data.frame(fit$data$empirical, fit$data$null.hyp[-1]),file.path(dest,"fit.csv"),row.names=FALSE)
write.csv(fit$metadata$local.pval,file.path(dest,"clusters.csv"),row.names=FALSE)
