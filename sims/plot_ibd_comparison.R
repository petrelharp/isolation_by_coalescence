flat_dir <- "run_013233"
biased_dir <- "run_018268"

read_divs <- function (dir, rep) {
    divs <- matrix(NA, nrow=640, ncol=640)
    divs[lower.tri(divs, diag=TRUE)] <- scan(file.path(dir, sprintf("pop_500000.trees.%d.divs.tsv", rep)), skip=1)
    divs[upper.tri(divs)] <- t(divs)[upper.tri(divs)]
    samps <- read.table(file.path(dir, sprintf("pop_500000.trees.%d.samples.tsv", rep)), header=TRUE)
    dists <- sqrt(outer(samps$x, samps$x, "-")^2 + outer(samps$y, samps$y, "-")^2)
    samps$pop <- as.numeric(cut(samps$x, 4)) + 4 * (as.numeric(cut(samps$y, 4)) - 1)
    sub_divs <- matrix(NA, nrow=16, ncol=16)
    for (i in 1:16) {
        for (j in (i:16)) {
            sub_divs[i,j] <- sub_divs[j,i] <- mean(divs[samps$pop == i, samps$pop == j], na.rm=TRUE)
        }
    }
    sub_xy <- expand.grid(x=2*(1:4), y=2*(1:4))
    sub_dists <- sqrt(outer(sub_xy$x, sub_xy$x, "-")^2 + outer(sub_xy$y, sub_xy$y, "-")^2)
    return(list(divs=divs, samples=samps, dists=dists, sub_divs=sub_divs, sub_dists=sub_dists))
}

flat <- read_divs(flat_dir, 0)
bias <- read_divs(biased_dir, 0)

ut <- upper.tri(flat$divs)

layout((1:2))
plot(flat$dists[ut], flat$divs[ut], pch=20, cex=0.5, col=adjustcolor('black', 0.5))
plot(bias$dists[ut], bias$divs[ut], pch=20, cex=0.5, col=adjustcolor('black', 0.5))
