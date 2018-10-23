flat_dirs <- c("run_013233", "run_032104", "run_023222")

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

all_subs <- list()
for (fdir in c(flat_dirs, "run_018268")) {
    for (repnum in c(0,1)) {
        x <- read_divs(fdir, repnum)
        ut <- upper.tri(x$divs, diag=TRUE)
        sub_ut <- upper.tri(x$sub_divs, diag=TRUE)
        all_subs <- c(all_subs, list(x$sub_divs[sub_ut]))
            png(file.path(fdir, sprintf("ibd.%d.png", repnum)), width=6*300, height=4*300, res=300)
            par(mar=c(4,4,1,1)+.1)
            plot(x$dists[ut], x$divs[ut] * 1e5, pch=20, 
                 xlab='geographic distance',
                 ylab='genetic distance (per 100Kb)',
                 cex=0.2, col=adjustcolor('black', 0.5))
            points(x$sub_dists[sub_ut], x$sub_divs[sub_ut] * 1e5, 
                   col='red', pch=20)
            dev.off()
    }
}

subdist <- x$sub_dists[sub_ut]
subdivs <- do.call(cbind, all_subs)

mean_div <- rowMeans(subdivs)

pdf(file="all_flat_divergences.pdf", width=5, height=5.5, pointsize=10)
    par(mar=c(4,4,1,1)+.1)
    matplot(mean_div * 1e5, subdivs * 1e5, col=rep(seq_along(flat_dirs), each=2), 
            pch=rep(1:2, length(flat_dirs)), asp=1, cex=0.5,
            xlab='mean across replicates (per 100Kb)',
            ylab='per-replicate divergence (per 100Kb)')
    abline(0,1)
    legend("topleft", col=c(1:3,1,1), legend=c(paste("sim", 1:3), paste("rep", 1:2)), pch=c(rep(1,4),2))
dev.off()

# partition the variance
subdir <- rep(flat_dirs, each=2)
rep_means <- do.call(cbind, tapply(1:ncol(subdivs), subdir, function (x) rowMeans(subdivs[,x])))
total <- sd(subdivs)^2
within_reps <- sd(subdivs - rep_means[,rep(1:length(flat_dirs), each=2)])^2 / total
between_reps <- sd(sweep(rep_means[,rep(1:length(flat_dirs), each=2)], 1, mean_div, "-"))^2 / total

cat(sprintf("percent variance explained by simulation replicate: %f\n", between_reps))
cat(sprintf("percent variance explained by sampling: %f\n", within_reps))
cat(sprintf("percent variance explained by residual noise: %f\n", 1 - within_reps - between_reps))


