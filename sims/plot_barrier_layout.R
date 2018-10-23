locs <- read.table("run_032600/pop_500000.trees.locs.tsv", skip=1, header=FALSE)
names(locs) <- c("x", "y", "z")
locs <- subset(locs, locs$x > 0.0)
samples <- read.table("run_032600/pop_500000.trees.0.samples.tsv", header=TRUE)


W <- 8.0
aspectRatio <- 5/3
wh <- W * sqrt(c(aspectRatio, 1/aspectRatio))
LOWER_Y <- (W / sqrt(aspectRatio) ) / 3
UPPER_Y <- 2 * (W / sqrt(aspectRatio) ) / 3
LOWER_LEFT <- 1 * (W * sqrt(aspectRatio) ) / 5
LOWER_RIGHT <- 3 * (W * sqrt(aspectRatio) ) / 5
UPPER_LEFT <- 2 * (W * sqrt(aspectRatio) ) / 5
UPPER_RIGHT <- 4 * (W * sqrt(aspectRatio) ) / 5
BARRIER_HALFWIDTH <- 0.2 * 3 / 2


pdf(file="run_032600/sample_locations_pretty.pdf", width=5, height=3.5, pointsize=10)
par(mar=c(2, 2, 2, 0)+.2, mgp=c(1,1,0))
    plot(locs$x, locs$y, asp=1, type='n',
         xlab='Position (x)', ylab='Position (y)',
         main='Individual locations',
         xaxt='n', yaxt='n')
    polygon(x=c(0, 0, wh[1], wh[1]), 
            y=c(0, wh[2], wh[2], 0),
            col=adjustcolor("yellow", 0.25))
    polygon(x=c(LOWER_LEFT, LOWER_LEFT, LOWER_RIGHT, LOWER_RIGHT),
            y=LOWER_Y + c(-1, +1, +1, -1) * BARRIER_HALFWIDTH,
            col=adjustcolor("red", 0.5),
            border=NA)
    polygon(x=c(UPPER_LEFT, UPPER_LEFT, UPPER_RIGHT, UPPER_RIGHT),
            y=UPPER_Y + c(-1, +1, +1, -1) * BARRIER_HALFWIDTH,
            col=adjustcolor("red", 0.5),
            border=NA)
    segments(x0=rep(0,4), x1=rep(wh[1], 4), y0=seq(0, wh[2], length.out=4))
    segments(x0=seq(0, wh[1], length.out=6), y0=rep(0, 6), y1=rep(wh[2], 6))
    points(locs$x, locs$y, pch=20, cex=0.5, col=adjustcolor("black", 0.5))
    points(samples$x, samples$y, pch=20, cex=2)
dev.off()
