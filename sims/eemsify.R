# convert our data to EEMS format

usage <- "eemsify.R DISTANCE_FILE SAMPLE_FILE OUT_FILEBASE"

args <- commandArgs(TRUE)
cat("args: ", args, "\n")
if (length(args) != 3) { stop(usage) }
DISTFILE <- args[1]
SAMPLEFILE <- args[2]
OUTBASE <- args[3]

DIFFFILE <- paste0(OUTBASE, ".diffs")
COORDFILE <- paste0(OUTBASE, ".coord")
OUTERFILE <- paste0(OUTBASE, ".outer")

x <- scan(DISTFILE, skip=1)
z <- read.table(SAMPLEFILE, header=TRUE, stringsAsFactors=FALSE)
z$individual <- ordered(z$individual, levels=unique(z$individual))
stopifnot(all(diff(as.numeric(z$individual)) >= 0))
stopifnot(all(diff(as.numeric(z$individual)) <= 1))
nr <- nrow(z)
d <- matrix(NA, nrow=nr, ncol=nr)
d[lower.tri(d, diag=TRUE)] <- x
d[upper.tri(d)] <- t(d)[upper.tri(d)]
dd <- do.call(rbind, tapply(1:nrow(d), z$individual, function (u) rowMeans(d[,u], na.rm=TRUE)))
ddd <- do.call(cbind, tapply(1:nrow(d), z$individual, function (u) rowMeans(dd[,u], na.rm=TRUE)))

write.table(ddd / mean(ddd), file=DIFFFILE, col.names=FALSE, row.names=FALSE)

ii <- match(levels(z$individual), z$individual)
write.table(z[ii,c("x","y")], file=COORDFILE, col.names=FALSE, row.names=FALSE)

x0 <- min(z$x)/2
x1 <- max(z$x)+x0
y0 <- min(z$y)/2
y1 <- max(z$y)+y0
outers <- matrix(c(x0, y0,
                   x0, y1,
                   x1, y1,
                   x1, y0,
                   x0, y0), ncol=2, byrow=TRUE)
write.table(outers, file=OUTERFILE, col.names=FALSE, row.names=FALSE)
