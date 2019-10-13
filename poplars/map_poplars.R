library(sp)
library(maps)
library(colorspace)
library(Matrix)
library(rgeos)
library(raster)
library(maptools)
library(rgdal)



info <- read.table("populus_info.tsv",header=TRUE)
centroids <- read.table("populus_group_centroids.tsv",header=TRUE)

pop_coords <- sp::SpatialPoints(as.matrix(info[,c("Longitude", "Latitude")]),
                            proj4string=sp::CRS("+proj=longlat"))

#make adjacency matrix
G_adj <- Matrix(0,nrow=9,ncol=9,sparse=TRUE)
G_adj <- as(G_adj,'dgCMatrix')
#make edges manually
G_adj[1,c(6,7)] <- 1
G_adj[2,c(4,5,7,8)] <- 1
G_adj[3,c(4,5,8,9)] <- 1
G_adj[4,c(2,3,8)] <- 1
G_adj[5,c(2,3,6)] <- 1
G_adj[6,c(1,5,7)] <- 1
G_adj[7,c(1,2,6)] <- 1
G_adj[8,c(2,3,4,9)] <- 1
G_adj[9,c(3,8)] <- 1

g_med <- scan("g_med.tsv")
gam_med <- scan("gam_med.tsv")


gr <- igraph::graph_from_adjacency_matrix(G_adj,mode="directed",weighted=TRUE)
igraph::E(gr)$curved <- TRUE
igraph::V(gr)$color <- cols <- adjustcolor(RColorBrewer::brewer.pal(9, "Set3"), 0.95)
edge_labels <- sprintf("%0.2f", g_med)

outpath <- "."


the.proj4 <- proj4string(pop_coords)
xbox <- as(extent(pop_coords),"SpatialPolygons")
crs(xbox) <- CRS( the.proj4 )
if (!file.exists("natural_earth_cleaned.RData")) {
    land <- spTransform( readOGR(file.path("natural_earth"),"ne_10m_land"), crs(the.proj4) )
    tmp <- land
    # this takes FOREVER
    tmp@polygons <- lapply(land@polygons, checkPolygonsHoles)
    tmp2 <- gBuffer(tmp, byid=TRUE, width=0)
    bigbox <- SpatialPolygons(
                list(Polygons(list(
                       Polygon(cbind(c(-155,-105)[c(1,1,2,2)],
                                     c(35,75)[c(1,2,2,1)]))), ID='box')),
                proj4string=CRS(the.proj4))
    land <- crop(tmp2, bigbox)
    save(land, file="natural_earth_cleaned.RData")
} else {
    load(file.path(.thisdir, "natural_earth_cleaned.RData"))
}
water <- rgeos::gDifference(bigbox, crop(land, bigbox))


source("igraph_plot.R")
cairo_pdf(filename=paste0(outpath,"/tmp_plot_map.pdf"),width=14,height=14)
  # set up map
    plot(pop_coords, pch=21, cex=2, 
         col=adjustcolor(c("blue", "red")[as.numeric(info$Species)], 0.75),
         bg=cols[as.numeric(info$groups)], lwd=3)
    # points(centroids[,1:2], pch=21, col=adjustcolor(cols[1:nlevels(info$groups)], 0.25), cex=5)
    maps::map(add=TRUE, col=grey(0.6))
  ip <- my_plot_igraph(gr,layout=as.matrix(centroids[,1:2]),
       edge.width=pmin(15, 3 * g_med + .1),
       edge.label=edge_labels,
       main=paste0("Graph Structure"), 
       vertex.size=200, vertex.label.cex=2,
       edge.arrow.size=2.5, edge.label.cex=1.9,
       edge.color=rgb(colorRamp(c("#DD8000","gray","#9999FF"))(tanh(g_med/2)),max=255),
       edge.label.color="#000000",
       add=TRUE, rescale=FALSE)
dev.off()

# manual adjustment to stop overlap
el_x <- ip$lc.x + c(
      0.0, 0.0, 0.0, 0.0, 0.3,-0.4, 0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,-0.3, 0.0,-0.4, 0.0, 0.0, 0.0, 0.0, 0.0)

el_y <- ip$lc.y + c(
      0.0, 0.0, 0.0, 0.0,-0.4, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0,-0.6, 0.0, 0.0, 0.0, 0.0, 0.0)

# now for REALS
cairo_pdf(filename=paste0(outpath,"/poplar_map_cover.pdf"),width=14,height=14)
    plot(pop_coords, pch=21, cex=3, 
         col=adjustcolor(c("blue", "red")[as.numeric(info$Species)], 0.75),
         bg=cols[as.numeric(info$groups)], lwd=3)
    # points(centroids[,1:2], pch=21, col=adjustcolor(cols[1:nlevels(info$groups)], 0.25), cex=5)
    # maps::map(add=TRUE, col=grey(0.6))
    plot(water, col=adjustcolor("blue", 0.25), add=TRUE)
    plot(gr,layout=as.matrix(centroids[,1:2]),
       edge.width=pmin(15, 3 * g_med + .1),
       vertex.size=200, labels=rep("x", 9), vertex.label.cex=0,
       edge.arrow.size=2.5, 
       edge.label.cex=0,
       edge.color=rgb(colorRamp(c("#DD8000","gray","#9999FF"))(tanh(g_med/2)),max=255),
       add=TRUE, rescale=FALSE)
  # text(el_x, el_y, labels=(edge_labels), cex=2,
  #      col=adjustcolor("black", 0.9))
dev.off()
