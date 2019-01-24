#script for inference on poplar data

require(gene.flow.inference)
require(Matrix)

fname <- "populus"
name <- "Populus"

outpath <- getwd()

#set RNG seed
seed <- sample(1000000000,1)

#read in pairwise genetic distances for individuals
dist <- read.table("populus_genetic_distance.tsv",header=TRUE)
#make into matrix
dist <- as.matrix(dist)

#read in info about individuals
info <- read.table("populus_info.tsv",header=TRUE)

#read in locations of groups
centroids <- read.table("populus_group_centroids.tsv",header=TRUE)

#find total number of individuals
t_ind <- length(dist[,1])

#find number of demes
n <- length(centroids$group)

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

#find number of movement parameters
ng <- length(G_adj@x)

#plot graph structure
gr <- igraph::graph_from_adjacency_matrix(G_adj,mode="directed",weighted=TRUE)
igraph::E(gr)$curved <- TRUE
igraph::V(gr)$color <- terrain.colors(n)
plot(gr,layout=as.matrix(centroids[,1:2]))

#find which individuals are in which groups
indiv <- lapply(1:n,FUN=function(x){which(info$groups==centroids$group[x])})

#find group mean genetic distances and standard errors
#initialize matrices and helping variables
Hs <- matrix(nrow=n,ncol=n) #mean genetic distance matrix (location to location)
Hs_se <- matrix(nrow=n,ncol=n) #matrix for standard error of the mean for values of Hs
dist_ij <- list() #list for subsets of full (genome to genome) genetic distance
k <- 0

#fill in lower triangle of Hs and Hs_se
for(i in 1:n){ #i, j, and k are different here as compared to above
  for(j in 1:i){ 
    k <- k+1
    #find indices (and distances) of both genomes of the individuals in the ith and jth locations
    #subset of matrix for each pair of locations
    #dist_ij[[k]] <- dist[start_ind[i]:(start_ind[i+1]-1),start_ind[j]:(start_ind[j+1]-1)]
    dist_ij[[k]] <- dist[indiv[[i]],indiv[[j]]]
    #remove data on and above diag if this is a self comparison #don't do this for real data
    #if(i==j) dist_ij[[k]][upper.tri(dist_ij[[k]],diag=TRUE)] <- NA
    #means for each pair
    Hs[i,j] <- mean(dist_ij[[k]],na.rm=TRUE)
    #standard error of the mean (#divide by sqrt of min # of ind in either location)
    Hs_se[i,j] <- sd(dist_ij[[k]],na.rm=TRUE)/sqrt(min(length(indiv[[i]]),length(indiv[[j]])))
    #print(length(dist_ij[[k]]))
    #print(length(locsi),locsj)
  }
}

#make it symmetric (replace upper triangle  with lower triangle)
Hs[upper.tri(Hs)] <- t(Hs)[upper.tri(Hs)]
Hs_se[upper.tri(Hs_se)] <- t(Hs_se)[upper.tri(Hs_se)]

#make H into vector
hs <- as.vector(Hs[upper.tri(Hs,diag=TRUE)])
#make standard error into a vector
hs_se <- as.vector(Hs_se[upper.tri(Hs_se,diag=TRUE)])

#do some plotting

#plot H matrix
cairo_pdf(filename=paste0(outpath,"/H_matrix_",fname,".pdf"),width=6,height=6)
  image(Hs,xaxt='n',yaxt='n',main=paste0("Genetic Distance Matrix (H): ",name))
dev.off()

#mcmc iterations
preburn_iter <- 2e6
burn_iter <- 15e6
iter <- 15e6

save(list=ls(),file=paste0(outpath,"/",fname,".RData"))

#run mcmc
time <- system.time(a <- run.mcmc(G_adj_known=TRUE,G_adj=G_adj,g_known=FALSE,const_coal=FALSE,H=Hs,h_se=hs_se,seed=seed,
                          preburn_iter=preburn_iter,burn_iter=burn_iter,iter=iter,noisy_H=FALSE,type="coal"))

save(list=ls(),file=paste0(outpath,"/",fname,".RData"))

cairo_pdf(filename=paste0(outpath,"/posterior_dists_",fname,".pdf"),width=8,height=7)
  boxplot(a$ans$g[10*(1:(iter/10)),order(a$ans$g_med)],outline=FALSE,main=paste0("Posterior Distributions: ",fname),
        names=paste0("g",order(a$ans$g_med)),xlab="Parameter Index",ylab="Parameter Value (rate)",las=2)
dev.off()

cairo_pdf(filename=paste0(outpath,"/grid_",fname,".pdf"),width=14,height=14)
  gr <- igraph::graph_from_adjacency_matrix(G_adj,mode="directed",weighted=TRUE)
  igraph::E(gr)$curved <- TRUE
  igraph::V(gr)$color <- terrain.colors(n)
  plot(gr,layout=as.matrix(centroids[,1:2]),edge.width=(a$ans$g_med+.2),
       edge.label=paste0("g",1:ng,"=",round(a$ans$g_med*1000)/1000),main=paste0("Graph Structure"))
dev.off()
  
#png(filename="g_trace%02d.png",width=800,height=600)
#  for(i in 1:ng){
#    plot(a$ans$g[(1:10000)*(iter/10000),i],type='l')
#  }
#dev.off()

#png(filename="gam_trace%02d.png",width=800,height=600)
#  for(i in 1:n){
#    plot(a$ans$gam[(1:10000)*(iter/10000),i],type='l')
#  }
#dev.off()

#set up for runnning it for longer

#save startpoints for continuing runs
init <- list()
init$g <- a$ans$g[iter,]
init$gam <- a$ans$gam[iter,]
init$sdG <- a$ans$sdG
init$sdgam <- a$ans$sdgam


time2 <- system.time(a2 <- run.mcmc(G_adj_known=TRUE,G_adj=G_adj,g_known=FALSE,const_coal=FALSE,H=Hs,h_se=hs_se,seed=seed,
                                   preburn_iter=preburn_iter,burn_iter=burn_iter,iter=iter,noisy_H=FALSE,
                                   continue=TRUE,g_init=init$g,gam_init=init$gam,sdG_init=init$sdG,sdgam_init=init$sdgam,
                                   type="coal"))

save(list=ls(),file=paste0(outpath,"/",fname,".RData"))

gr <- igraph::graph_from_adjacency_matrix(G_adj,mode="directed",weighted=TRUE)
igraph::E(gr)$curved <- TRUE
igraph::V(gr)$color <- terrain.colors(n)
plot(gr,layout=as.matrix(centroids[,1:2]),edge.width=(a2$ans$g_med+.2),
     edge.label=paste0("g",1:ng,"=",round(a2$ans$g_med*1000)/1000),main=paste0("Graph Structure"))

#png(filename="g2_trace%02d.png",width=800,height=600)
#  for(i in 1:ng){
#    plot(a2$ans$g[(1:10000)*(iter/10000),i],type='l')
#  }
#dev.off()

g <- rbind(a$ans$g,a2$ans$g)
gam <- rbind(a$ans$gam,a2$ans$gam)
lllh <- c(a$ans$lllh,a2$ans$lllh)

g_med <- matrixStats::colMedians(g)
gam_med <- matrixStats::colMedians(gam)

#save medians
save(file="g_med.txt",g_med)
save(file="gam_med.txt",gam_med)

cairo_pdf(filename=paste0(outpath,"/grid_",fname,".pdf"),width=14,height=14)
  gr <- igraph::graph_from_adjacency_matrix(G_adj,mode="directed",weighted=TRUE)
  igraph::E(gr)$curved <- TRUE
  igraph::V(gr)$color <- terrain.colors(n)
  plot(gr,layout=as.matrix(centroids[,1:2]),edge.width=(g_med+.2),
       edge.label=paste0("g",1:ng,"=",round(g_med*1000)/1000),main=paste0("Graph Structure"))
dev.off()

cairo_pdf(filename=paste0(outpath,"/posterior_dists_g_",fname,".pdf"),width=8,height=7)
  boxplot(g[100*(1:(30e6/100)),order(g_med)],outline=FALSE,main=paste0("Posterior Distributions (g): ",name),
        names=paste0("g",order(g_med)),xlab="Parameter Index",ylab="Parameter Value (rate)",las=2)
dev.off()

cairo_pdf(filename=paste0(outpath,"/posterior_dists_gam_",fname,".pdf"),width=8,height=7)
  boxplot(gam[100*(1:(30e6/100)),],outline=FALSE,main=paste0("Posterior Distributions (gamma): ",name),
        names=paste0("gam",1:n),xlab="Parameter Index",ylab="Parameter Value (rate)",las=2)
dev.off()
