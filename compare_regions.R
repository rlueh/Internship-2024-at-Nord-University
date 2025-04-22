source("functions.R")
source("create_plots.R")
require("Biostrings")
require("ape")
require("httr")
require("jsonlite")

organisms <- list.files(path="analytics2")

scaffold_folders <- lapply(organisms, function(organism){
    folders <- list.files(path=paste0("analytics2/",organism), pattern="scaffold", full.names=TRUE)
    list("folders"=folders,"organism_name"=organism)
})

##remove Hypophtalmichthys nobilis, somethings causes NaN for k3 & k4, !!fix later!!
scaffold_folders <- scaffold_folders[-17]


data <- list()
for(i in 1:length(scaffold_folders)){
    print(paste0("Iteration: ", i))
    files <- lapply(scaffold_folders[[i]]$folders, function(folders){
        print(paste0("Folder:",folders))
        regions <- read.csv(paste0(folders,"/regions.csv"), check.names=FALSE)
        k1 <- read.csv(paste0(folders,"/kmer_spectra_k1.csv"), check.names=FALSE)
        k1 <- k1[,-1]
        k2 <- read.csv(paste0(folders,"/kmer_spectra_k2.csv"), check.names=FALSE)
        k2 <- k2[,-1]
        k3 <- read.csv(paste0(folders,"/kmer_spectra_k3.csv"), check.names=FALSE)
        k3 <- k3[,-1]
        k4 <- read.csv(paste0(folders,"/kmer_spectra_k4.csv"), check.names=FALSE)
        k4 <- k4[,-1]
        list("organism"=scaffold_folders[[i]]$organism_name,"regions"=regions,"spectra_k1"=k1,"spectra_k2"=k2,"spectra_k3"=k3,"spectra_k4"=k4)
    })
##    data$stats[i] <- files
##    data$organism[i] <- scaffold_folders[[i]]$organism_name
    data[[i]] <- files
}

get.sp <- function(data){
    sapply(data, function(x){
        sp <- x[[1]]$organism
    })
}

sp_orders <- taxonomy.get(get.sp(data))
sp <- unique(k1.mgd$sp)

sp_no_akrab <- readLines("sp_lost_akrab.txt")
sp_no_akrab <- strsplit(sp_no_akrab, ",")[[1]]



## spectrum defines the kmer spectra to extract
## and is the name of the lst.
merge.kmer.counts <- function(spectrum, data=data, orders=orders, freq=TRUE){
    tmp <- do.call(rbind, lapply(data, function(x){
        ## x is for a given species
        ##browser()
        sp <- x[[1]]$organism
        sp <- gsub("_"," ",sp)
         if(!is.na(unlist(orders)[sp])){
             sp_order <- unlist(orders)[sp]
             sp_order <- unname(sp_order)
         } else{
             sp_order <- NA
         }
         if(sp_order%in%sp_no_akrab){
             akrab <- FALSE
         } else{
             akrab <- TRUE
         }
         if(is.na(sp_order)){
             akrab <- NA
         }
        nr <- sapply(x, function(y){ nrow( y[[spectrum]] )})
        do.call(rbind, lapply(1:length(x[nr > 0]), function(i){
            ##browser()
            y <- x[[i]]
            ## y is a scaffold
            if(freq == TRUE){
                y[[ spectrum ]] <- y[[ spectrum ]] / rowSums(y[[spectrum]])
            }
            data.frame( sp, sp_order, akrab,  y[["regions"]][,c('size', 'ent')], sc=paste0("sc_", i), y[[ spectrum ]])
        }))
    }))
    tmp[ !is.na(tmp$size), ]
}


k1.mgd <- merge.kmer.counts("spectra_k1", data=data, orders=sp_orders)
k2.mgd <- merge.kmer.counts("spectra_k2", data=data, orders=sp_orders)
k3.mgd <- merge.kmer.counts("spectra_k3", data=data, orders=sp_orders)
k4.mgd <- merge.kmer.counts("spectra_k4", data=data, orders=sp_orders)

k2.mgd.sp <- do.call(rbind, tapply(1:nrow(k2.mgd), k2.mgd$sp, function(i){
    colMeans(k2.mgd[i,7:ncol(k2.mgd)])
}))

k3.mgd.sp <- do.call(rbind, tapply(1:nrow(k3.mgd), k3.mgd$sp, function(i){
    colMeans(k3.mgd[i,7:ncol(k3.mgd)])
}))

k4.mgd.sp <- do.call(rbind, tapply(1:nrow(k4.mgd), k4.mgd$sp, function(i){
    colMeans( k4.mgd[i,7:ncol(k4.mgd)] )
}))

order <- vector()
akrab <- vector()
for(i in 1:nrow(k2.mgd.sp)){
   order[i] <- unlist(sp_orders)[rownames(k2.mgd.sp)[i]]
   if(!is.na(order[i])){
       akrab[i] <- !order[i]%in%sp_no_akrab
   } else{
       akrab[i] <- NA
   }
}


k4.mgd.sp <- data.frame(k4.mgd.sp)
k4.mgd.sp <- cbind(k4.mgd.sp, order, akrab)

k2.mgd.sp <- data.frame(k2.mgd.sp)
k2.mgd.sp <- cbind(k2.mgd.sp, order, akrab)

k3.mgd.sp <- data.frame(k3.mgd.sp)
k3.mgd.sp <- cbind(k3.mgd.sp, order, akrab)

######


sp.cols <- hsv(0.75 * 1:length(sp) / length(sp), 0.8, 0.8, 0.2)
names(sp.cols) <- sp

with(k1.mgd, plot( C+G, ent, cex=0.5, col=sp.cols[sp] ))

barplot(with(k1.mgd, tapply(1:length(sp), sp, length)), las=2)

with(k1.mgd, plot(k4.pca$x[,1], C+G, col=sp.cols[sp]))
with(k1.mgd, plot(k4.pca$x[,2], C+G, col=sp.cols[sp]))

comp.lm <- apply(k4.pca$x, 2, function(x){
    lm(x ~ k1.mgd$sp)
})

k2.mgd.sorted <- k2.mgd[k2.mgd$size<512,]
k2.pca.sc.sorted <- prcomp( k2.mgd.sorted[, 7:ncol(k2.mgd.sorted)], scale=TRUE )
pdf("analysis_results/pca_sc_k2_less_values.pdf", width=12, height=8)
plot(k2.pca.sc.sorted)
with(k2.pca.sc.sorted, plot(x[,1], x[,2], col=ifelse(is.na(k2.mgd.sorted[rownames(x),"akrab"]),adjustcolor("grey", alpha.f=0.2),ifelse(k2.mgd.sorted[rownames(x),"akrab"]==TRUE,adjustcolor("green", alpha.f=0.2),adjustcolor("red",alpha.f=0.2))), main="DIM1 to DIM2"), cex=0.5)
with(k2.pca.sc.sorted, plot(x[,2], x[,3], col=ifelse(is.na(k2.mgd.sorted[rownames(x),"akrab"]),adjustcolor("grey", alpha.f=0.2),ifelse(k2.mgd.sorted[rownames(x),"akrab"]==TRUE,adjustcolor("green", alpha.f=0.2),adjustcolor("red",alpha.f=0.2))), main="DIM2 to DIM3"), cex=0.5)
with(k2.pca.sc.sorted, plot(x[,4], x[,5], col=ifelse(is.na(k2.mgd.sorted[rownames(x),"akrab"]),adjustcolor("grey", alpha.f=0.2),ifelse(k2.mgd.sorted[rownames(x),"akrab"]==TRUE,adjustcolor("green", alpha.f=0.2),adjustcolor("red",alpha.f=0.2))), main="DIM4 to DIM5"), cex=0.5)
dev.off()

#######
set.seed("222")
k2.pca.sc.kmeans <- kmeans(k2.pca.sc$x[,], centers=10, iter.max=10)
kmeans_clusters <- tapply(1:length(k2.pca.sc.kmeans$cluster), k2.pca.sc.kmeans$cluster, eval)

cluster_data <- lapply(kmeans_clusters, function(cluster){
        k2.mgd[cluster,]
})

pdf("analysis_results/compare_clusters.pdf", height=8, width=14)
for(i in 1:length(cluster_data)){
    par(mfrow=c(2,2))
    cluster <- cluster_data[[i]]
    cluster_mean <- colMeans(cluster[,7:ncol(cluster)])
    barplot(cluster_mean, main=paste0("Cluster ", i))
    akrab_counts <- c("true"=length(which(cluster$akrab==TRUE)),"false"=length(which(cluster$akrab==FALSE)))    
    cluster_orders <- unique(cluster$sp_order)
    order_amounts <- unlist(lapply(cluster_orders, function(x){
        length(which(cluster$sp_order==x))
    }))
    barplot(order_amounts, names.arg=cluster_orders, las=2, main="Amount of regions per species")
    barplot(akrab_counts, names.arg=c("has aKrab", "lost aKrab"), col=c("green", "red"), main="Amount of regions from species who have or lost aKrab")
    barplot(cluster$size)
}
dev.off()
    
cluster_mean <- do.call(rbind, tapply(1:nrow(cluster), cluster, function(i){
    colMeans(cluster[i,7:ncol(cluster)])
}))

cluster_mean <- colMeans(cluster[,7:ncol(cluster)])

barplot(cluster_mean)

barplot(unlist(cluster[,7:ncol(cluster)]))


points(k2.pca.sc$x[unlist(tmp.vals["1"]),1:2], col="red")


plot(k2.pca.sc$x[,1:2], col=as.numeric(k2.pca.sc.kmeans$cluster))
plot(k2.pca.sc$x[,3:4], col=as.numeric(k2.pca.sc.kmeans$cluster))
plot(k2.pca.sc$x[,4:5], col=as.numeric(k2.pca.sc.kmeans$cluster))

tmp.ind1 <- as.numeric(names(tmp$cluster[tmp$cluster==1]))

points(k2.pca.sc$x[tmp.ind1,1:2], col="red")


test.val <- k2.mgd[1:100,]
test.pca <- prcomp(test.val[,7:ncol(test.val)], scale=TRUE)

test.kmean <- kmeans(test.pca$x[,1:2], centers=4, iter.max=10)
plot(test.pca$x[,1:2], col=as.numeric(test.kmean$cluster))
test.ind1 <- as.numeric(names(test.kmean$cluster[test.kmean$cluster==1]))
test.ind2 <- as.numeric(names(test.kmean$cluster[test.kmean$cluster==2]))

points(test.pca$x[test.ind1,1:2], col="brown")
points(test.pca$x[test.ind2,1:2], col="purple")



plot(k2.pca.sc$x[,1:2], col=ifelse(k2.pca.sc$x[,1:2]%in%tmp.cluster6,"green","red"))


########
k2.mgd.dst <- dist(k2.mgd[,7:ncol(k2.mgd)])

test <- unique(k2.mgd$sp)


k2.sc <- sapply(unique(k2.mgd$sp), function(x){
    k2.sp <- k2.mgd[k2.mgd$sp==x,]
    k2.sc.means <- sapply(unique(k2.sp$sc), function(k){
        colMeans(k2.sp[k2.sp$sc==k,7:ncol(k2.sp)])
    })
    data.frame("sp"=x,k2.sc.means)
})

k2.sc <- do.call(cbind, k2.sc)
k2.sc <- t(k2.sc)



k2.sc.dst <- dist(k2.sc)

k2.sc.dst <- dist(k2.sc[[1]][,2:ncol(k2.sc[[1]])])

k2.sc.nj <- nj(k2.sc.dst)
plot(k2.sc.nj)

                
########


k2.pca.sc <- prcomp( k2.mgd[, 7:ncol(k2.mgd)], scale=FALSE)
pdf("analysis_results/pca_sc_k2.pdf", width=12, height=8)
plot(k2.pca.sc)
with(k2.pca.sc, plot(x[,1], x[,2], col=ifelse(is.na(k2.mgd[rownames(x),"akrab"]),adjustcolor("grey", alpha.f=0.2),ifelse(k2.mgd[rownames(x),"akrab"]==TRUE,adjustcolor("green", alpha.f=0.2),adjustcolor("red",alpha.f=0.2))), main="DIM1 to DIM2"), cex=0.5)
with(k2.pca.sc, plot(x[,2], x[,3], col=ifelse(is.na(k2.mgd[rownames(x),"akrab"]),adjustcolor("grey", alpha.f=0.2),ifelse(k2.mgd[rownames(x),"akrab"]==TRUE,adjustcolor("green", alpha.f=0.2),adjustcolor("red",alpha.f=0.2))), main="DIM2 to DIM3"), cex=0.5)
with(k2.pca.sc, plot(x[,4], x[,5], col=ifelse(is.na(k2.mgd[rownames(x),"akrab"]),adjustcolor("grey", alpha.f=0.2),ifelse(k2.mgd[rownames(x),"akrab"]==TRUE,adjustcolor("green", alpha.f=0.2),adjustcolor("red",alpha.f=0.2))), main="DIM4 to DIM5"), cex=0.5)
dev.off()

k3.mgd.sorted <- k3.mgd[k3.mgd$size<512,]
k3.pca.sc.sorted <- prcomp( k3.mgd.sorted[, 7:ncol(k3.mgd.sorted)], scale=TRUE )
pdf("analysis_results/pca_sc_k3_less_values.pdf", width=12, height=8)
plot(k3.pca.sc.sorted)
with(k3.pca.sc.sorted, plot(x[,1], x[,2], col=ifelse(is.na(k3.mgd.sorted[rownames(x),"akrab"]),adjustcolor("grey", alpha.f=0.2),ifelse(k3.mgd.sorted[rownames(x),"akrab"]==TRUE,adjustcolor("green", alpha.f=0.2),adjustcolor("red",alpha.f=0.2))), main="DIM1 to DIM2"), cex=0.5)
with(k3.pca.sc.sorted, plot(x[,2], x[,3], col=ifelse(is.na(k3.mgd.sorted[rownames(x),"akrab"]),adjustcolor("grey", alpha.f=0.2),ifelse(k3.mgd.sorted[rownames(x),"akrab"]==TRUE,adjustcolor("green", alpha.f=0.2),adjustcolor("red",alpha.f=0.2))), main="DIM2 to DIM3"), cex=0.5)
with(k3.pca.sc.sorted, plot(x[,4], x[,5], col=ifelse(is.na(k3.mgd.sorted[rownames(x),"akrab"]),adjustcolor("grey", alpha.f=0.2),ifelse(k3.mgd.sorted[rownames(x),"akrab"]==TRUE,adjustcolor("green", alpha.f=0.2),adjustcolor("red",alpha.f=0.2))), main="DIM4 to DIM5"), cex=0.5)
dev.off()


k3.pca.sc <- prcomp( k3.mgd[, 7:ncol(k3.mgd)], scale=TRUE )
pdf("analysis_results/pca_sc_k3.pdf", width=12, height=8)
plot(k3.pca.sc)
with(k3.pca.sc, plot(x[,1], x[,2], col=ifelse(is.na(k3.mgd[rownames(x),"akrab"]),adjustcolor("grey", alpha.f=0.2),ifelse(k3.mgd[rownames(x),"akrab"]==TRUE,adjustcolor("green", alpha.f=0.2),adjustcolor("red",alpha.f=0.2))), main="DIM1 to DIM2"), cex=0.5)
with(k3.pca.sc, plot(x[,2], x[,3], col=ifelse(is.na(k3.mgd[rownames(x),"akrab"]),adjustcolor("grey", alpha.f=0.2),ifelse(k3.mgd[rownames(x),"akrab"]==TRUE,adjustcolor("green", alpha.f=0.2),adjustcolor("red",alpha.f=0.2))), main="DIM2 to DIM3"), cex=0.5)
with(k3.pca.sc, plot(x[,4], x[,5], col=ifelse(is.na(k3.mgd[rownames(x),"akrab"]),adjustcolor("grey", alpha.f=0.2),ifelse(k3.mgd[rownames(x),"akrab"]==TRUE,adjustcolor("green", alpha.f=0.2),adjustcolor("red",alpha.f=0.2))), main="DIM4 to DIM5"), cex=0.5)
dev.off()

k4.mgd.sorted <- k4.mgd[k4.mgd$size<512,]
k4.pca.sc.sorted <- prcomp( k4.mgd.sorted[, 7:ncol(k4.mgd.sorted)], scale=TRUE )
pdf("analysis_results/pca_sc_k4_less_values.pdf", width=12, height=8)
plot(k4.pca.sc.sorted)
with(k4.pca.sc.sorted, plot(x[,1], x[,2], col=ifelse(is.na(k4.mgd.sorted[rownames(x),"akrab"]),adjustcolor("grey", alpha.f=0.2),ifelse(k4.mgd.sorted[rownames(x),"akrab"]==TRUE,adjustcolor("green", alpha.f=0.2),adjustcolor("red",alpha.f=0.2))), main="DIM1 to DIM2"), cex=0.5)
with(k4.pca.sc.sorted, plot(x[,2], x[,3], col=ifelse(is.na(k4.mgd.sorted[rownames(x),"akrab"]),adjustcolor("grey", alpha.f=0.2),ifelse(k4.mgd.sorted[rownames(x),"akrab"]==TRUE,adjustcolor("green", alpha.f=0.2),adjustcolor("red",alpha.f=0.2))), main="DIM2 to DIM3"), cex=0.5)
with(k4.pca.sc.sorted, plot(x[,4], x[,5], col=ifelse(is.na(k4.mgd.sorted[rownames(x),"akrab"]),adjustcolor("grey", alpha.f=0.2),ifelse(k4.mgd.sorted[rownames(x),"akrab"]==TRUE,adjustcolor("green", alpha.f=0.2),adjustcolor("red",alpha.f=0.2))), main="DIM4 to DIM5"), cex=0.5)
dev.off()


k4.pca.sc <- prcomp( k4.mgd[, 7:ncol(k4.mgd)], scale=TRUE )
pdf("analysis_results/pca_sc_k4.pdf", width=12, height=8)
plot(k4.pca.sc)
with(k4.pca.sc, plot(x[,1], x[,2], col=ifelse(is.na(k4.mgd[rownames(x),"akrab"]),adjustcolor("grey", alpha.f=0.2),ifelse(k4.mgd[rownames(x),"akrab"]==TRUE,adjustcolor("green", alpha.f=0.2),adjustcolor("red",alpha.f=0.2))), main="DIM1 to DIM2"), cex=0.5)
with(k4.pca.sc, plot(x[,2], x[,3], col=ifelse(is.na(k4.mgd[rownames(x),"akrab"]),adjustcolor("grey", alpha.f=0.2),ifelse(k4.mgd[rownames(x),"akrab"]==TRUE,adjustcolor("green", alpha.f=0.2),adjustcolor("red",alpha.f=0.2))), main="DIM2 to DIM3"), cex=0.5)
with(k4.pca.sc, plot(x[,4], x[,5], col=ifelse(is.na(k4.mgd[rownames(x),"akrab"]),adjustcolor("grey", alpha.f=0.2),ifelse(k4.mgd[rownames(x),"akrab"]==TRUE,adjustcolor("green", alpha.f=0.2),adjustcolor("red",alpha.f=0.2))), main="DIM4 to DIM5"), cex=0.5)
dev.off()

comp.sc.lm <- apply(k4.pca.sc$x, 2, function(x){
    lm(x ~ k1.mgd$sp)
})

comp.sc.lm.sum <- lapply( comp.sc.lm, function(x){ summary(x)$coefficients })

## are regions dominated by di or tetranucleotides?
hist( apply( k4.mgd[,5:ncol(k4.mgd)], 1, max ))

hist( apply( k2.mgd[,5:ncol(k2.mgd)], 1, max ))

hist( apply( k4.mgd[,5:ncol(k4.mgd)], 1, function(x){ sum(sort(x, decreasing=TRUE)[1:4])}))

## make a tree from distances of the kmers
## Don't do this, 38,000 by 38,000 distances is too much!
## k4.dst <- dist( as.matrix(k4.mgd[,5:ncol(k4.mgd)]) )



##k2.mgd.sp.pca <- prcomp(k2.mgd.sp, center=FALSE)
k2.mgd.sp.pca <- prcomp(k2.mgd.sp[1:nrow(k2.mgd.sp),1:(ncol(k2.mgd.sp)-2)], scale=TRUE)
plot(k2.mgd.sp.pca)
##with(k2.mgd.sp.pca, plot(x[,1], x[,2], col=topo.colors(length(sp))))

pdf("analysis_results/pca_sp_k2.pdf", width=12, height=8)
plot(k2.mgd.sp.pca)
##with(k2.mgd.sp.pca, plot(x[,1], x[,2], col=ifelse(is.na(k2.mgd.sp[rownames(x),"order"]),"grey",ifelse(k2.mgd.sp[rownames(x),"akrab"]==TRUE,"green","red"))))
with(k2.mgd.sp.pca, plot(x[,1], x[,2], col=ifelse(is.na(k2.mgd.sp[rownames(x),"akrab"]),"grey",ifelse(k2.mgd.sp[rownames(x),"akrab"]==TRUE,"green","red"))))
legend(x=-0.24, y=0.16, legend=c("NA", "aKrab", "no aKrab"), fill=c("grey","green","red"))
dev.off()

with(k2.mgd.sp.pca, identify(x[,1], x[,2], rownames(x)))

##k3.mgd.sp.pca <- prcomp(k3.mgd.sp, center=FALSE)
k3.mgd.sp.pca <- prcomp(k3.mgd.sp[1:nrow(k3.mgd.sp),1:(ncol(k3.mgd.sp)-2)], center=FALSE)
plot(k3.mgd.sp.pca)

pdf("analysis_results/pca_sp_k3.pdf", width=12, height=8)
plot(k3.mgd.sp.pca)
with(k3.mgd.sp.pca, plot(x[,1], x[,2], col=ifelse(is.na(k3.mgd.sp[rownames(x),"order"]),"grey",ifelse(k3.mgd.sp[rownames(x),"akrab"]==TRUE,"green","red"))))
legend(x=-0.11, y=0.16, legend=c("NA", "aKrab", "no aKrab"), fill=c("grey","green","red"))
dev.off()

##with(k3.mgd.sp.pca, plot(x[,1], x[,2], col=topo.colors(length(sp))))
##with(k3.mgd.sp.pca, identify(x[,1], x[,2], rownames(x)))


k4.mgd.sp.pca <- prcomp(k4.mgd.sp[1:nrow(k4.mgd.sp),1:(ncol(k4.mgd.sp)-2)], center=FALSE )
plot(k4.mgd.sp.pca)

pdf("analysis_results/pca_sp_k4.pdf", width=12, height=8)
plot(k4.mgd.sp.pca)
with(k4.mgd.sp.pca, plot(x[,1], x[,2], col=ifelse(is.na(k4.mgd.sp[rownames(x),"order"]),"grey",ifelse(k4.mgd.sp[rownames(x),"akrab"]==TRUE,"green","red"))))
legend(x=-0.05, y=0.16, legend=c("NA", "aKrab", "no aKrab"), fill=c("grey","green","red"))
##with(k4.mgd.sp.pca, plot(x[,1], x[,2]))
##with(k4.mgd.sp.pca, identify(x[,1], x[,2], rownames(x)))
dev.off()



plot.dims <- function(pc, col, id.leg=NULL){
    x <- matrix(1:ncol(pc$x),nrow=nrow(pc$x),ncol=ncol(pc$x), byrow=TRUE)
    plot(x, pc$x, col=col)
}

plot.dims(k2.mgd.sp.pca, col=ifelse(is.na(k2.mgd.sp$akrab), "grey", 1+k2.mgd.sp$akrab))

k2.cols <- ncol(k2.mgd.sp)-2

x2.dist <- dist(k2.mgd.sp[,1:k2.cols])

x2.dist.nj <- nj(x2.dist)

x2.dist.nj$tip.label <- paste(x2.dist.nj$tip.label, k2.mgd.sp$order)

plot(x2.dist.nj, cex=1)



k4.mgd.sp.dst <- dist( k4.mgd.sp )
k4.mgd.sp.nj <- nj( k4.mgd.sp.dst )
plot(k4.mgd.sp.nj)


sp_tree <- read.tree("species.nwk")

plot.phylo(sp_tree, edge.color=rainbow(8))
plot.phylo(sp_tree, node.color=rainbow(2))
plot.phylo(sp_tree, tip.color=topo.colors(length(sp)))

##Save lengths of all regions
region_length <- list()
for(i in 1:length(data)){
    length <- lapply(data[[i]], function(scaffold){
        scaffold$region$size
    })
    region_length[[i]] <- length
}

plot(y=region_length, x=1:length(region_length))

plot(unlist(region_length[[1]]))


##Compare k-mer frequencies
pdf("results/test.pdf")
par(mfrow=c(3,3))
for(i in 1:length(data)){
    for(j in 1:length(data[[i]])){
        col_n <- ncol(data[[i]][[j]]$spectra_k2)
        csum <- colSums(data[[i]][[j]]$spectra_k2)
        sum_freq <- csum/sum(csum)
        barplot(sum_freq, main=paste(data[[i]][[j]]$organism_name,"Scaff.",j))
    }
}
dev.off()


##Compare k-mer frequencies, with sums for all scaffolds of a given genome
pdf("analysis_results/kmer_spectra2.pdf", height=10, width=14)
par(mfrow=c(2,2))
for(i in 1:length(data)){
    cum_sum_k2 <- 0
    cum_sum_k3 <- 0
    cum_sum_k4 <- 0
    csum_nucfreq <- 0
    ##reg_sizes <- vector(length=length(data[[i]]))
    reg_sizes <- list()
##    reg_ent <- list()
    for(j in 1:length(data[[i]])){
        print(paste0("I:",i,"//J:",j))
        csum_k2 <- colSums(data[[i]][[j]]$spectra_k2)
        cum_sum_k2 <- cum_sum_k2 + csum_k2
        csum_k3 <- colSums(data[[i]][[j]]$spectra_k3)
        cum_sum_k3 <- cum_sum_k3 + csum_k3
        csum_k4 <- colSums(data[[i]][[j]]$spectra_k4)
        cum_sum_k4 <- cum_sum_k4 + csum_k4
        csum_nucfreq <- csum_nucfreq + colSums(data[[i]][[j]]$spectra_k1)
##        reg_ent[[j]] <- data[[i]][[j]]$region$ent
        reg_sizes[[j]] <- data[[i]][[j]]$region$size
    }
    organism_name <- data[[i]][[1]]$organism
    organism_name <- gsub("_"," ",organism_name)
    order_name <- sp_orders[organism_name]
##    has_akrab <- order_name%in%sp_no_akrab
    has_akrab <- k2.mgd.sp[organism_name,"akrab"]
    barplot(cum_sum_k2, main=paste0(organism_name, "\nk2 spectrum"), col.main=ifelse(is.na(has_akrab), "grey", ifelse(has_akrab, "green","red")))
    barplot(cum_sum_k3, main=paste0(organism_name, "\nk3 spectrum"), col.main=ifelse(is.na(has_akrab), "grey", ifelse(has_akrab, "green","red")))
    barplot(cum_sum_k4, main=paste0(organism_name, "\nk4 spectrum"), col.main=ifelse(is.na(has_akrab), "grey", ifelse(has_akrab, "green","red")))
##    plot(unlist(reg_ent), main=paste0(data[[i]][[1]]$organism, "\nregion entropy"))
    barplot(unlist(reg_sizes), main=paste0(organism_name, "\nreg size"), col.main=ifelse(is.na(has_akrab), "grey", ifelse(has_akrab, "green","red")))
    nuc_freq <- csum_nucfreq/sum(csum_nucfreq)
    freq_str <- paste0("Nucleotide Frequency: ",names(nuc_freq[1]),": ",round(nuc_freq[1], digits=2)," // ",names(nuc_freq[2]),": ",round(nuc_freq[2], digits=2),
                       " // ",names(nuc_freq[3]),": ",round(nuc_freq[3], digits=2)," // ",names(nuc_freq[4]),": ",round(nuc_freq[4], digits=2))
    mtext(text=freq_str, side=1, adj=0, line=1, cex=1)
    csum_nucfreq <- 0
}
dev.off()


##Generate PCA's for different kmer spectra, as well as for all species

org_names <- vector()
all_k2 <- matrix(nrow=length(data), ncol=16)
all_k3 <- matrix(nrow=length(data), ncol=64)
all_k4 <- matrix(nrow=length(data), ncol=256)
for(i in 1:length(data)){
    cum_sum_k2 <- 0
    cum_sum_k3 <- 0
    cum_sum_k4 <- 0
    org_names[[i]] <- data[[i]][[1]]$organism
    for(j in 1:length(data[[i]])){
        csum_k2 <- colSums(data[[i]][[j]]$spectra_k2)
        cum_sum_k2 <- cum_sum_k2 + csum_k2
        csum_k3 <- colSums(data[[i]][[j]]$spectra_k3)
        cum_sum_k3 <- cum_sum_k3 + csum_k3
        csum_k4 <- colSums(data[[i]][[j]]$spectra_k4)
        cum_sum_k4 <- cum_sum_k4 + csum_k4
    }
    all_k2[i,] <- cum_sum_k2
    all_k3[i,] <- cum_sum_k3
    all_k4[i,] <- cum_sum_k4
}
colnames(all_k2) <- names(cum_sum_k2)
rownames(all_k2) <- org_names
colnames(all_k3) <- names(cum_sum_k3)
rownames(all_k3) <- org_names
colnames(all_k4) <- names(cum_sum_k4)
rownames(all_k4) <- org_names

pca_k2 <- prcomp(all_k2)
pca_k3 <- prcomp(all_k3)
pca_k4 <- prcomp(all_k4)




pdf("results/PCAs.pdf", height=10, width=14)
par(mfrow=c(1,2))


pdf("results/test_ent.pdf", height=10, width=14)
par(mfrow=c(2,2))
for(i in 1:length(data)){
    reg_ent <- list()
    reg_size <- list()
    for(j in 1:length(data[[i]])){
        reg_ent[[j]] <- data[[i]][[j]]$region$ent
        reg_size[[j]] <- data[[i]][[j]]$region$size
    }
    plot(unlist(reg_ent), main=paste0(data[[i]][[1]]$organism, "\nregion entropy"))
    plot(x=unlist(reg_size), y=unlist(reg_ent), main=paste0(data[[i]][[1]]$organism, "\nReg entropy to size"))
}
dev.off()







##Test stuff

testdf <- data.frame(x=c(1,2,3,4,5,6,7,8,9,1,2,3),
                     grp=rep(c("group 1", "group 2",
                               "group 3", "group 4"),
                             each=3),
                     subgroup=LETTERS[1:3])

testdf <- reshape(testdf, idvar="subgroup",
                  timevar="grp",
                  direction="wide")

row.names(testdf) <- testdf$subgroup

tstdf <- testdf[ , 2:ncol(testdf)]

colnames(testdf) <- c("subgroup", "group1", "group2", "group3", "group4")

test <- as.matrix(testdf)

pca <- prcomp(data[[1]][[1]]$spectra_k2)
