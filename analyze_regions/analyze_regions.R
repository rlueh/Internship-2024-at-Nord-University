source("nuc_complexity/nuc_complexity.R")
source("functions.R")
source("create_plots.R")
require("Biostrings")


##k is kmer sized used to calcualte entropy for individual regions
w_size=512
k=2

assemblies_info <- read.csv("../assembly_info/assemblies_info.csv")

##Only use assemblies of the Chromosome level
assemblies_info <- assemblies_info[which(assemblies_info$Assembly.Level=="Chromosome"), ]

assemblies_info <- assemblies_info[which(assemblies_info$Assembly.Stats.Contig.N50>1e6),]

dna_assemblies <- read.assemblies(assemblies_info[1:50,])

##identifier <- paste0(Sys.Date(), "-", format(Sys.time(), "%H-%M"))
for(i in 28:length(dna_assemblies$dna)){
##for(i in 2:2){
    ##Only use contigs that are bigger than 1e6
    dna_assemblies$dna[[i]] <- dna_assemblies$dna[[i]][which(width(dna_assemblies$dna[[i]])>1e6)]
    ##Remove a few indexes for testing purposes
##    if(length(dna_assemblies$dna[[i]]) > 3){
##        dna_assemblies$dna[[i]] <- dna_assemblies$dna[[i]][10]        
##    }
    print(paste("DNA Assembly:",i))
##    organism_dir <- paste0("analytics/",dna_assemblies$organism[i], "_", identifier)
    organism_dir <- paste0("../play_02/analytics2/",dna_assemblies$organism[i])
    organism_dir <- sub(" ", "_", organism_dir)
    dir.create(organism_dir)
    n_matches <- vmatchPattern("N", dna_assemblies$dna[[i]])
    matches_i <- startIndex(n_matches)
    reg_ignore <- match.findrange(matches_i,512)
    data  <- entropy.analyze(dna_assemblies$dna[[i]], w_size, k, reg_ignore)
    info_df <- data.frame("organism"=dna_assemblies$organism[i],"file_path"=names(dna_assemblies$dna[i]),"index"=i)
    write.csv(info_df, paste0(organism_dir,"/info.csv"))
    ##Generate CSV files with data
    print("##Writing to files##")
    for(j in 1:length(data)){
        print(paste("Data loop:",j))
        scaffold_dir <- paste0(organism_dir, "/", "scaffold_", j)
        dir.create(scaffold_dir)
        write.csv(data[[j]]$spec_k1, paste0(scaffold_dir, "/", "kmer_spectra_k1.csv"))
        write.csv(data[[j]]$spec_k2, paste0(scaffold_dir, "/", "kmer_spectra_k2.csv"))
        write.csv(data[[j]]$spec_k3, paste0(scaffold_dir, "/", "kmer_spectra_k3.csv"))
        write.csv(data[[j]]$spec_k4, paste0(scaffold_dir, "/", "kmer_spectra_k4.csv"))
        write.csv(data[[j]]$regions, paste0(scaffold_dir, "/", "regions.csv"))
        ##Draw plots for each chromosome
        print("Drawing plot")
        png(paste0(scaffold_dir,"/entropy_plot.png"), type="Xlib", width=1200, height=600)
        plot.entropy(data[[j]]$entdata)
        dev.off()
    }
}    



##Testing parallelization
library("parallel")
library("foreach")
library("doParallel")

dyn.load("~/introns/play_02/nuc_complexity/src/nuc_complexity.so", sep="/")

n_matches <- vmatchPattern("N", dna_assemblies$dna[[1]])
matches_i <- startIndex(n_matches)
reg_ignore <- match.findrange(matches_i,512)

##dna_assemblies$dna[[1]] <- dna_assemblies$dna[[i]][1:2]        

testlist <- entropy.analyze.p(dna_assemblies$dna[[1]], w_size, k, reg_ignore)

k=2
window_s=512

dna_data <- dna_assemblies$dna[[1]]

cluster <- makeCluster(2)
registerDoParallel(cluster)
datalist <- foreach(i = 1:length(dna_data)) %dopar% {
    print(paste("Assembly Contig:",i))
    cmplx <- kmer.complexity(seq=as.character(dna_data[[i]]), k=k, window.size=window_s, olap=TRUE)
    cmplx_sma <- region.smoothing(cmplx)
    thrsh  <- entropy.threshold(cmplx_sma, sens=TRUE, reg_ignore=reg_ignore[[i]])
    regions <- vector.regions(cmplx_sma, thrsh, window_s, reg_ignore=reg_ignore[[i]])
    regions_b <- vector.borders(regions, cmplx_sma, thrsh)
    ## Sometimes regions contain regions that go out of bounds of the string.
    ## If first region starts at 0, set starting point to window size
    if(length(which(regions_b$lborder==0))>0){
        ##regions_b <- regions_b[-(which(regions_b$lborder==0)), ]
        regions_b[1,]$lborder <- window_s
    }
    strngs <- sequence.getchars(dna_data[[i]], regions_b)
    ent <- entropy.calculate(strngs, k=k)
    regions_b <- cbind(regions_b,"ent"=ent)
    regions_b <- cbind(regions_b,"seq"=strngs)
    spec_k1 <- oligonucleotideFrequency(DNAStringSet(strngs), width=1)
    spec_k2 <- oligonucleotideFrequency(DNAStringSet(strngs), width=2)
    spec_k3 <- oligonucleotideFrequency(DNAStringSet(strngs), width=3)
    spec_k4 <- oligonucleotideFrequency(DNAStringSet(strngs), width=4)
    print("ALL STATS GENERATED")
    datalist[[i]]$regions <- regions_b
    datalist[[i]]$seq <- strngs
    datalist[[i]]$spec_k1 <- spec_k1
    datalist[[i]]$spec_k2 <- spec_k2
    datalist[[i]]$spec_k3 <- spec_k3
    datalist[[i]]$spec_k4 <- spec_k4
    datalist[[i]]$threshold <- thrsh
    datalist[[i]]$entdata <- cmplx
}
stopCluster(cl=cluster)






####################
###################
dnastr <- dna_assemblies$dna[[2]][2]
cmplx <- kmer.complexity(seq=as.character(dnastr), k=2, window.size=512, olap=TRUE)

## n-regions have complexity close to 0. This causes problems for my algorithm.
cmplx_sma <- region.smoothing(cmplx)

png("TESTPLOT.png", type="Xlib", width=1200, height=600)
plot.entropy(cmplx)
dev.off()

n_matches <- vmatchPattern("N", dna_assemblies$dna[[2]])
matches_i <- startIndex(n_matches)

reg_ignore <- match.findrange(matches_i,512)

thrsh <- entropy.threshold(cmplx_sma, sens=TRUE, reg_ignore=reg_ignore[[2]])

regions <- vector.regions(cmplx_sma, thrsh, 512, reg_ignore[[2]])

regions_b <- vector.borders(regions, cmplx_sma, thrsh)

##Throw away last one, this one is usually problematic (Exceeding the limits of the DNA Str for example)
##regions_b <- regions_b[-length(regions_b$lborder),]


if(length(which(regions_b$lborder==0))>0){
    ##regions_b <- regions_b[-(which(regions_b$lborder==0)), ]
    regions_b[1,]$lborder <- 1
}





##Careful when regions_b first index has lborder of 0. the case in [2][2] assembly 
regions_b <- regions_b[-1,]

strngs <- sequence.getchars(dnastr[[1]], regions_b)

l_freq <- letterFrequency(strngs, letters="N")
strngs <- strngs[-(which(l_freq>0))]
