window_s <- 512
k <- 2
synth <- paste0(rand.nseq(5000), kmer.rep("AT", 200), rand.nseq(5000), kmer.rep("TCA", 100), rand.nseq(2500), kmer.rep("GTAC", 100), rand.nseq(5000)) 
synthcmplx <- kmer.complexity(seq=synth, k=k, window.size=window_s, olap=TRUE)
png("synth_example_1")
plot.entropy(synthcmplx)
ent_stats <- paste0("window size:", window_s, "\nk:", k)
mtext(text=ent_stats, cex=1, side=3, line=0, adj=1, col="black", family="mono")
dev.off()

synth <- paste(rand.nseq(1000), kmer.rep("CG", 250), rand.nseq(1000), sep="")
synthcmplx <- kmer.complexity(seq=synth, k=2, window.size=512, olap=TRUE)
plot.entropy(synthcmplx)

synth2 <- paste(rand.nseq(10000), kmer.rep("CG", 450), rand.nseq(10000), kmer.rep("TGGT", 650), rand.nseq(2000), sep="")
synthcmplx <- kmer.complexity(seq=synth2, k=2, window.size=1024, olap=TRUE)
plot.entropy(synthcmplx)

synth3 <- paste(rand.nseq(12000), kmer.rep("CG", 100), rand.nseq(12300), kmer.rep("TGGT", 228), rand.nseq(1500), kmer.rep("AA", 100), sep="")
synthcmplx <- kmer.complexity(seq=synth3, k=2, window.size=1024, olap=TRUE)
plot.entropy(synthcmplx)

n <- c(0.35, 0.35, 0.15, 0.15)

synth4 <- paste(rand.nseq(12000, n), kmer.rep("CG", 800), rand.nseq(12300, n), kmer.rep("TGGT", 850), rand.nseq(1500), kmer.rep("AA", 100), rand.nseq(2000), sep="")
synthcmplx <- kmer.complexity(seq=synth4, k=2, window.size=1024, olap=TRUE)
plot.entropy(synthcmplx)

synth5 <- kmer.mutate(DNAString(synth4), 10000)
synthcmplx <- kmer.complexity(seq=as.character(synth5), k=2, window.size=1024, olap=TRUE)
plot.entropy(synthcmplx)

synth6 <- paste(rand.nseq(12000), rand.nseq(15000, f=c(0.5, 0.25, 0.125, 0.125)), kmer.rep("TGT", 550), rand.nseq(1500), rand.nseq(2000), sep="")
synthcmplx <- kmer.complexity(seq=synth6, k=2, window.size=1024, olap=TRUE)
plot.entropy(synthcmplx) 

synth7 <- paste(rand.nseq(12000, n), kmer.rep("CG", 800), rand.nseq(12300, n), kmer.rep("TGGT", 850), rand.nseq(1500), kmer.rep("AA", 100), rand.nseq(2000), kmer.rep("TGC", 254), rand.nseq(8500), sep="")
synth7 <- kmer.mutate(DNAString(synth7), 10000)
synthcmplx <- kmer.complexity(seq=as.character(synth7), k=2, window.size=1024, olap=TRUE)
plot.entropy(synthcmplx)

synth8 <- paste(rand.nseq(12000, n), kmer.rep("CG", 800), kmer.rep("TGC", 1000), kmer.rep("GCTA", 500), rand.nseq(12300, n), kmer.rep("TGGT", 850), rand.nseq(1500), kmer.rep("AA", 100), rand.nseq(2000), kmer.rep("TGC", 254), rand.nseq(8500), sep="")
synth8 <- kmer.mutate(DNAString(synth8), 10500)
synthcmplx <- kmer.compleity(seq=as.character(synth8), k=2, window.size=1024, olap=TRUE)
plot.entropy(synthcmplx)

synth9 <- paste0(rand.nseq(400, n), kmer.rep("CG", 100), rand.nseq(450, n))
synth9 <- kmer.mutate(DNAString(synth9), 200)
synthcmplx <- kmer.complexity(seq=as.character(synth9), k=2, window.size=256, olap=TRUE)
plot.entropy(synthcmplx)




######################
dna_data_i <- 32
window_s <- 512
cmplx <- kmer.complexity(seq=tel_dna_char[[i]], k=2, window.size=window_s, olap=TRUE)

cmplx <- kmer.complexity(seq=tel_dna_char[[i]], k=2, window.size=window_s, olap=TRUE)
regions <- vector.regions(cmplx[[1]], 2, window_s)
regions_b <- vector.borders(regions, cmplx[[1]], 2)
#######

dna_data_i <-  12
window_s <- 512
cmplx <- kmer.complexity(seq=tel_dna_char[[1]], k=4, window.size=512, olap=TRUE)

regions <- vector.regions(cmplx[[1]], 2, window_s)
regions_b <- vector.borders(regions, cmplx[[1]], 2)
########
