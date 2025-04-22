source("nuc_complexity/nuc_complexity.R")
source("functions.R")
source("create_plots.R")
require("Biostrings")


##Index to work with
dna_data_i <- 26
dna_data <- read.fasta(dna_data_i)
organism_name <- dna_data$organism_name[1]
##index 4, w=512, k=2 as example for weird plot

ent <- kmer.complexity(seq=dna_data$dna[1], k=2, window.size=512, olap=TRUE)

window_s <- 512
k <- 2
generate.entropy.plots(dna_data, dna_data_i, n=2, k=k, window.size=window_s)

data <- entropy.analyze(dna_data, 512, k=2)


plot.sizefreq(data, k=2, organism_name)


k=2

plot(data[[1]]$regions$size, data[[1]]$scores[,k-1])
identify(x=data[[1]]$regions$size, y=data[[1]]$scores[,k-1])

####
synth_seq <- kmer.rep("TGG", 50)
seq_enc <- fft.encode(synth_seq, center=FALSE)


####

print(paste("Organism:", organism_name))

seq <- paste(rep("AGCTCAACACACACACACACACGTGT", 100), collapse="")
seq <- paste(rep("CGAATTT", 100), collapse="")


seq_fft <- fft.get(data[[1]]$seq)


generate.fft.plots(seq_fft, data[[1]], organism_name)



generate.region.plots(paste0(organism_name, "_data_k2_w256", ".pdf"), ent, data[[1]], dna_data$dna[1])


seq_synth <- kmer.rep("AT", 300)
generate.fft.mutate.plots(seq_synth, "fft_sum_col_plot_AT_peaks.pdf", center=FALSE, use_peaks=TRUE, plot_cols=FALSE)

seq_synth <- kmer.rep("TGG", 200)
generate.fft.mutate.plots(seq_synth, "fft_sum_col_plot_TGG_peaks_new.pdf", center=FALSE, use_peaks=TRUE, plot_cols=FALSE)

seq_synth <- kmer.rep("A", 500)
generate.fft.mutate.plots(seq_synth, "fft_sum_col_plot_A_new.pdf", center=FALSE)



mut_arr <- vector(length=10)
prob_s <- 0
prob_i <- 0
prob_d <- 0
for(i in 1:10){
    mut_arr[i] <- as.character(seq.mutate(DNAString(seq_synth), prob_s=prob_s, prob_i=prob_i, prob_d=prob_d))
    prob_s <- prob_s + 0.1
    prob_i <- prob_i + 0.01
    prob_d <- prob_d + 0
}

pdf("recovery_comparison_singlecol_2.pdf", width=14, height=6)
for(i in 1:length(mut_arr)){
    fft <- fft.get.single(as.character(mut_arr[i]))
    ##fft <- fft.get(as.character(mut_arr[i]))
    ##fft <- fft[[1]]
    extract_seq <- paste(rev.fft(fft, qnt=0.98, fft_single=TRUE), collapse="")
    fft.plot(fft, q.col=0.98)
    lines(Mod(fft[,1]), col="red")
    lines(Mod(fft[,2]), col="blue")
    lines(Mod(fft[,3]), col="purple")
    lines(Mod(fft[,4]), col="yellow")
    legend(x="topright", legend=c("rowSums", "A", "C", "G", "T"), fill=c("black", "red", "blue", "purple", "yellow"), cex=0.8)
    print_seq <- sequence.insertlines(as.character(mut_arr[i]), l_max=100)
    print_seq_ex <- sequence.insertlines(extract_seq, l_max=100)
    mtext(text=print_seq, cex=0.4, side=3, line=0, adj=0.2, col="red", family="mono")
    mtext(text=print_seq_ex, cex=0.4, side=3, line=0, adj=1, col="blue", family="mono")
    stats_text <- paste0("Ent:", entropy.calculate(mut_arr[i], k=4), "\n", "Diff:", adist(seq_synth, extract_seq))
    mtext(text=stats_text, cex=0.7, side=3, line=0, adj=0, col="black", family="mono")
}
dev.off()




synth_fft <- fft.get(as.character(seq_synth))
par(mfrow=c(2,2))
plot(Mod(synth_fft[[1]][,1]))
plot(Mod(synth_fft[[1]][,2]))
plot(Mod(synth_fft[[1]][,3]))
plot(Mod(synth_fft[[1]][,4]))

########################

synth_AT <- kmer.rep("AT", 300)

fft_norm <- fft.get(synth_AT)
plot(rowSums(Mod(fft_norm[[1]])))

mut_AT <- as.character(seq.mutate(DNAString(synth_AT), prob_s=0.5, prob_i=0.015, prob_d=0.001))

mut_AT <- as.character(seq.mutate(DNAString(synth_AT), prob_s=0.9, prob_i=0.045, prob_d=0.05))

fft_AT <- fft.get(mut_AT)
fft_AT <- fft_AT[[1]]

fft.plot(fft_AT)

fft_sma <- stats::filter(Mod(fft_AT), filter=c(1/10, 1/5, 1/5, 1/10 ), sides=1, circular=TRUE)
plot(rowSums(fft_sma), type="l")

smooth_AT <- smoother(rowSums(Mod(fft_AT)), n=20,  plot=FALSE, method="supersmooth")
##smooth_AT <- smoother(smooth_AT, n=20, plot=FALSE)
plot(smooth_AT, type="l")

mod_fft <- rowSums(Mod(fft_AT))

span <- 5/length(mod_fft)
smooth_AT <- supsmu(x=1:length(mod_fft), y=mod_fft, span=span)
smooth_AT <- smooth_AT$y
plot(smooth_AT, type="l")


smooth_AT <- convolve(mod_fft, dnorm(-20:20, mean=0, sd=5), type='filter')

plot(smooth_AT, type="l")

plot(mod_fft, type="l")
lines(1:(length(mod_fft)-40)+21, smooth_AT, col="red")

fft_diff <- diff(smooth_AT, lag=1, differences=1)

lines(2:(length(mod_fft)-40)+21, fft_diff, col="blue")



peaks_AT <- fft.peaks(fft_AT)

peaks_AT <- append(1,peaks_AT)

rev_AT <- paste(rev.fft.peaks(fft_AT, peaks_AT), collapse="")
rev_AT_old <- paste(rev.fft(fft_AT, qnt=0.98), collapse="")

###
synth_TGG <- kmer.rep("TGG", 300)
fft_TGG_n <- fft.get(as.character(synth_TGG))
fft_TGG_sums_n <- rowSums(Mod(fft_TGG_n[[1]]))
plot(fft_TGG_sums_n, type="l")
##fft_sums_[301] = 1448.528


##Mutate sequence
mutate_TGG <- as.character(seq.mutate(DNAString(synth_TGG), prob_s=0.4, prob_i=0.015, prob_d=0.001))

fft_TGG <- fft.get(mutate_TGG)

##Manually select peak area
##fft_sums[304] = 509.6775
##fft_sums[612] = 509.6775
##sum(fft_sums[304:308]) = 1413.734 // fft_sums[608:612]
##606:613 / 303:310

mod_TGG <- rowSums(Mod(fft_TGG[[1]]))
plot(mod_TGG, type="l")
                   
fft.plot(fft_TGG[[1]])

##!!
peaks <- c(1, 303:314, 605:616)

peaks <- fft.peaks(fft_TGG[[1]], max_n=4)

peaks <- fft.peaks(fft_TGG[[1]])

peaks <- append(1,peaks)

rev_peaks <- paste(rev.fft.peaks(fft_TGG[[1]], peaks), collapse="")

adist(synth_TGG, rev_peaks)


rev_seq <- paste(rev.fft(fft_TGG[[1]], qnt=0.98), collapse="")

seq_fft <- fft.get(as.character(seq_synth_m), center=TRUE)
fft_sums <- rowSums(Mod(seq_fft[[1]]))
plot(fft_sums, type="l")

fft_sma <- stats::filter(Mod(seq_fft[[1]]), filter=c(1/10, 1/5, 1/5, 1/10 ), sides=1, circular=TRUE)
plot(rowSums(fft_sma), type="l")



fft.plot(seq_fft[[1]], q.col=0.98)

fft.plot(seq_fft[[1]], peaks=peaks)

fft.plot(tmp, peaks=peaks)

fft.plot(tmp, q.col=0.98)

fft.plot(seq_fft[[1]], q.col=0.98)

fft_mod <- Mod(seq_fft[[1]])



##fft_sma <- stats::filter(fft_mod, filter = rep(1, 5) / 5, circular=FALSE)
sma_TGG <- stats::filter(Mod(fft_TGG[[1]]), filter=c(1/10, 1/5, 1/5, 1/10 ), sides=1, circular=TRUE)
##fft_sma <- lowess(rowSums(fft_mod), f=1/20)
plot(rowSums(sma_TGG), type="l")


source("smoother.R")

smooth_TGG <- smoother(rowSums(Mod(fft_TGG[[1]])), n=20,  plot=FALSE, method="supersmooth")
plot(smooth_TGG, type="l")


smooth_TGG <- convolve(mod_TGG, dnorm(-20:20, mean=0, sd=5), type='filter')

plot(smooth_TGG, type="l")


plot(mod_TGG, type="l")
lines(1:(length(mod_TGG)-40)+21, smooth_TGG, col="red")



fft.plot(smooth_TGG)
##fft_sma <- decompose(fft_mod)
##fft_sums <- rowSums(fft_sma)



pos <- fft.peaks(fft_sma$y)
plot(rowSums(fft_sma), type="l")



#######################

source("functions.R")
sma <- region.smoothing(ent)


par(mar=c(5, 4, 10, 2)+0.1)
r <- (data[[1]]$reg[26,]$lborder-10000):(data[[1]]$reg[26,]$stop+10000)
plot(r, sma[r], type="l")
plot.drawseq(seq=dna_data$dna[1], reg=data[[1]]$reg[26,], ent=sma, offset=800)


###


pdf("FFT_clone_shape.pdf", width=14, height=6)
generate.fft.shapes(paste0(kmer.rep("AT", 50), "A", kmer.rep("AT", 50)), "AT(50)-A-AT(50)")
generate.fft.shapes(paste0(kmer.rep("AT", 50), "T"), "AT(50)-T")
generate.fft.shapes(paste0(kmer.rep("AT", 50), "T", kmer.rep("AT", 50), "C", kmer.rep("AT", 10)), "AT(50)-T-AT(50)-C-AT(10)")
generate.fft.shapes(paste0(kmer.rep("AT", 50), "T", kmer.rep("AT", 50), "T", kmer.rep("AT", 20)), "AT(50)-T-AT(50)-T-AT(20)")
dev.off()


pdf("entropy_to_size_k6_w4096.pdf")
plot(data[[1]]$reg$size, data[[1]]$ent, xlab="Size", ylab="Entropy", main="k6, w4096")
dev.off()
