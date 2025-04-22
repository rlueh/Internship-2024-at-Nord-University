generate.entropy.plots <- function(dna_data, dna_data_i,  n,  k, window.size){
    organism_name <- dna_data$organism_name[1]
    dir.create(paste0(organism_name, "_entropy"))
    for(i in 1:n){
        png(paste0(organism_name,"_entropy/",organism_name,"_", i, ".png"), type="Xlib", width=2000, height=500)
        ent <- kmer.complexity(dna_data$dna[i], k=k, window.size=window_s, olap=TRUE)
        plot.entropy(ent, main=organism_name)
        ent_stats <- paste0("window size:", window_s, "\nk:", k, "\nindex:", dna_data_i)
        mtext(text=ent_stats, cex=1, side=3, line=0, adj=1, col="black", family="mono")
        dev.off()
    }
}

generate.fft.plots <- function(seq_fft, data, organism_name){
    pdf(paste0(organism_name,"_fft.pdf"), width=14, height=6)
    par(mar=c(5, 4, 10, 2)+0.1)
    for (i in 1:length(seq_fft)){
        fft.plot(seq_fft[[i]], data$seq[i], data$reg[i,], data$score[i,], data$ent[i], q.col=0.9)
    }
    dev.off()
}

generate.region.plots <- function(file_name, ent, data, seq_chromosome){
##    sma <- region.smoothing(ent)
    pdf(file_name, width=14, height=6)
    par(mar=c(5, 4, 10, 2)+0.1)
    for (i in 1:length(data$seq)){
        r <- (data$reg[i,]$lborder-1000):(data$reg[i,]$stop+1000)
         plot(r, ent[[1]][r], type="l")
        ## plot(r, sma[r], type="l")
       
        plot.drawseq(seq=seq_chromosome, reg=data$reg[i,], ent=ent[[1]], offset=800)
        ent_text <- paste0("Ent:",data$ent[i])
        mtext(text=ent_text, cex=0.7, side=3, line=0, adj=0, col="black", family="mono")
    }
    dev.off()
}

generate.fft.shapes <- function(seq, descr){
    synth_fft <- fft.get(as.character(seq))
    fft_sum <- rowSums(Mod(synth_fft[[1]]))
    plot(fft_sum, type="l", col=1)
    mtext(text=descr, cex=1, side=3, line=0, adj=0, col="black", family="mono")
}

plot.entropy.length <- function(reg, ent){
    
}

generate.fft.mutate.plots <- function(seq_synth, name,  plot_cols=TRUE, center=TRUE, use_peaks=FALSE){
    prob_s <- 0
    prob_i <- 0
    prob_d <- 0
    pdf(name, width=14, height=6)
    for(i in 1:11){
        seq_synth_m  <- seq.mutate(DNAString(seq_synth), prob_s=prob_s, prob_i=prob_i, prob_d=prob_d)    
        ##synth_fft <- fft.get.single(as.character(seq_synth_m), center=center)
        synth_fft <- fft.get(as.character(seq_synth_m), center=center)
        synth_fft <- synth_fft[[1]]
        if(use_peaks){
            peaks <- fft.peaks(synth_fft)
##            print("Peaks retrieved")
##            print(peaks)

            extract_seq_maxval <- paste(rev.fft(synth_fft, qnt=0.98), collapse="")
            if(!is.null(peaks)){
                peaks <- append(1, peaks)
                extract_seq <- paste(rev.fft.peaks(synth_fft, peaks), collapse="")
            } else{
                extract_seq <- paste("Peaks could not be found")
            }
        } else{
            ##extract_seq <- paste(rev.fft(synth_fft[[1]], qnt=0.98), collapse="")
            extract_seq <- paste(rev.fft(synth_fft, qnt=0.98, fft_single=TRUE), collapse="")
        }
        ##fft_sum <- rowSums(Mod(synth_fft[[1]]))
        fft_sum <- rowSums(Mod(synth_fft))
        if(plot_cols){
            if(use_peaks){
                fft.plot(synth_fft, peaks=peaks)
            } else{
                ##fft.plot(synth_fft[[1]], q.col=0.98)
                fft.plot(synth_fft, q.col=0.98)
            }
            lines(Mod(synth_fft[,1]), col="red")
            lines(Mod(synth_fft[,2]), col="blue")
            lines(Mod(synth_fft[,3]), col="purple")
            lines(Mod(synth_fft[,4]), col="yellow")
            legend(x="topright", legend=c("rowSums", "A", "C", "G", "T"), fill=c("black", "red", "blue", "purple", "yellow"), cex=0.8)
        } else{
            ##fft.plot(synth_fft[[1]], q.col=0.98)
            if(use_peaks){
                fft.plot(synth_fft, peaks=peaks)
                fft_smooth <- convolve(rowSums(Mod(synth_fft)), dnorm(-20:20, mean=0, sd=4), type='filter')
                lines(1:(length(rowSums(Mod(synth_fft)))-40)+21, fft_smooth, col="blue")
                abline(h=mean(fft_smooth), col="green")
            } else {
                fft.plot(synth_fft, q.col=0.98)
            }
        }
        
        print_seq <- sequence.insertlines(as.character(seq_synth_m), l_max=75)
        print_seq_ex <- sequence.insertlines(extract_seq, l_max=75)
        mtext(text=print_seq, cex=0.4, side=3, line=0, adj=0.2, col="red", family="mono")
        mtext(text=print_seq_ex, cex=0.4, side=3, line=0, adj=1, col="blue", family="mono")
        stats_text <- paste0("Ent:", entropy.calculate(seq_synth_m, k=4), "\n", "prob_s:", prob_s, "\n", "prob_i:", prob_i, "\n", "prob_d:", prob_d, "\n", "adist_peaks:",
                             adist(seq_synth,extract_seq), "\n", "adist_maxval: ", adist(seq_synth, extract_seq_maxval))
        mtext(text=stats_text, cex=0.7, side=3, line=0, adj=0, col="black", family="mono")
        prob_s <- prob_s + 0.05
        prob_i <- prob_i + 0.005
        prob_d <- prob_d + 0.005
    }
    dev.off()
}

## generate.fft.mutate.plots <- function(seq_synth, name,  plot_cols=TRUE, center=TRUE, use_peaks=FALSE){
##     prob_s <- 0
##     prob_i <- 0
##     prob_d <- 0
##     pdf(name, width=14, height=6)
##     for(i in 1:21){
##         seq_synth_m  <- seq.mutate(DNAString(seq_synth), prob_s=prob_s, prob_i=prob_i, prob_d=prob_d)    
##         synth_fft <- fft.get.single(as.character(seq_synth_m), center=center)
##         if(use_peaks){
##             peaks <- fft.peaks(synth_fft[[1]], n=3)
##             extract_seq <- paste(rev.fft.peaks(synth_fft[[1]], peaks), collapse="")
##         } else{
##             ##extract_seq <- paste(rev.fft(synth_fft[[1]], qnt=0.98), collapse="")
##             extract_seq <- paste(rev.fft(synth_fft, qnt=0.98), collapse="")
##         }
##         ##fft_sum <- rowSums(Mod(synth_fft[[1]]))
##         fft_sum <- rowSums(Mod(synth_fft))
##         if(plot_cols){
##             ##par(mfrow=c(1,2))
##             mar <- c(0.5,0.5,1,0.5)
##             par(fig=c(0, 0.2, 0, 0.5), mar=mar)
##             ##plot(colSums(Mod(synth_fft[[1]])), type="p", col=1)
##            ##plot(Mod(synth_fft[[1]][,1]), ylab="", yaxt="n", xlab="", xaxt="n", main="A", type="l")
##             plot(Mod(synth_fft[,1]), ylab="", yaxt="n", xlab="", xaxt="n", main="A", type="l")
##             par(fig=c(0.2, 0.4, 0, 0.5), mar=mar, new=TRUE)
##             ##plot(Mod(synth_fft[[1]][,2]), ylab="", yaxt="n", xlab="", xaxt="n", main="C", type="l")
##             plot(Mod(synth_fft[,2]), ylab="", yaxt="n", xlab="", xaxt="n", main="C", type="l")
##             par(fig=c(0, 0.2, 0.5, 1), mar=mar, new=TRUE)
##             ##plot(Mod(synth_fft[[1]][,3]), ylab="", yaxt="n", xlab="", xaxt="n", main="G", type="l")
##             plot(Mod(synth_fft[,3]), ylab="", yaxt="n", xlab="", xaxt="n", main="G", type="l")
##             par(fig=c(0.2, 0.4, 0.5, 1), mar=mar, new=TRUE)
##             ##plot(Mod(synth_fft[[1]][,4]), ylab="", yaxt="n", xlab="", xaxt="n", main="T", type="l")
##             plot(Mod(synth_fft[,4]), ylab="", yaxt="n", xlab="", xaxt="n", main="T", type="l")
             
##             par(fig=c(0.4, 1, 0, 1), mar=c(4, 4, 4, 1), new=TRUE)

##             if(use_peaks){
##                 fft.plot(synth_fft[[1]], peaks=peaks)
##             } else{
##                 ##fft.plot(synth_fft[[1]], q.col=0.98)
##                 fft.plot(synth_fft, q.col=0.98)
##             }
##          } else{
##              ##fft.plot(synth_fft[[1]], q.col=0.98)
##              fft.plot(synth_fft, q.col=0.98)
##          }

        
##         print_seq <- sequence.insertlines(as.character(seq_synth_m), l_max=75)
##         print_seq_ex <- sequence.insertlines(extract_seq, l_max=75)
##         mtext(text=print_seq, cex=0.4, side=3, line=0, adj=0.2, col="red", family="mono")
##         mtext(text=print_seq_ex, cex=0.4, side=3, line=0, adj=1, col="blue", family="mono")
##         stats_text <- paste0("Ent:", entropy.calculate(seq_synth_m, k=4), "\n", "prob_s:", prob_s, "\n", "prob_i:", prob_i, "\n", "prob_d:", prob_d)
##         mtext(text=stats_text, cex=0.7, side=3, line=0, adj=0, col="black", family="mono")
##         prob_s <- prob_s + 0.05
##         prob_i <- prob_i + 0.005
##         prob_d <- prob_d + 0.005
##     }
##     dev.off()
## }

