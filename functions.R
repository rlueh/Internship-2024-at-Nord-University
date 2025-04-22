

read.assemblies <- function(assemblies_info, contig_size=1e6){
     tel_dna <- sapply(assemblies_info$file, function(x){
         dna_set <- readDNAStringSet(x)
         ##Remove contigs that are smaller than contig_size
         ##dna_set <- dna_set[which(width(dna_set)>=contig_size)]
         ##dna_set <- gsub("N","",dna_set)
         dna_set
     })
     ##Clean up ambiguity character "N"
     #tel_dna <- gsub("N","",tel_dna)
     list("dna"=tel_dna, "organism"=assemblies_info$Organism.Name)vers
}

match.findrange <- function(matches_i, w_size){
    regions_N <- vector(mode="list", length=length(matches_i))
    for(i in 1:length(matches_i)){
        ##print(paste("HEAD, I: ", i))
        match_diff <- diff(matches_i[[i]], lag=1)
        start <- matches_i[[i]][1]
        v_start <- vector()
        v_end <- vector()
        count <- 1
        if(!is.null(match_diff)){
            for(j in 1:length(match_diff)){
               if(match_diff[j]>1){
                   ##print(paste("MATCH_DIFF J: ",j))
                   ##ranges[count] <- c(start, end)
                   v_start[count] <- start
                   v_end[count] <- matches_i[[i]][j]
                   ##print(paste("START:",start,"//END:",end))
                   start <- matches_i[[i]][j+1]
                   count <- count+1
                 } else if(j==length(match_diff)){
                     ##print("##IN##")
                     v_start[count] <- start
                     v_end[count] <- matches_i[[i]][j]
                 }
            }
        }
        regions_N[[i]] <- data.frame("start"=v_start,"end"=v_end)
    }
    regions_N
}

entropy.analyze.p <- function(dna_data, window_s, k, reg_ignore){
    cluster <- makeCluster(2)
    registerDoParallel(cluster)
    datalist <- foreach(i = 1:length(dna_data), .export=ls(globalenv()), .packages="Biostrings") %dopar% {
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
        
         spec_k1 <- oligonucleotideFrequency(DNAStringSet(strngs), width=1, as.prob=TRUE)
         spec_k2 <- oligonucleotideFrequency(DNAStringSet(strngs), width=2, as.prob=TRUE)
         spec_k3 <- oligonucleotideFrequency(DNAStringSet(strngs), width=3, as.prob=TRUE)
         spec_k4 <- oligonucleotideFrequency(DNAStringSet(strngs), width=4, as.prob=TRUE)
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
     print("##Returning datalist##")
    stopCluster(cl=cluster)
    datalist
}
        


entropy.analyze <- function(dna_data, window_s, k, reg_ignore){
    datalist <- vector(mode="list", length=length(dna_data))
    for(i in 1:length(dna_data)){
        print(paste("Assembly Contig:",i))
        cmplx <- kmer.complexity(seq=as.character(dna_data[[i]]), k=k, window.size=window_s, olap=TRUE)
        ##thrsh  <- entropy.threshold(cmplx[[1]], sens=TRUE)
        ##regions <- vector.regions(cmplx[[1]], thrsh, window_s)
        ##regions_b <- vector.borders(regions, cmplx[[1]], thrsh)

        ##cmplx_sma <- stats::filter(cmplx[[1]], rep(1,12*10)/(12*10), circular=TRUE)
        ##Change smoothing function to take parameters for kernel
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
        ## Handle scenarios where no regions are found
        if(nrow(regions_b)!=0){
            strngs <- sequence.getchars(dna_data[[i]], regions_b)
            isEmpty <- FALSE
        } else {
            isEmpty <- TRUE
            strngs <- DNAStringSet("")
            regions_b[1,] = NA
        }
        ent <- entropy.calculate(strngs, k=k)
        regions_b <- cbind(regions_b,"ent"=ent)
        regions_b <- cbind(regions_b,"seq"=strngs)

        if(isEmpty){
            regions_b[0,]
        }
        
        spec_k1 <- oligonucleotideFrequency(DNAStringSet(strngs), width=1)
        spec_k2 <- oligonucleotideFrequency(DNAStringSet(strngs), width=2)
        spec_k3 <- oligonucleotideFrequency(DNAStringSet(strngs), width=3)
        spec_k4 <- oligonucleotideFrequency(DNAStringSet(strngs), width=4)
        print("ALL STATS GENERATED")
        datalist[[i]]$regions <- regions_b
        datalist[[i]]$seq <- strngs
        ##datalist[[i]]$scores <- t(scores)
        datalist[[i]]$spec_k1 <- spec_k1
        datalist[[i]]$spec_k2 <- spec_k2
        datalist[[i]]$spec_k3 <- spec_k3
        datalist[[i]]$spec_k4 <- spec_k4
        datalist[[i]]$threshold <- thrsh
        datalist[[i]]$entdata <- cmplx
    }
    print("##Returning datalist##")
    datalist
}

## Using "taxize" package would be better, but installing currently not possible
taxonomy.get <- function(species){
    api_key <- "a998370518e456c17a14720b557ec9922107"
    species <- gsub("_", " ", species)
    order_names <- sapply(species, function(sp){
        args <- URLencode(sp)
        req <- paste0("https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/",args,"?returned_content=METADATA")
        result <- GET(req,
                      add_headers("api-key"=api_key))
        print("result": result)
        result_json <- fromJSON(content(result,as="text"))
        print("result JSON:",  result_json)
        lineage <- result_json$taxonomy_nodes$taxonomy$lineage
        
        ##Translate lineage IDs to names and extract order
        lineage <- lineage[[1]][-1]
        args <- URLencode(paste(lineage, collapse=","))
        req <- paste0("https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/",args,"/name_report")
        result <- GET(req,
                      add_headers("api-key"=api_key))
        
        result_json <- fromJSON(content(result,as="text"))
        print(sp)
##        if(sp=="Ilyophis sp. JC 2022 1"){
##            browser()
##        }
        ##Find order in json
        if(!is.null(result_json$reports$taxonomy$rank)){
            result <- which(sapply(result_json$reports$taxonomy$rank, function(json) is.element("ORDER", json)))
            result <- result_json$reports$taxonomy$current_scientific_name$name[result]
        } else{
            result <- NA
        }
        
        result
     })
}


vector.regions  <- function(v, thrsh, beg, reg_ignore){
    r <- (1+beg):length(v)
    l <- r-1
    
    start <- which(v[l] > thrsh & v[r] <= thrsh)
    stop <- which(v[l] <= thrsh & v[r] > thrsh)
    ## Regions will be wrong if beg or end are lower than trhsh

    isFixed_start <- FALSE
    ## Fix start condition
    if(start[1] > stop[1]){
        start <- c(beg, start)
        isFixed_start <- TRUE
    }

    ## Fix stop condition
    if(stop[length(stop)] < start[length(start)]){
        stop <- c(stop, length(v))
    }
    
    ## Start & Stop should be same length, if not the logic is wrong
    if(length(stop) != length(start)){
        stop("Bugger! The lengths are not the same")
    }
    
    ## Add beg to start & stop to get full size of low-complexity region
    ## In some cases, when stop is initially below window size, start pos will get fixed to window size
    ## resulting in negative region size, check that and ignore first start position if true
    if(isFixed_start){
        start[2:length(start)] <- start[2:length(start)]+beg
    } else {
        start <- start+beg
    }
    stop <- stop+beg

    start_n <- vector()
    stop_n <- vector()
    ## modify the start and ends so that they don't overlap with the regions
    ## for(i in 1:nrow(reg_ignore)){
    ##     ## define overlaps;
    ##     ol1 <- reg_ignore$end[i] - start
    ##     ol2 <- stop - reg_ignore$start[i]
    ##     ol.b <- pmin( ol1, ol2 ) > 0
    ##     start_n[i] <- start[which(ol.b)]
    ##     stop_n[i] <- stop[which(ol.b)]
    ## }
    
    ## Check if caught regions are in an area with N's, if so ignore the region
    if(nrow(reg_ignore)!=0){
        count <- 1
        for(i in 1:nrow(reg_ignore)){
            for(j in 1:length(start)){
                is_olap <- start[j]:stop[j]%in%reg_ignore$start[i]:reg_ignore$end[i]
                if(length(which(is_olap==TRUE)>0)){
                    start_n[count] <- j
                    stop_n[count] <- j
                    count <- count+1
                }
            }
        }
        if(length(start_n) || length(stop_n) > 0){
            start <- start[-start_n]
            stop <- stop[-stop_n]
        }
    }
    
    data.frame("start"=start, "stop"=stop, "size"=1 + stop - start)
}

vector.borders <- function(v, ent, thrsh){

   # ent_med <- median(ent)

    ent_mean <- mean(ent)/1.02

    isEmpty <- FALSE
    if(nrow(v)==0){
        isEmpty <- TRUE
        v[1,] = NA
    }
    
    v$lborder = 0
    v$rborder = 0

    if(!isEmpty){
        
        ## I would instead do:
        ## borders <- matrix(nrow=nrow(v), ncol=2)
        ## colnames(border) <- c('r', 'l')
        lborder <- 0
        rborder <- 0
        
        for (reg in 1:nrow(v)){
            i <- v[reg,]$start
            ## typical to write this:
            ## v$start[reg]
            ## v[reg,'start']
            for (i in v[reg,]$start:1){
                ## print(paste("i: ", i))
                if (ent[i] >= ent_mean){
                    lborder <- i-1
                    v[reg,]$lborder = lborder
                    break
                }
            }

            for (i in v[reg,]$start:length(ent)){
                if (ent[i] >= ent_mean){
                    rborder <- i+1
                    v[reg,]$rborder = rborder
                    break
                }
            }

            v[reg,]$size <- (v[reg,]$stop-v[reg,]$lborder)+1
        }

        ##Remove regions with duplicated left-borders (Sometimes regions get caught multiple times)
        v <- v[!duplicated(v$lborder), ]
        ##Remove last index, this one is usually problematic
        v <- v[-length(v$lborder),]
    } else{
        v <- v[0,]
    }
    v
}

##Vecorizing by just using "subseq" or "substring" is not working.
##Using rep() on DNAString object is extending the string, not duplicating it into a StringSet.
##Would have to first convert to character, then use rep, etc.
sequence.getchars <- function(dna_seq, reg, offset=0){
    ## str_extracts <- sapply(1:nrow(reg), function(i){
    ##     subseq(dna_seq, reg$lborder[i]-offset, reg$stop[i]-offset)
    ## })
##    DNAStringSet(strngs)
    DNAStringSet(sapply(1:nrow(reg), function(i){
        subseq(dna_seq, reg$lborder[i]-offset, reg$stop[i]-offset)
    }))
}

#Give a BioStrings DNAStringSet Object
sequence.getfreq <- function(dna_str, k=7, prob=TRUE,  oligonuc=TRUE, order=FALSE){
    if(oligonuc){
        freq <- oligonucleotideFrequency(dna_str, k, as.prob=prob)
    } else {
        freq <- dinucleotideFrequency(dna_str, as.prob=prob)
    }
    if(order){
        freq <- freq[order(freq, decreasing=TRUE)]
    }
     freq
}

entropy.threshold <- function(v, sens=FALSE, sens_t=2, reg_ignore){
    ##Remove regions with N from vector before calculating quantiles
    b <- rep(TRUE, length(v))
    invisible(apply(reg_ignore, 1, function(x){
        b[x[1]:x[2]] <<- FALSE
    }))
    v <- v[b]
    if(sens){
        pr <- seq(0,0.001,0.001)
    }else{
        pr <- seq(0,0.01,0.01)
    }
    quant <- quantile(v, probs=pr)
    quant[sens_t]
}

entropy.calculate <- function(seq, k){
    ##seq <- DNAString(seq)
    ## f <- oligonucleotideFrequency(seq, width=k, as.prob=TRUE)
    ## f <- f[f>0]
    ## -sum(f * log2(f))

    f <- oligonucleotideFrequency(seq, width=k, as.prob=TRUE)
    ##f <- f[f>0]
    f.l <- log2(f)
    f.l[f==0] <- 0
    -rowSums(f*f.l)
}

fft.get <- function(seq, center=FALSE){
    seq_fft <- lapply(seq, function(x){
        fft(fft.encode(x, center=center))
       })
}

fft.get.single <- function(seq, center=FALSE){
    seq_enc <-  fft.encode(seq, center=center)
    seq_fft <- apply(seq_enc, 2, function(x){
        fft(x)
    })
    seq_fft
}

fft.peaks <- function(fft){
    fft_rsum <- rowSums(Mod(fft))
    
    ##fft_sma <- smoother(fft_rsum, n=20, plot=FALSE, method="supersmooth")
    ##fft_sma <- smoother(fft_sma, n=20, plot=FALSE)
    ##span <- 20/length(fft_rsum)
    ##fft_sma <- supsmu(x=1:length(fft_rsum), y=fft_rsum, span=span)
    ##fft_sma <- fft_sma$y

    ##Alternative smoothing, using convolve
    ##fft_sma <- convolve(fft_rsum, dnorm(-20:20, mean=0, sd=5), type='open'), type="l", col="red")

    fft_sma  <- convolve(fft_rsum, dnorm(-20:20, mean=0, sd=4), type='filter')
   
    
    ##Remove NA's by setting to 0
    fft_sma[which(is.na(fft_sma))] <- 0
    ##Remove first peak
    fft_sma[1:20] <- 0
    
    
    
    ##Loop over the whole graph, save highest values    
##    fft_median <- median(fft_sma)
##    fft_median <- fft_median+fft_median/10
    fft_mean <- mean(fft_sma)
    val_max <- 0
    index <- 0
    i_max <- vector()
    peaks_found = FALSE
    i <- 1
    fft_mean <- mean(fft_sma)
    l_pos <- vector()
    r_pos <- vector()
    while(!peaks_found){
        for(j in 1:length(fft_sma)){
            ##Exclude values around peak area by checking for saved indces for left/right position, so that peak is not caught multiple times
            if(val_max < fft_sma[j] && !j%in%i_max && !j%in%(unlist(Map(":", l_pos, r_pos)))){
                val_max <- fft_sma[j]
                index <- j
            }
        }
        ##When the caught peak is not 40% bigger than the mean of the values, set peaks_found to true and abandon peak search
        if((val_max-fft_mean)/fft_mean*100 <= 40){
            peaks_found <- TRUE
        } else{
            ##Add peak index to Array of indeces. Then check for left-&right end of peak to get the whole area
            val_max <- 0
            i_max[i] <- index
            
            ##check for left pos
            for(j in i_max[i]:1){
                ##Check for j+1 and j+2, in case there are small "bumps" in the peak. Probably not the best solution.
                if(fft_sma[j] >= fft_sma[j-1] || fft_sma[j] >= fft_sma[j-2]){
                ##if(fft_sma[j] >= fft_mean){
                ##if(fft_sma[j] >= fft_sma[j-1]){
                    l_pos[i] <- j
                } else{
                    break
                }
            }
            ##check for right pos
            for(j in i_max[i]:length(fft_sma)){
                ##Check for j+1 and j+2, in case there are small "bumps" in the peak. Probably not the best solution.
                if(fft_sma[j] >= fft_sma[j+1] || fft_sma[j] >= fft_sma[j+2]){
                ##if(fft_sma[j] >= fft_mean){
                ##if(fft_sma[j] >= fft_sma[j+1]){
                    r_pos[i] <- j                   
                } else{
                    break
                }
            }
            ##Remove NA's for arrays with border indices. Obscured peaks are soemtimes caught in multiple intersections
            ##Causing NA values to be generated whenever there is no left or right border to be found for the previous peak index. 
            l_pos <- l_pos[-(which(is.na(l_pos)))]
            r_pos <- r_pos[-(which(is.na(r_pos)))]
            i <- i+1
        }
              
    }
##    print(i_max)
    
    pos <- vector(length=length(i_max))
    
##    print(paste("l_pos:", l_pos, "\n"))
##    print(paste("r_pos:", r_pos, "\n"))

    l_pos <- l_pos+20
    r_pos <- r_pos+20
    
    unlist(Map(":", l_pos, r_pos))
}



fft.plot <- function(seq_fft, seq=NULL, reg=NULL, score=NULL, ent_score=NULL, col=1, q.col=NULL, l_max=250, peaks=NULL){
    fft_sum <- rowSums(Mod(seq_fft))
    plot(fft_sum, type="l", col=col)
    if(!is.null(q.col)){
        b <- (fft_sum >= quantile(fft_sum, q.col))
        points(which(b), fft_sum[b], col="red")
    }

    if(!is.null(peaks)){
        ##b <- (peaks==fft_sum)
        points(peaks, fft_sum[peaks], col="red")
    }

    if(!is.null(seq) && !is.null(reg) && !is.null(score) && !is.null(ent_score)){
        seq <- sequence.insertlines(seq, l_max=l_max)
        stats <- paste(c("Size:", reg$size, "\n", "Ent:", ent_score, "\n",  "Start:", reg$start, "\n",  "Stop:", reg$stop, "\n",
                         "k2:", score["k2"], "\n",  "k3:", score["k3"], "\n",  "k4:", score["k4"]), collapse="")
        mtext(text=stats, cex=0.4, side=3, line=0, adj=1, family="mono") 
        mtext(text=seq, cex=0.3, side=3, line=0, adj=0, col="red", family="mono")
    }
}


kmer.score <- function(k1, k2){
    score <- sqrt( sum( (k1-k2)^2 ))
    score
}

repeat.freq <- function(k){
    f <- rep(0, 4^k)
    f[1:k] <- 1/k
    f
}

region.getscores <- function(reg, k){
    scores <- vector(mode='numeric', length=k-1)
    k_names <- vector(mode='character', length=k-1)
    for (i in 2:k){
         freq <- oligonucleotideFrequency(reg, as.prob=TRUE, width=i)
         kmer_freq <- repeat.freq(i)
         k_names[i-1] <- paste0("k",i)
         freq <- freq[order(freq, decreasing=TRUE)]
         kmer_freq <- kmer_freq[order(kmer_freq, decreasing=TRUE)]
         scores[i-1] <- kmer.score(kmer_freq[1:i], freq[1:i])
    }
    names(scores) <- k_names
    scores
}

region.summarize <- function(strngs, regions){
    data <- sapply(strngs, region.getscores, k=2)
    cbind( regions[,c("size", "lborder", "stop")], t(data))
}

region.smoothing <- function(ent){
    sma <- stats::filter(ent[[1]], filter = rep(1, 10 * 12) / (12*10), circular=TRUE)
    sma
}

plot.sizefreq <- function(reg_data, k, organism_name){
    pdf(paste0(organism_name,".pdf"))
        for(i in 1:length(reg_data)){
           plot(reg_data[[i]]$regions$size, reg_data[[i]]$scores[,k-1], xlab="size", ylab=paste(paste0("k",k),"score"), main=paste(paste0("k",k), organism_name))
        }
    dev.off()
}

plot.drawseq <- function(seq, reg, ent, offset=1000){
    length_max <- 150

    seq_b <- substr(seq, reg$lborder-offset, reg$lborder)
    seq_lcomp <- substr(seq, reg$lborder, reg$stop)

    seq_b <- sequence.insertlines(seq_b)
    seq_lcomp <- sequence.insertlines(seq_lcomp)
    
    posx_b <- c(reg$lborder-offset, reg$lborder)
    posy_b <- ent[posx_b]
    points(x=posx_b, y=posy_b, col="blue", pch=19)

    posx_lcomp <- c(reg$lborder, reg$stop)
    posy_lcomp <- ent[posx_lcomp]
    points(x=posx_lcomp, y=posy_lcomp, col="red", pch=19)

    mtext(text=seq_lcomp, cex=0.3, side=3, line=0, adj=1, col="red", family="mono")
    mtext(text=seq_b, cex=0.3, side=3, line=0, adj=0.3, col="blue", family="mono")
}

sequence.insertlines <- function(seq, l_max=150){
    if (nchar(seq) > l_max){
        start <- seq(0, nchar(seq), l_max)
        stop <- seq(l_max, nchar(seq), l_max)
        start <- start[-length(start)]
        
        ## Catch the rest of the string
        start <- c(start, stop[length(stop)])
        stop <- c(stop, nchar(seq))
        
        seq <- substring(seq, start, stop)
        seq <- paste(seq, collapse="\n")
    }
    seq 
}



rand.nseq <- function(length, f=rep(0.25, 4)){
    paste(sample(c("A", "C", "G", "T"), size=length, replace=TRUE, prob=f), collapse="")
}

kmer.rep <- function(kmer, length){
    paste(rep(kmer, length), collapse="")
}

mutate.substitute <- function(seq, f=rep(0.25,4)){
    alphabet <- c("A", "T", "C", "G")
    substitute <- sample(alphabet, size=1, prob=f)
    pos <- sample(1:nchar(seq), size=1)
    replaceLetterAt(seq, pos, substitute)
}

mutate.insert <- function(seq, f=rep(0.25, 4)){
    alphabet <- c("A", "T", "C", "G")
    insert <- sample(alphabet, size=1,  prob=f)
    pos <- sample(1:nchar(seq), size=1)
    replaceAt(seq, pos, insert)
}

mutate.delete <- function(seq){
    pos <- sample(1:nchar(seq), size=1)
    subseq(seq, pos, pos) <- NULL
    seq
}

seq.mutate <- function(seq, prob_s=0, prob_i=0, prob_d=0){
    prob <- seq(0.001,1,0.001)
    for(i in 1:nchar(seq)){
        if(sample(prob, size=1) <= prob_s){
##            print("SUBSTITUTED")
            seq <- mutate.substitute(seq)
        }
        if(sample(prob, size=1) <= prob_i){
##            print("INSERTED")
##            print(paste0("\n", "prob:", prob_i, "\n"))
            seq <- mutate.insert(seq)
        }
        if(sample(prob, size=1) <= prob_d){
##            print("DELETED")
##            print(paste0("\n", "prob:", prob_d, "\n"))
            seq <- mutate.delete(seq)
        }
    }
    seq
}


## takes a normal character vector
## ins and del, make single nucleotide insertions,
## but it would probably be better to sample from
## some distribution to get the length of the indels
## seq.mutate <- function(seq, sub.f=0.01, ins.f=sub.f, del.f=ins.f){
##     n <- nchar(seq)
##     sub.pos <- sample(1:nchar(seq), sub.f*n)
##     ins.pos <- sort(sample(3:nchar(seq) - 1, ins.f*n))
##     del.pos <- sort(sample(3:nchar(seq) - 1, del.f*n))
##     nucs <- c("A", "C", "T", "G")
##     sub.n <- sample(nucs, length(sub.pos), replace=TRUE)
##     ins.n <- sample(nucs, length(ins.pos), replace=TRUE)
##     if(length(sub.pos)){
##         seq <- paste(rbind( substring( seq, c(1, sub.pos+1), c(sub.pos-1, nchar(seq)) ),
##                            c(sub.n, "")), collapse="")
##     }
##     seq.sp <- substring( seq, c(c(1, ins.pos+1)), c(ins.pos, nchar(seq)))
##     seq <- paste( rbind(seq.sp, c(ins.n, "")), collapse="" )
##     paste( substring( seq, c(1, del.pos+1), c(del.pos-1, nchar(seq)) ), collapse="")
## }

rev.fft <- function(ft, qnt=0.9, do.plot=TRUE, fft_single=FALSE){
    m <- rowSums(Mod(ft))
    ##m[1] <- 0
    q <- quantile(m, qnt)
    b <- m <= q
    ft[ b, ] <- 0
    ft[q, ] <- 0
##  f2 <- fft( ft, inverse=TRUE )
    if(fft_single){
        ## f1 <- fft(ft[,1], inverse=TRUE)
        ## f2 <- fft(ft[,2], inverse=TRUE)
        ## f3 <- fft(ft[,3], inverse=TRUE)
        ## f4 <- fft(ft[,4], inverse=TRUE)

        fft_r <- apply(ft, 2, function(x){
            fft(x, inverse=TRUE)
        })
        
    } else{
        fft_r <- fft(ft, inverse=TRUE)
    }
    r2 <- Re(fft_r) ## should give the same as Mod, as inverse fft?
    nucs <- colnames(ft)
    nucs[ apply(r2, 1, which.max) ]
}

## rev.fft.peaks <- function(ft, peaks){
##     m <- rowSums(Mod(ft))
##     peaks <- m[peaks]
##     b <- m!=peaks
##     ft[b, ] <- 0
##     f2 <- fft(ft, inverse=TRUE)
##     r2 <- Re(f2)
##     nucs <- colnames(ft)
##     nucs[ apply(r2, 1, which.max) ]
## }

rev.fft.peaks <- function(ft, peaks){
    ft[-peaks, ] <- 0
    f2 <- fft(ft, inverse=TRUE)
    r2 <- Re(f2)
    nucs <- colnames(ft)
    nucs[ apply(r2, 1, which.max) ]
}


fft.get.repeat <- function(seq){
    seq_fft <- fft.get(seq)
    mod_fft <- Mod(seq_fft[[1]])
    thrsh <- quantile(mod_fft, 0.99)
    
    seq_fft[[1]][mod_fft <= thrsh] = 0
    fft_r <- fft(seq_fft[[1]], inverse=TRUE)
    #fft_r[mod_fft < thrsh] = 0
    seq_rec <- colnames(fft_r)[apply(Re(fft_r), 1, which.max)]

    seq_rec
}
























