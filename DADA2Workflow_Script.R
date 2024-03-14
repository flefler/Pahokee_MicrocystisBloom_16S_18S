######## PACKAGES & SUCH ###########

#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("dada2")
#BiocManager::install("DECIPHER")
#BiocManager::install("Biostrings")

#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#devtools::install_github("jbisanz/qiime2R")
#install.packages('ggthemes', dependencies = TRUE)
#install.packages('gg_ordiplot')
#install.packages("remotes")
#remotes::install_github("microbiome/microbiome")
#remotes::install_github("gauravsk/ranacapa")
#devtools::install_github("david-barnett/microViz")
#install.packages("remotes")
#remotes::install_github("jfq3/ggordiplots")
library(readxl)
library(dada2)
library(microViz)
library(gclus)
library(ggrepel)
library(ggforce)
library(mixOmics)
library(BiodiversityR)# also loads vegan
library(ranacapa)
library(microbiome)
library(tidyverse)
library(phyloseq)
library(ggsci)
library(ggthemes)
library(ade4)
library(adegraphics)
library(adespatial)
library(vegan3d)
library(MASS)
library(ellipse)
library(FactoMineR)
library(rrcov)
library(missMDA)
library(Amelia)
library(BiodiversityR)
library(ggordiplots)
library(DECIPHER)
library(phyloseq)
library(Biostrings)
library(ggplot2)
library("ape")
###########################################################################
#####IMPORT READS
path <- "~/Dropbox (UFL)/Laughinghouse_Lab/PROJECTS/West_Pahokee/Pahokee_Bacteria" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

#####BEGIN HERE TO EXAMINE QUALITY OF READS

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Inspect read quality profiles
plotQualityProfile(fnFs[1:2]) #Forward Reads

plotQualityProfile(fnRs[1:2]) #Reverse Reads

#####BEGIN HERE FOR DADA2 PROCESS
###Filter and Trim

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Filter forward and reverse reads. TruncLen is determined by inspecting quality profiles and is dependent on sequences.
#In this example, forward reads are truncated at bp 285 and reverse reads are truncated at bp 275.
#maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE and maxEE=2 are the default parameters.
#The maxEE parameter sets the maximum number of "expected errors" allowed in a read, which is a better filter than simply averaging quality scores.
#Primers
#FWD <- "GTGYCAGCMGCCGCGGTAA"  ## CHANGE ME # this is 515F
#REV <- "CCGYCAATTYMTTTRAGTTT"  ## CHANGE ME # this is 926YR


#the data looks great, therefore were are not going to trun anything
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(225,221),
                     trimLeft = c(19,20), #These are the lengths of the forward and reverse primers
                     maxN=0, maxEE=c(1,1), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, verbose = TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)


#####Determining Error Rates#####
#Four options for learning error rates with NovaSeq data
#https://github.com/ErnakovichLab/dada2_ernakovichlab/tree/split_for_premise

#Option 1 from JacobRPrice alter loess arguments (weights and span and enforce monotonicity) benjjneb/dada2#1307
loessErrfun_mod1 <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # Gulliem Salazar's solution
        # https://github.com/benjjneb/dada2/issues/938
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

# check what this looks like
errF_1 <- learnErrors(
  filtFs,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod1,
  verbose = TRUE
)
errR_1 <- learnErrors(
  filtRs,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod1,
  verbose = TRUE
)


#Option 2 enforce monotonicity only. Originally recommended in: benjjneb/dada2#791
loessErrfun_mod2 <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # Gulliem Salazar's solution
        # https://github.com/benjjneb/dada2/issues/938
        # mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}


# check what this looks like
errF_2 <- learnErrors(
  filtFs,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod2,
  verbose = TRUE
)
errR_2 <- learnErrors(
  filtRs,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod2,
  verbose = TRUE
)

#Option 3 alter loess function (weights only) and enforce monotonicity From JacobRPrice benjjneb/dada2#1307
loessErrfun_mod3 <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # Gulliem Salazar's solution
        # https://github.com/benjjneb/dada2/issues/938
        # mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        
        # only change the weights
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot))
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

# check what this looks like
errF_3 <- learnErrors(
  filtFs,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod3,
  verbose = TRUE
)
# check what this looks like
errR_3 <- learnErrors(
  filtRs,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod3,
  verbose = TRUE
)

#Option 4 Alter loess function arguments (weights and span and degree, also enforce monotonicity) From Jonalim’s comment in benjjneb/dada2#1307

loessErrfun_mod4 <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # jonalim's solution
        # https://github.com/benjjneb/dada2/issues/938
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),degree = 1, span = 0.95)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

# check what this looks like
errF_4 <- learnErrors(
  filtFs,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod4,
  verbose = TRUE
)
errR_4 <- learnErrors(
  filtRs,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod4,
  verbose = TRUE
)

#There are four options, none of which will yield “ideal” error plots. Instead look for the solution where the black line is continuously decreasing 
#(i.e. as quality scores improve on the x-axis the predicted error rate (y-axis) goes down) 
#and for plots that have points that mostly align with the black lines, although you will likely have some points along 0 on the y-axis.


# Original default recommended way (not optimal for NovaSeq data!)
errF_plot <- plotErrors(errF, nominalQ = TRUE)
errR_plot <- plotErrors(errR, nominalQ = TRUE)

errF_plot
errR_plot

# Trial 1 (alter span and weight in loess, enforce montonicity)
errF_plot1 <- plotErrors(errF_1, nominalQ = TRUE)
errR_plot1 <-plotErrors(errR_1, nominalQ = TRUE)

errF_plot1
errR_plot1

# Trial 2 (only enforce monotonicity - don't change the loess function)
errF_plot2 <- plotErrors(errF_2, nominalQ = TRUE)
errR_plot2 <-plotErrors(errR_2, nominalQ = TRUE)

errF_plot2
errR_plot2

# Trial 3 (alter loess (weights only) and enforce monotonicity)
errF_plot3 <- plotErrors(errF_3, nominalQ = TRUE)
errR_plot3 <-plotErrors(errR_3, nominalQ = TRUE)

errF_plot3
errR_plot3

# Trial 4 (alter loess (span, weight, and degree) and enforce monotonicity)
errF_plot4 <- plotErrors(errF_4, nominalQ = TRUE)
errR_plot4 <-plotErrors(errR_4, nominalQ = TRUE)

errF_plot4
errR_plot4









#Applying sample inference algorithm
dadaFs <- dada(filtFs, err=errF_4, multithread=TRUE, pool=TRUE) #Forward Reads
dadaRs <- dada(filtRs, err=errR_4, multithread=TRUE, pool=TRUE) #Reverse Reads


dadaFs[[1]] #Identifies amount of sequence variants from forward reads.
dadaRs[[1]] #Identifies amount of sequence variants from reverse reads.

###Merge paired reads

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

###Construct sequence tables

seqtab <- makeSequenceTable(mergers) #Creates amplicon sequence table
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

##Optional: remove non-target-length sequences from your sequence table

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 381:394] #get amplicons of the targeted length

###Remove chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE) #combines a left-segment and a right-segment from two more abundant "parent" sequences.
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab2) #determines the amount of chimeras in merged sequence reads. 1-n = %chimeras.


###Track reads through the pipeline. Determine the number of reads that made it through each step in the pipeline

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
track

##Note: Outside of filtering, there should be no step in which a majority of reads are lost.
##If too many reads were lost return to filtering step.

#####END OF DADA2 PROCESSING
########################################################################
#####ASSIGN TAXONOMY
###IDTAXA Classifier
library(phyloseq)
library(tidyverse)
#remotes::install_github("vmikk/metagMisc")

#Format ASV Table

mertab_collapse = collapseNoMismatch(seqtab.nochim)

taxa <- dada2::assignTaxonomy(mertab_collapse, minBoot = 80, "~/Dropbox (UFL)/Laughinghouse_Lab/CyanoSeq/v_1.2/CyanoSeq1.2_BLCC_SILVA138.1_dada2.fastq.gz", multithread = T, verbose = F)

Pahokee_Metadata <- read.csv("~/Dropbox (UFL)/Laughinghouse_Lab/PROJECTS/West_Pahokee/Metadata/Pahokee_Metadata.csv", row.names = 1)

ps_SilvaCyanoseq <- phyloseq(otu_table(mertab_collapse, taxa_are_rows=FALSE), 
                             sample_data(Pahokee_Metadata),
                             tax_table(taxa))

ps_SilvaCyanoseq_clean = ps_SilvaCyanoseq %>%
  subset_taxa((Class!="Chloroplast") | is.na(Class)) %>%
  subset_taxa((Class !="Cyanobacteriia") | is.na(Class)) %>%
  subset_taxa((Kingdom != "Eukaryota") | is.na(Kingdom)) %>%
  subset_taxa((Family != "Mitochondria") | is.na(Family)) %>%
  subset_taxa((Phylum != "NA") | is.na(Phylum))

ps_SilvaCyanoseq_clean <- microbiomeutilities::add_refseq(ps_SilvaCyanoseq_clean,tag="ASV")
#8759 taxa

otutable = (ps_SilvaCyanoseq_clean@otu_table) %>% as.data.frame()
taxtable = (ps_SilvaCyanoseq_clean@tax_table) %>% as.data.frame()

vegan::rarecurve(otutable, step=50, cex=0.5)

ps_rarefied <- rarefy_even_depth(ps_SilvaCyanoseq_clean,
                                 rngseed = 42069,
                                 sample.size = min(sample_sums(ps_SilvaCyanoseq_clean)),
                                 replace = FALSE)
ps_rarefied
#8583 taxa


#This is where i ended december 18th 2023
saveRDS(ps_rarefied, "~/Dropbox (UFL)/Laughinghouse_Lab/PROJECTS/West_Pahokee/Pahokee_Jan19.rds")
save.image(file="~/Dropbox (UFL)/Laughinghouse_Lab/PROJECTS/West_Pahokee/Pahokee_Jan19th.RData") 


#####Below this is trash####

#Filter out low abundance taxa
ps_bac_filt = tax_filter(ps_SilvaCyanoseq_clean,
                         prev_detection_threshold = .2, #ASV's must occur >0.1% of total realtavie abundance
                         #min_prevalence =  #must occur in at least TWO samples we have 49 samples, 2/49 = 0.04, not using this since so many diverse environments 
                         min_total_abundance = 100 #THE ASV MUST OCCUR AT LEAST 100 TIMES IN TOTAL ACROSS ALL SAMPLES
)
ps_bac_filt

q = as.data.frame(ps_bac_filt@tax_table)

#Reduction from 9678 to 1184 ASVs

saveRDS(ps_bac_filt, "~/Dropbox (UFL)/Laughinghouse_Lab/PROJECTS/West_Pahokee/Pahokee_Aug24th.rds")
save.image(file="~/Dropbox (UFL)/Laughinghouse_Lab/PROJECTS/West_Pahokee/Pahokee_Aug24th.RData") 

load("/Users/flefler/Dropbox (UFL)/Laughinghouse_Lab/PROJECTS/West_Pahokee/Pahokee_Aug24th.RData")



##BACTERIA
ps_bac = ps_SilvaCyanoseq_clean %>%
  subset_taxa((Class !="Cyanophyceae") | is.na(Class))

saveRDS(ps_bac, "~/Dropbox (UFL)/Laughinghouse_Lab/PROJECTS/West_Pahokee/Pahokee_Baconly_unfiltered.rds")


#Filter out low abundance taxa
ps_bac_filt = tax_filter(ps_bac,
                         prev_detection_threshold = .2 #ASV's must occur >0.1% of total realtavie abundance
                         #min_prevalence =  #must occur in at least TWO samples we have 49 samples, 2/49 = 0.04, not using this since so many diverse environments 
                         #min_total_abundance = 100 #THE ASV MUST OCCUR AT LEAST 100 TIMES IN TOTAL ACROSS ALL SAMPLES
)

#Reduction from 1304 to 1067 ASVs

library(dada2)

#Steps to make the phylogenetic tree
#Very computationally intensive

library(phangorn)

alignment <- AlignSeqs(DNAStringSet(ps_bac_filt@refseq), anchor=NA) 
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))



#function for picking outgroup/Root https://john-quensen.com/r/unifrac-and-tree-roots/
pick_new_outgroup <- function(tree.unrooted){
  require("magrittr")
  require("data.table")
  require("ape") # ape::Ntip
  # tablify parts of tree that we need.
  treeDT <-
    cbind(
      data.table(tree.unrooted$edge),
      data.table(length = tree.unrooted$edge.length)
    )[1:Ntip(tree.unrooted)] %>%
    cbind(data.table(id = tree.unrooted$tip.label))
  # Take the longest terminal branch as outgroup
  new.outgroup <- treeDT[which.max(length)]$id
  return(new.outgroup) }


my.tree <- phy_tree(fitGTR$tree)
out.group <- pick_new_outgroup(my.tree)
out.group ## [1] ASV441 Look and see what it is in the ASV fasta file

new.tree <- ape::root(my.tree, outgroup=out.group, resolve.root=TRUE)
ps_bac_filt@phy_tree <- phy_tree(new.tree)
phy_tree(ps_bac_filt)

saveRDS(ps_bac_filt, "~/Dropbox (UFL)/Laughinghouse_Lab/PROJECTS/West_Pahokee/Pahokee_Bac_only_wTree.rds")

##CYANOBACTERIA
ps_cyano = ps_SilvaCyanoseq_clean %>%
  subset_taxa((Class =="Cyanophyceae") | is.na(Class))

saveRDS(ps_cyano, "~/Dropbox (UFL)/Laughinghouse_Lab/PROJECTS/West_Pahokee/Pahokee_Baconly_unfiltered.rds")


#Filter out low abundance taxa
ps_cyano_filt = tax_filter(ps_cyano,
                         prev_detection_threshold = .2 #ASV's must occur >0.1% of total realtavie abundance
                         #min_prevalence =  #must occur in at least TWO samples we have 49 samples, 2/49 = 0.04, not using this since so many diverse environments 
                         #min_total_abundance = 100 #THE ASV MUST OCCUR AT LEAST 100 TIMES IN TOTAL ACROSS ALL SAMPLES
)

#Reduction from 1304 to 1067 ASVs

library(dada2)

#Steps to make the phylogenetic tree
#Very computationally intensive

library(phangorn)

alignment <- AlignSeqs(DNAStringSet(ps_cyano_filt@refseq), anchor=NA) 
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))



#function for picking outgroup/Root https://john-quensen.com/r/unifrac-and-tree-roots/
pick_new_outgroup <- function(tree.unrooted){
  require("magrittr")
  require("data.table")
  require("ape") # ape::Ntip
  # tablify parts of tree that we need.
  treeDT <-
    cbind(
      data.table(tree.unrooted$edge),
      data.table(length = tree.unrooted$edge.length)
    )[1:Ntip(tree.unrooted)] %>%
    cbind(data.table(id = tree.unrooted$tip.label))
  # Take the longest terminal branch as outgroup
  new.outgroup <- treeDT[which.max(length)]$id
  return(new.outgroup) }


my.tree <- phy_tree(fitGTR$tree)
out.group <- pick_new_outgroup(my.tree)
out.group ## [1] ASV441 Look and see what it is in the ASV fasta file

new.tree <- ape::root(my.tree, outgroup=out.group, resolve.root=TRUE)
ps_cyano_filt@phy_tree <- phy_tree(new.tree)
phy_tree(ps_cyano_filt)

saveRDS(ps_cyano_filt, "~/Dropbox (UFL)/Laughinghouse_Lab/PROJECTS/West_Pahokee/Pahokee_cyano_only_wTree.rds")




save.image(file="~/Dropbox (UFL)/Laughinghouse_Lab/PROJECTS/West_Pahokee/Pahokee_Bac.RData") 










