

# Analysing one primer pairs in a sequencing run --------------------------

#This will be an example script when two different pairs of primers are in one run. You can use the first half of the script if only one primer pair is used

# Loading data into R and dada2 -------------------------------------------

#change directory to home
cd ~ 
  #make a directory in home called "name?"
mkdir name?
#change directory to "name?"
cd name?


# Pipeline analysis for Metabardcoding - packages -------------------------

library(dada2)
library(Biostrings)
library(ShortRead)
library(phyloseq)
library(vegan)
library(ggplot2)


# Loading data --------------------------------------------------------------------

input.path <- #where data is
output.path <- #working location in Linux (name?)

list.files(input.path) #should see a list of paired fastq sequencing files


# Inputting the forward and reverse reads ---------------------------------

fnFs <- sort(list.files(input.path, pattern = "file name_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(input.path, pattern = "file name_R2.fastq", full.names = TRUE))


# Identifying and removing primer -----------------------------------------

#Primer sequences 
FWD <- "sequence"
REV <- "sequence"

#Primer orientation checking
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString 
  #objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

FWD.orients <- allOrients(FWD)
FWD.orients # print all orientations of the forward primer to the console
REV.orients <- allOrients(REV)
REV.orients # print all orientations of the reverse primer to the console

fnFs.filtN <- file.path(output.path, "filtN", basename(fnFs)) 
fnRs.filtN <- file.path(output.path, "filtN", basename(fnRs))

filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

#For this function we have specified:
  
#the path to our raw forward read fatsqs fnFs
#the path to write our N-trimmed forward fastqs fnFs.filtN
#the path to our raw reverse read fatsqs fnRs
#the path to write our N-trimmed reverse fastqs fnRs.filtN
#and the function-specific options for maxN and multithread

#The following function will use all possible primer combinations to count the number of times a primer is found in the forward and reverse read in each orientation.

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

#[[1]] runs the function on the first file in the list/this number can be increase
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]])) 

# Cutadapt ----------------------------------------------------------------

#Specify the path to the cutadapt software so R knows where to find and run it 

cutadapt <- "location of cutadapt"

#Specofy the path to the output files where we want to put the cutadapt output files 

path.cut <- file.path(output.path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

#Specify the options needed by cutadapt

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC)

#In addition to trimming the primers from the reads we will also specify a couple of extra useful parameters.
#--discard-untrimmed this tells cutadapt to discard any read where the primers haven't been trimmed off. 
#This is especially important for our data as our files at the moment contain sequences amplified using both MiFish-U and 12S-V5 primer sets.
#We only want to keep sequences matching the MiFish-U primer set for this analysis.--minimum-length 60 discard reads shorter than 60bp. 
#This will remove unexpected short reads and help speed up further analysis.

for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, 
                             # -n 2 required to remove FWD and REV
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i], # input files
                             "--discard-untrimmed",
                             "--minimum-length 60"))
}

#Check whether all the primers have been removed using the primerHits function specified earlier

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

#matric should show all 0's - meaning we have no primers remaining in the file 


# Checking the quality of your data ---------------------------------------

#To run first import the cutadapt files and extract the sample names and then run the plotQualityProfile function

# Specify the paths and file names the forward and reverse primer cleaned files 
cutFs <- sort(list.files(path.cut, pattern = "_L001_R1_001.fastq", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_L001_R2_001.fastq", full.names = TRUE))

# Extract sample names
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][2]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

# check the quality for the first file
plotQualityProfile(cutFs[1:1])

#check the quality for the reverse file 
plotQualityProfile(cutRs[1:1])

#to check the quality of additional fastq files i.e second and third 
plotQualityProfile(cutFs[2:3])

plotQualityProfile(cutRs[2:3])


# Cleaning the data -------------------------------------------------------

#Set the path to a directory to store the filtered output files called "filtered
filtFs <- file.path(path.cut, "../filtered", basename(cutFs))
filtRs <- file.path(path.cut, "../filtered", basename(cutRs))

#Now run filterAndTrim. This time we use the standard filtering parameters:
#maxN=0 After truncation, sequences with more than 0 Ns will be discarded. (DADA2 requires sequences contain no Ns)
#truncQ = 2 Truncate reads at the first instance of a quality score less than or equal to 2
#rm.phix = TRUE Discard reads that match against the phiX genome
#maxEE=c(2, 2) After truncation, reads with higher than 2 "expected errors" will be discarded
#minLen = 60 Remove reads with length less than 60 (note these should have already been removed by cutadapt)
#multithread = TRUE input files are filtered in parallel

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), 
                     truncQ = 2, minLen = 60, rm.phix = TRUE, compress = TRUE, 
                     multithread = TRUE)
out


# Identification of ASVs & Generate an error model ------------------------

#Model the error rates using both the forward and reverse reads 
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

#plot errors
plotErrors(errF, nominalQ = TRUE)

#The error rates for each possible transition (e.g. A→C, A→G) are shown
#Red line - expected based on the quality score. (These are plotted when nominalQ = TRUE is included in the plot command)
#Black line - estimate
#Black dots - observed

#Dereplication 
exists <- file.exists(filtFs) 
# check that all the samples are still present after filtering
derepFs <- derepFastq(filtFs[exists], verbose=TRUE)
derepRs <- derepFastq(filtRs[exists], verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names[exists]
names(derepRs) <- sample.names[exists]

#Inferrence of ASVs
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)
#Important that you want to run the error and inference steps on datasets generated from a single Illumina run,
#different runs can have different error structures


# Merging paired end reads ------------------------------------------------

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

#Making our ASV matrix 
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#Chimera detection and removal
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                    multithread=TRUE, verbose=TRUE)

dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

#range of ASV lengths
table(nchar(getSequences(seqtab.nochim)))

#Sequence tracking sanity check 
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), 
               sapply(dadaRs, getN), 
               sapply(mergers, getN), 
               rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", 
                     "denoisedF", "denoisedR", 
                     "merged", "nonchim")
rownames(track) <- sample.names
track


# Assigning taxonomy ------------------------------------------------------

taxa <- assignTaxonomy(seqtab.nochim, 
                       ".fasta", #fasta reference database i.e. /pub14/tea/nsc206/NEOF/metabarcoding_workshop/taxonomy/12S-V5_Reference_Database_taxonomy.fasta
                       multithread=TRUE, verbose = T)

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)


# Further analysis  -------------------------------------------------------

#Load table which contains meta data
meta<-read.csv(".csv", 
               row.names = 1)

#How to remove samples from data i.e. too little reads, negative controls etc
seqtab.rmcontrol<-seqtab.nochim[-c(25,31),] #25 and 31 are sample names - these will change with your data
meta.rmcontrol<-meta[-c(25,31),]

#Remove sameples with failed/low sequence numbers - all examples of sample IDs and read depth
#Check number of sequences per sample
rowSums(seqtab.rmcontrol)
# S40, S41, S50 all have low sequence counts. So we will remove these rows.
seqtab.rmlow<-seqtab.rmcontrol[-c(13,14,23),]
meta.rmlow<-meta.rmcontrol[-c(13,14,23),]

# Print the minimum sequence number in one sample.
min(rowSums(seqtab.rmlow))
# The lowest number of sequences in one sample is now 26667.


# Rarefaction curves ------------------------------------------------------

rarecurve(seqtab.rmlow, step=100, col=meta.rmlow$COLOUR, lwd=2, ylab="ASVs", label=F)
# add a vertical line to represent the fewest sequences in any sample
abline(v=(min(rowSums(seqtab.rmlow))))


# Alpha diversity ---------------------------------------------------------

seqtab.rmlow2<-t(as.data.frame(seqtab.rmlow))
phylo_asv <- otu_table(seqtab.rmlow2, taxa_are_rows=TRUE)
phylo_tax <- tax_table(taxa)
phylo_samples <- sample_data(meta.rmlow)

phylo <- phyloseq(phylo_asv, phylo_tax, phylo_samples)

sample_names(phylo)
rank_names(phylo)
sample_variables(phylo)

#Two alpha diversity metrics to plot - Shannon's diversity index and Simpson diversity index 
plot_richness(phylo, 
              measures=c("Shannon", "Simpson"), 
              color = "SITE")
plot_richness(phylo, x="SITE", measures=c("Shannon", "Simpson"), 
              color = "SITE") + geom_boxplot()


# Beta diversity ----------------------------------------------------------

ps.prop <- transform_sample_counts(phylo, function(otu) otu/sum(otu))

ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="SITE", title="Bray NMDS")

bray.dist<-distance(ps.prop, method="bray")
sampledf <- data.frame(sample_data(phylo))
adonis2(bray.dist ~ SITE, data = sampledf)

plot_bar(ps.prop, fill = "Family")+
  geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack")+
  facet_grid(~SITE, scales = "free", space = "free")


# Analysing a second primer pair ----------------------------------------

#It can be common to use multiple genes or multiple markers of the same gene in metabarcoding studies, depending on taxa
#being targetted, multiple primer pairs can complement each other by increasing taxonomic coverage or resolution 

# Loading data 2 --------------------------------------------------------------------

#these maybe in the same location
input.path <- #where data is
output.path <- #working location in Linux (name?)
  
  list.files(input.path) #should see a list of paired fastq sequencing files


# Inputting the forward and reverse reads 2 ---------------------------------

#maybe the same reads just looking at a smaller amplicon
fnFs <- sort(list.files(input.path, pattern = "file name_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(input.path, pattern = "file name_R2.fastq", full.names = TRUE))


# Identifying and removing primer 2 -----------------------------------------

#Primer sequences of second pair
FWD2 <- "sequence"
REV2 <- "sequence"

#Primer orientation checking
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString 
  #objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}


FWD2.orients <- allOrients(FWD2)
FWD2.orients # print all orientations of the forward primer to the console
REV2.orients <- allOrients(REV2)
REV2.orients # print all orientations of the reverse primer to the console

fnFs.filtN <- file.path(output.path, "filtN", basename(fnFs)) 
fnRs.filtN <- file.path(output.path, "filtN", basename(fnRs))

filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

#The following function will use all possible primer combinations to count the number of times a primer is found in the forward and reverse read in each orientation.

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD2.ForwardReads = sapply(FWD2.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD2.ReverseReads = sapply(FWD2.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV2.ForwardReads = sapply(REV2.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV2.ReverseReads = sapply(REV2.orients, primerHits, fn = fnRs.filtN[[1]]))

#[[5]] runs the function on the first five files in the list/this number can be increased further
rbind(FWD2.ForwardReads = sapply(FWD2.orients, primerHits, fn = fnFs.filtN[[5]]), 
      FWD2.ReverseReads = sapply(FWD2.orients, primerHits, fn = fnRs.filtN[[5]]), 
      REV2.ForwardReads = sapply(REV2.orients, primerHits, fn = fnFs.filtN[[5]]), 
      REV2.ReverseReads = sapply(REV2.orients, primerHits, fn = fnRs.filtN[[5]]))

#Specify the path to the cutadapt software so R knows where to find and run it 

cutadapt <- "location of cutadapt"

#Specofy the path to the output files where we want to put the cutadapt output files - cut2 specifies the second primer pair

path.cut2 <- file.path(output.path, "cutadapt2")
if(!dir.exists(path.cut2)) dir.create(path.cut2)
fnFs.cut2 <- file.path(path.cut2, basename(fnFs))
fnRs.cut2 <- file.path(path.cut2, basename(fnRs))

#Specify the options needed by cutadapt

FWD2.RC <- dada2:::rc(FWD2)
REV2.RC <- dada2:::rc(REV2)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags2 <- paste("-g", FWD, "-a", REV2.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags2 <- paste("-G", REV, "-A", FWD2.RC)

#In addition to trimming the primers from the reads we will also specify a couple of extra useful parameters.
#--discard-untrimmed this tells cutadapt to discard any read where the primers haven't been trimmed off. 
#This is especially important for our data as our files at the moment contain sequences amplified using both MiFish-U and 12S-V5 primer sets.
#We only want to keep sequences matching the MiFish-U primer set for this analysis.--minimum-length 60 discard reads shorter than 60bp. 
#This will remove unexpected short reads and help speed up further analysis.

for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags2, R2.flags2, "-n", 2, 
                             # -n 2 required to remove FWD and REV
                             "-o", fnFs.cut2[i], "-p", fnRs.cut2[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i], # input files
                             "--discard-untrimmed",
                             "--minimum-length 60"))
}

#Check whether all the primers have been removed using the primerHits function specified earlier

rbind(FWD2.ForwardReads = sapply(FWD2.orients, primerHits, fn = fnFs.cut2[[1]]), 
      FWD2.ReverseReads = sapply(FWD2.orients, primerHits, fn = fnRs.cut2[[1]]), 
      REV2.ForwardReads = sapply(REV2.orients, primerHits, fn = fnFs.cut2[[1]]), 
      REV2.ReverseReads = sapply(REV2.orients, primerHits, fn = fnRs.cut2[[1]]))

#matric should show all 0's - meaning we have no primers remaining in the file 


# Checking the quality of your data 2 ---------------------------------------

#To run first import the cutadapt files and extract the sample names and then run the plotQualityProfile function

# Specify the paths and file names the forward and reverse primer cleaned files 
cutFs2 <- sort(list.files(path.cut2, 
                          pattern = "_L001_R1_001.fastq", full.names = TRUE))
cutRs2 <- sort(list.files(path.cut2, 
                          pattern = "_L001_R2_001.fastq", full.names = TRUE))

# Extract sample names
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][2]
sample.names <- unname(sapply(cutFs2, get.sample.name))
head(sample.names)

# check the quality for the first file
plotQualityProfile(cutFs2[1:1])

#check the quality for the reverse file 
plotQualityProfile(cutRs2[1:1])

#to check the quality of additional fastq files i.e second and third 
plotQualityProfile(cutFs2[2:3])

plotQualityProfile(cutRs2[2:3])


# Cleaning the data 2 -------------------------------------------------------

#Set the path to a directory to store the filtered output files called "filtered
filtFs2 <- file.path(path.cut2, "../filtered2", basename(cutFs2))
filtRs2 <- file.path(path.cut2, "../filtered2", basename(cutRs2))

#Now run filterAndTrim. This time we use the standard filtering parameters:
#maxN=0 After truncation, sequences with more than 0 Ns will be discarded. (DADA2 requires sequences contain no Ns)
#truncQ = 2 Truncate reads at the first instance of a quality score less than or equal to 2
#rm.phix = TRUE Discard reads that match against the phiX genome
#maxEE=c(2, 2) After truncation, reads with higher than 2 "expected errors" will be discarded
#minLen = 60 Remove reads with length less than 60 (note these should have already been removed by cutadapt)
#multithread = TRUE input files are filtered in parallel

out2 <- filterAndTrim(cutFs2, filtFs2, cutRs2, filtRs2, maxN = 0, maxEE = c(2, 2), 
                      truncQ = 2, minLen = 60, rm.phix = TRUE, compress = TRUE, 
                      multithread = TRUE)
out2


# Identification of ASVs & Generate an error model 2 ------------------------

#Model the error rates using both the forward and reverse reads 
errF2 <- learnErrors(filtFs2, multithread = TRUE)
errR2 <- learnErrors(filtRs2, multithread = TRUE)

#plot errors
plotErrors(errF2, nominalQ = TRUE)

#The error rates for each possible transition (e.g. A→C, A→G) are shown
#Red line - expected based on the quality score. (These are plotted when nominalQ = TRUE is included in the plot command)
#Black line - estimate
#Black dots - observed

#Dereplication 
exists2 <- file.exists(filtFs2) # check that all the samples 
#are still present after filtering
derepFs2 <- derepFastq(filtFs2[exists2], verbose=TRUE)
derepRs2 <- derepFastq(filtRs2[exists2], verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs2) <- sample.names[exists2]
names(derepRs2) <- sample.names[exists2]

#Inferrence of ASVs
dadaFs2 <- dada(derepFs2, err = errF2, multithread = TRUE)
dadaRs2 <- dada(derepRs2, err = errR2, multithread = TRUE)
#Important that you want to run the error and inference steps on datasets generated from a single Illumina run,
#different runs can have different error structures


# Merging paired end reads 2 ------------------------------------------------

mergers2 <- mergePairs(dadaFs2, derepFs2, dadaRs2, derepRs2, verbose=TRUE)

#Making our ASV matrix 
seqtab2 <- makeSequenceTable(mergers2)
dim(seqtab2)

#Chimera detection and removal
seqtab2.nochim <- removeBimeraDenovo(seqtab2, method="consensus", 
                                     multithread=TRUE, verbose=TRUE)

dim(seqtab2.nochim)

sum(seqtab2.nochim)/sum(seqtab2)

#range of ASV lengths
table(nchar(getSequences(seqtab2.nochim)))

#Sequence tracking sanity check 
getN <- function(x) sum(getUniques(x)) 

track2 <- cbind(out2, sapply(dadaFs2, getN), 
                sapply(dadaRs2, getN), 
                sapply(mergers2, getN), 
                rowSums(seqtab2.nochim))
colnames(track2) <- c("input", "filtered", 
                      "denoisedF", "denoisedR", 
                      "merged", "nonchim")
rownames(track2) <- sample.names
track2


# Assigning taxonomy 2 ------------------------------------------------------
#in this example it is using 12S - the locqtion of reference database will change
taxa12S <- assignTaxonomy(seqtab2.nochim, 
                          "/pub14/tea/nsc206/NEOF/metabarcoding_workshop/taxonomy/12S-V5_Reference_Database.fasta",
                          multithread=TRUE, verbose = T)

taxa12S.print <- taxa12S
rownames(taxa12S.print) <- NULL
head(taxa12S.print)


# Further analysis 2  -------------------------------------------------------

#Load table which contains meta data
meta<-read.csv(".csv", 
               row.names = 1)

#How to remove samples from data i.e. too little reads, negative controls etc
seqtab2.rmcontrol<-seqtab2.nochim[-31,] #31 are sample names - these will change with your data
meta2.rmcontrol<-meta[-31,]

#Remove sameples with failed/low sequence numbers - all examples of sample IDs and read depth
#Check number of sequences per sample
rowSums(seqtab2.rmcontrol)
# S16, S25, S40 and S43 all have low sequence counts. So we will remove these rows.
seqtab2.rmlow<-seqtab2.rmcontrol[-c(1,4,13,16),]
meta2.rmlow<-meta2.rmcontrol[-c(1,4,13,16),]

# Print the minimum sequence number in one sample.
min(rowSums(seqtab2.rmlow))
# The lowest number of sequences in one sample is now 26841.


# Rarefaction curves 2 ------------------------------------------------------

rarecurve(seqtab2.rmlow, step=100, col=meta2.rmlow$COLOUR, lwd=2, ylab="ASVs", label=F)
# add a vertical line to represent the fewest sequences in any sample
abline(v=(min(rowSums(seqtab2.rmlow))))


# Alpha diversity 2 ---------------------------------------------------------

seqtab2.rmlow2<-t(as.data.frame(seqtab2.rmlow))
phylo_asv2 <- otu_table(seqtab2.rmlow2, taxa_are_rows=TRUE)
phylo_tax2 <- tax_table(taxa12S)
phylo_samples2 <- sample_data(meta2.rmlow)

phylo2 <- phyloseq(phylo_asv2, phylo_tax2, phylo_samples2)

sample_names(phylo2)
rank_names(phylo2)
sample_variables(phylo2)

#Two alpha diversity metrics to plot - Shannon's diversity index and Simpson diversity index 
plot_richness(phylo2, measures=c("Shannon", "Simpson"), 
              color = "SITE")
plot_richness(phylo2, x="SITE", measures=c("Shannon", "Simpson"), 
              color = "SITE") + geom_boxplot()


# Beta diversity 2 ----------------------------------------------------------

ps.prop2 <- transform_sample_counts(phylo2, function(otu) otu/sum(otu))

ord.nmds.bray2 <- ordinate(ps.prop2, method="NMDS", distance="bray")
plot_ordination(ps.prop2, ord.nmds.bray2, 
                color="SITE", title="Bray NMDS 12S-V5")

#Plot the proportion of ASV sequences within each sample that belong to different taxonomic families
plot_bar(ps.prop2, fill = "Family")+
  geom_bar(aes(color=Family, fill=Family), 
           stat="identity", position="stack")+
  facet_grid(~SITE, scales = "free", space = "free")

#Remove the negative control and re-plot the NMDS
#To remove specific samples from a phyloseq object
SamplesToRemove <- c("S62", "S63", "S64", "S65", "S66", "S67")
phylo2.nofilt <- subset_samples(phylo2, !(sample_names(phylo2) %in% SamplesToRemove))

#Transform the data into proportions and re-plot NMDS
ps.prop.nofilt <- transform_sample_counts(phylo2.nofilt, function(otu) otu/sum(otu))
ord.nmds.bray.nofilt <- ordinate(ps.prop.nofilt, method="NMDS", distance="bray")
plot_ordination(ps.prop.nofilt, ord.nmds.bray.nofilt, color="SITE", 
                title="Bray NMDS 12S-V5 lake samples")

#Bray-Curtis distances using the distance function and perform a PERMANOVA using adonis function
#from Vegan to check whether there is any significant different separation of samples by site 
bray.dist.nofilt<-distance(ps.prop.nofilt, method="bray")
sampledf.nofilt <- data.frame(sample_data(phylo2.nofilt))
adonis2(bray.dist.nofilt ~ SITE, data = sampledf.nofilt)


#the differences between both scripts is just adding a nummber 2 aat each section to split it into a different analysis. If you are looking at greater than two primer pairs you can adjust these numbers
#i.e. FWD3 - REV3 & cutFs3 - cutRs3 ... FWD4 .... FWD5 ... cutFs6 etc 



