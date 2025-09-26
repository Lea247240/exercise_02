# Task - instalovat seqinr, BiocManager, Biostrings 

#install.packages("seqinr") 
library(seqinr)

#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

# BiocManager::install("Biostrings")

library(Biostrings)

#....................................................................................
# Task 1 
#Load the DNA sequence fishes.fna.gz using functions from the seqinr package 
#and the Biostrings package. Note the differences between the created variables.


# find path
#setwd('V:/MPA_PRG/exercise_02')# PC VUT
setwd('D:/VUT/4-5rocnik/moje/MPA-PRG/exercise_02')

# load data fishes.fna.gz
# with Biostring
dna <- readDNAStringSet('fishes.fna.gz', format = 'fasta')

# pomocí seqinr
seq <- read.fasta(gzfile("fishes.fna.gz"), seqtype="DNA")

# Task 2.............................................................................
# Next, focus on the Biostrings package. Practice working with loaded data

length(dna) # pocet vsech sekvenci - Check the number of loaded sequences
width(dna[1]) # delka aktualni sekvece [bp] - Determine the lengths of each sequence
names(dna) # jmeno sekvence (jeji hlavicka) - View the sequence names (FASTA headers)

seq1 <- dna[1] # prvni sekvence - Assign the first sequence including the name to the variable seq1
seq1

seq1_sequence <- dna[1] # prirazeni 1 sequence - Assign the first sequence without the name to the variable seq1_sequence
seq1_sequence

seq1_string <- toString(dna[1]) # dna na string (retezec) - Assign the first sequence as a vector of characters to the variable seq1_string
seq1_string

help("XStringSet") # help - Learn more about the XStringSet class and the Biostrings package


# Task 3 ..............................................................................
#Globally align the two selected sequences using the BLOSUM62 matrix, a gap opening
#cost of -1 and a gap extension cost of 1.

BiocManager::install("pwalign")
library(pwalign)

seq1 <- dna[1]
seq2 <- dna[2]

data("BLOSUM62")

globalAlign <- pairwiseAlignment(seq1, seq2, substitutionMatrix = "BLOSUM62", gapOpening = -1, gapExtension = -1, type = "global") # zarovnani s Blosum 62
globalAlign


# Task 4 ...............................................................................
#Practice working with regular expressions

name_list <- c("anna", "jana", "kamil", "norbert", "pavel", "petr", "stanislav", "zuzana") # list jmen

#Search for name "jana":
grep("jana", name_list, perl=TRUE) # dá ti indexy kde je to v seznamu ta jana = 2

#Search for all names containing letter "n" at least once:
grep("n+", name_list, perl = TRUE) # indexy, kde se vyskytuje alespoň jedno "n" za sebou (nn, nnn, …)

#Search for all names containing letters "nn":
grep("n{2}", name_list, perl = TRUE) # dve nn za sebou 

#Search for all names starting with "n":
grep("^n", name_list, perl = TRUE)# zacina na n

#Search for names "Anna" or "Jana":
grep("Anna|Jana", name_list, perl = TRUE) # indexy, kde je "Anna" NEBO "Jana" (case-sensitive → velká A) - nic v seznamu jsou mala pismena

#Search for names starting with "z" and ending with "a":
grep("^z.*a$", name_list, perl = TRUE) # zacina z za tim cokoliv pak kolikrat chceme a koncito a 


# Task 5..................................................................
#Load an amplicon sequencing run from 454 Junior machine fishes.fna.gz.
#Get a sequence of a sample (avoid conditional statements), that is tagged by forward
#and reverse MID ACGAGTGCGT.
#How many sequences are there in the sample?

sek <- readDNAStringSet('fishes.fna.gz', format = 'fasta')

# mělo by se najít něco okolo 4351
# ^AGGCT.CAATG$
mid_f = "ACGAGTGCGT" 
mid_b = "ACGCACTCGT"

sek
pattern <- grep("^ACGAGTGCGT.*ACGCACTCGT$", sek, perl=TRUE)# nalezeni vzoru
pattern
sek[98] # sparavne

together <- paste0("^",mid_f,".*",mid_b,"$") # to je vlastne toto:"^ACGAGTGCGT.*ACGCACTCGT$"
pattern1 <- grep(together, sek, perl=TRUE)# nalezeni vzoru
pattern
sek[98] # sparavne

# Task 6......................................................................
#Create a function Demultiplexer() for demultiplexing of sequencing data.
#Input:
  #a string with path to fasta file
  #a list of forward MIDs
  #a list of reverse MIDs
  #a list of sample labels
#Output:
  # - fasta files that are named after the samples and contain sequences of the
  #sample without MIDs (perform MID trimming)
  # - table named report.txt containing the labels of the samples and the
  #number of sequences each sample has
  # - Check the functionality again on the fishes.fna.gz file, the list of samples and
  #MIDs can be found in the corresponding table fishes_MIDs.csv.

#..........
#path <- setwd("V:/MPA_PRG/exercise_02") #- skolni pc
path <- setwd('D:/VUT/4-5rocnik/moje/MPA-PRG/exercise_02') #- doma pc

# nacteni souboru csv
library(Biostrings)
sek <- readDNAStringSet('fishes.fna.gz', format = 'fasta')
MID <- read.csv("fishes_MIDs.csv", sep = ";")

# převedení sloupcu z tabulky na seznam 
f_MIDs <- as.character(MID$FBarcodeSequence)
r_MIDs <- as.character(MID$RBarcodeSequence)
labels_list <- as.character(MID$SampleID)

# .................................................... zjistit jestli tam vubec něco takové je 
#grep("^ACGAGTGCGT", sek) # f_MIDs[1] je ve vícero
#grep("ACGAGTGCGT$", sek) # r_MIDs[1] je v 91460
#..............................................................................................


# vytvoreni listu kam se ulozi orizle sekvence (dlouhy jako SampleID)
save_trimmed <- vector("list", length(labels_list))
#pridame nazvy SampleID ke kterym prodame nase orizle sekvence
names(save_trimmed) <- as.character(labels_list) 


# forcykl
for (i in 1:length(f_MIDs)){
  fwd <- f_MIDs[i]
  rev <- r_MIDs[i]
  label <- labels_list[i]
  together <- paste0("^",fwd,".*",rev,"$")
  match <- grep(together, sek, perl=TRUE)
  print(match)
  if (length(match) > 0){
    #print(label)
    #print(sek[match])
    
    terminned <- substr(sek[match], nchar(fwd) + 1, nchar(sek[match]) - nchar(rev)) 
    # nchar(fwd) - delka f_mid
    # nchar(rev) - delka r_mid
    # nchar(sek[match]) - delka aktualni sekvence
    print(terminned)
    
    save_trimmed[[label]] <- c(save_trimmed[[label]], terminned) # prirazeni orizle sekvence k nasvu v listu save_trimmed
    
    }
}

# vytvoření FASTA souborů pro každý vzorek
for (label in labels_list) {
  label_str <- as.character(label)
  seqs <- save_trimmed[[label_str]]
  if (length(seqs) > 0) {
    writeXStringSet(DNAStringSet(seqs), paste0(label_str, ".fasta"))
  }
}

# vytvoření report.txt
report <- data.frame(
  SampleID = labels_list,
  Number_of_sequences = sapply(save_trimmed, length)
)
write.table(report, "report.txt", sep = "\t", row.names = FALSE, quote = FALSE)

cat("Demultiplexing completed. FASTA files and report.txt created.\n")







  

Demultiplexer <- function(path, forward_mids, reverse_mids, labels_list){
  setwd(path)
  #together <- paste0(path, "/fishes_MIDs.csv")
  MID2 <- read.csv("fishes_MIDs.csv", sep = ";")
  sek <- readDNAStringSet('fishes.fna.gz', format = 'fasta')
  
  return (MID2)
}

Demultiplexer(path = 'D:/VUT/4-5rocnik/moje/MPA-PRG/exercise_02',forward_mids = forward_mids, reverse_mids = reverse_mids, labels_list = labels_list ) 

#....................................................................................

grep("^ACGAGTGCGT.*ACGCACTCGT$", sek)
rev <- as.character(MID$RBarcodeSequence)
rev[1]

seq <- DNAStringSet(rev[1])
seq
revers_comp <- reverseComplement(seq)
revers_comp # pak zpátky na string a dát to místo toho revs
