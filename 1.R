# Task1

install.packages("seqinr")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")

#library(package)
library(Biostrings)

# load data fishes.fna.gz
# with Biostring

# načtení cesty
setwd('V:/MPA_PRG/exercise_02')


library(Biostrings)
library(seqinr)

#pomocí biostrings
dna <- readDNAStringSet('fishes.fna.gz', format = 'fasta')

# pomocí seqinr
seq <- read.fasta("fishes.fna.gz")

# task 2......................................................

length(dna) # delka
width(dna[1]) # velikost
names(dna) # jmeno

seq1 <- dna[1]
seq1

seq1_sequence <- dna[1] # sequence
seq1_sequence

seq1_string <- toString(dna[1]) # dna na string (retezec)
seq1_string

help("XStringSet") # help

# task 3 ............
BiocManager::install("pwalign")
seq1 <- dna[1]
seq2 <- dna[2]
data("BLOSUM62")
globalAlign <- pairwiseAlignment(seq1, seq2, substitutionMatrix = "BLOSUM62", gapOpening = -1, gapExtension = -1)
globalAlign


# task 4 .........
name_list <- c("anna", "jana", "kamil", "norbert", "pavel", "petr", "stanislav", "zuzana")
grep("jana", name_list, perl=TRUE) # dá ti popřadí kde je to v seznamu ta jana
grep("n+", name_list, perl = TRUE) # od n a zatím n kolikrát chceme od 1 az nekonecno
grep("n{2}", name_list, perl = TRUE) # dve nn za sebou 
grep("^n", name_list, perl = TRUE)# zacina na n
grep("Anna|Jana", name_list, perl = TRUE) # ???
grep("^z.*a$", name_list, perl = TRUE) # zacina z za tim jenou cokoliv pak kolik chceme pak a koncito a 

# task 5...................
sek <- readDNAStringSet('fishes.fna.gz', format = 'fasta')
mid_f <- "AGGCT"
mid_b <- "CAATG"

# mělo by se najít něco okolo 4351
# ^AGGCT.CAATG$
# mid_f = ACGAGTGCGT
# mid_b = ACGCACTCGT
sek
pattern <- grep("^ACGAGTGCGT.*ACGCACTCGT$", sek, perl=TRUE)
pattern
sek[98] # asi sparavne
# together <- paste0("^",midf,".*",midb,"$")


# task 6.................
MID <- read.csv("fishes_MIDs.csv", sep = ";")
path <- setwd("V:/MPA_PRG/exercise_02")
  
  
Demultiplexer <- function(path, forward_mids, reverse_mids, labels_list){
  
  MID <- read.csv("path/fishes_MIDs.csv", sep = ";")
  return MID
}
  

