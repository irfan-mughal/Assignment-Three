library("seqinr")
library("R.utils")
library("rBLAST")
library("ape")
library("ORFik")
library("Biostrings")
source("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/mutblast_functions.R")


#Question 1 - Download whole set of E.coli gene DNA sequences

download.file("ftp://ftp.ensemblgenomes.org/pub/bacteria/release-42/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/cds/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz",
              destfile = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz")

# decompressing the file
R.utils::gunzip("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz", overwrite = TRUE)
# Creating a blast database
makeblastdb("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa", dbtype = "nucl", "-parse_seqids")


#Question 2 - Download sample fasta sequences
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/sample.fa",
              destfile = "sample.fa")
# read into r script
ec <- read.fasta("sample.fa")
mySeq <- ec [[69]]
str(mySeq)
#calculating the GC 
seqinr::GC(mySeq)
#length of myseq in basepairs
seqinr::getLength(mySeq)


#Question 3 - Create Blast function with provided R function and do blast searches
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/mutblast_functions.R",
              destfile = "mutblast_functions.R")


# Use Blast to identify ecoli sequense
bmut <- myblastn_tab (myseq = mySeq, db ="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
head(bmut)
str(bmut)
mysseqid <- as.character(bmut$sseqid)
mysseqid
# top 3 hit genes with percent identity, E-values and bit scores
hits <- as.character(bmut$sseqid[1:3])
hits

#Question 4 - Find number of mismatches  between original and mutated seq?
mutator(mySeq, 50)
Myseq_mut <- mutator(mySeq, 50)
myblastn_tab( Myseq_mut, db = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")


#Question 5 - Perform mutating and blasting sequence to determine number of proprtion of sites need to be altered
Result <- function(myseq = mySeq, nmut = nmut){
  seqmut <- mutator(mySeq, nmut)
  res <- myblastn_tab(myseq = seqmut, db = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
  if (is.null(res) ) {
    myres = 0} else {myres=1
    }
  return(myres)
}
replicate(n = 100, expr = Result(myseq = mySeq, nmut =100))
mean(replicate(n=100, expr = Result(mySeq, 100)))
str(mean)
finalres <- function(n){
  mean(replicate(n, expr = Result(mySeq, 100)))
}
n <- c(10, 100, 200, 300, 400)
sapply(n, finalres)

#Ques 6 -  Plot a chart with increasing proprotion of mutated bases

plot(n, sapply(n,finalres), xlab = "number of sites", ylab = "proportion of mutated bases", main = "Increased proportion of mutated bases", type = 'o')

