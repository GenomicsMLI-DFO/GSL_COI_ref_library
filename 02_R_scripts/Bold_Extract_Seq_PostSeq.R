
# Info --------------------------------------------------------------------

# Code to improve the first version of NWA reference sequence db
# Post-seq list of taxa came from a comparison of real DNAdata
# 
# Audrey Bourret
# 2021-01-14
#


# Library -----------------------------------------------------------------

if(!require(readxl)){install.packages("readxl")}
library(readxl)

if(!require(tidyverse)){install.packages("tidyverse")}
library(tidyverse)

#if(!require(rentrez)){install.packages("rentrez")}
#library(rentrez)

if(!require(bold)){install.packages("bold")}
library(bold)

library(Hmisc) # for %nin%

library(Biostrings)
library(msa)
library(ape)

source("02_R_scripts/00_Functions/fct_BOLD_v2.R")
source("02_R_scripts/00_Functions/fct_DNAseq.R")

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Data --------------------------------------------------------------------

list.files("01_Raw_data/01_SpLists/")

postseqSP.join <- read.csv("01_Raw_data/01_SpLists/PostSeqSP_revisedByTaxo_validateWORMS_final.csv", 
                          fileEncoding = "ISO-8859-1")

postseq.addSP.join <- read.csv("01_Raw_data/01_SpLists/PostSeq.AddSP_revisedByTaxo_validateWORMS_final.csv", 
                          fileEncoding = "ISO-8859-1")

postseq.genus.join <- read_excel("01_Raw_data/01_SpLists/PostSeqSP_revisedByTaxo.xlsx",
                         sheet = "BLAST-GENUS")


postseq.genus.join <- postseq.genus.join %>% left_join(postseqSP.join %>% select(genus, family) %>% unique(),
                       by = c("Genus" = "genus") )

head(postseqSP.join)
#head(addSP.join)
head(postseq.genus.join)

postseq.list <- postseqSP.join %>% filter(WithinGSL %in% c(1:2)) %>% 
                             pull(Species) %>% as.character()
# Check for duplicates
postseq.list[postseq.list %>% duplicated]

# Small stats - overall before and after taxo revision
prioriSP.join %>% mutate(GSL = ifelse(WithinGSL %in% c(1:2), "YES", "NO")) %>% 
                  group_by(RequestedBy, Categorie, GSL) %>% 
                  summarise(N.row.check = n(),
                            N.species = length(unique(Species)),
                            N.family = length(unique(family)),
                            N.order = length(unique(order)))

# Small stats - taxo revision
#prioriSP.join %>% mutate(GSL = ifelse(WithinGSL %in% c(1:2), "YES", "NO"),
#                         PROB = ifelse(TaxonomicProblems == 1 | MorphoIDProblems == 1, 1, 0)) %>%
#  filter(GSL == "YES") %>% 
#  group_by(Categorie, PROB) %>% 
#  summarise(N.row.check = n(),
#            N.species = length(unique(Species)),
#            N.family = length(unique(family)),
#            N.order = length(unique(order)))


prioriSP.join %>% filter(WithinGSL %in% c(1:2)) %>% pull(Species) 
  
#add.list <- addSP.join %>% filter(WithinGSL %in% c(1:2),
#                                  Step == "Addition") %>% 
#                                  pull(Species) %>% as.character()

postseq.add.postgenus.list <- postseq.addSP.join %>% filter(WithinGSL %in% c(1:2),
                                  Step == "BLAST-GENUS") %>% 
                                pull(Species) %>% as.character()  


postseq.BlackList <- read_excel( "01_Raw_data/01_SpLists/PostseqSP_revisedByTaxo.xlsx",
                         sheet = "BLAST-BlackList")

cat.list <- bind_rows(postseqSP.join %>% select(Name = Species, Categorie, family),
                      postseq.addSP.join %>% select(Name = Species, Categorie, family),
                      postseq.genus.join %>% select(Name = Genus, Categorie, family) )
cat.list %>% head()

# Saved results -----------------------------------------------------------

# load(file.path("01_Raw_data/02_Sequences/02_Post-sequencing/postseq_raw_2021-01-14.Rdata"))
# load(file.path("01_Raw_data/02_Sequences/02_Post-sequencing/postseq_bin_2021-01-14.Rdata"))
# 
# load(file.path("01_Raw_data/02_Sequences/02_Post-sequencing/postseq.wrongtaxo_raw_2021-01-16.Rdata"))
# load(file.path("01_Raw_data/02_Sequences/02_Post-sequencing/postseq.wrongtaxo_bin_2021-01-19.Rdata"))
# 
# load(file.path("01_Raw_data/02_Sequences/02_Post-sequencing/postseq.genus_raw_2021-01-19.Rdata"))
# load(file.path("01_Raw_data/02_Sequences/02_Post-sequencing/postseq.genus_bin_2021-01-19.Rdata"))
# 
# load(file.path("01_Raw_data/02_Sequences/02_Post-sequencing/postseq.add.postgenus_raw_2020-01-20.Rdata"))
# load(file.path("01_Raw_data/02_Sequences/02_Post-sequencing/postseq.add.postgenus_bin_2020-01-20.Rdata"))
# 
# Meta.seq.total <- bind_rows(postseq.seq$Meta.seq,
#                             postseq.wrongtaxo.seq$Meta.seq,
#                             #postseq.add.seq$Meta.seq,
#                             postseq.genus.seq$Meta.seq,
#                             postseq.add.postgenus.seq$Meta.seq)
# 
# # Stats on what was loaded
# names(priori.res)
# 
# nrow(postseq.res$Meta.filt) 
# nrow(postseq.wrongtaxo.res$Meta.filt) 
# nrow(postseq.add.postgenus.res$Meta.filt)
# nrow(postseq.genus.res$Meta.filt)
# 
# 

# First round of BOLD search ----------------------------------------------


postseq.res <- BOLD.sratch(postseq.list, taxo.level = "Species", min.l.pb = 650)

postseq.seq <- BOLD.seq(postseq.res$Meta.filt, min.l.pb = 650)

names(postseq.res)
names(postseq.seq)

postseq.seq$DNA %>% names()

#save(list = c("postseq.res", "postseq.seq"), file = file.path("01_Raw_data/postseq_raw_2021-01-14.Rdata"))


# Stats

# Number of species that were looked
postseq.res["Summary"]$Summary %>% nrow()  # Chlamys was already look

# Number of species retained, after filtration
postseq.res["Summary"]$Summary %>% filter(n.BOLD.filt != 0) %>% nrow()
postseq.res["Summary"]$Summary %>% filter(n.BOLD.filt.NWA != 0) %>% nrow()

# N BIN
postseq.res["Summary"]$Summary %>% group_by(n.BIN.raw) %>% summarise(N = n())
postseq.res["Summary"]$Summary %>% group_by(n.BIN.raw.NWA) %>% summarise(N = n())


postseq.seq["Summary"]$Summary %>% filter(n.BIN.final != 0) %>% nrow()
postseq.seq["Summary"]$Summary %>% filter(n.BIN.final.NWA != 0) %>% nrow()


priori.seq["Summary"]$Summary %>% group_by(n.BIN.final) %>% summarise(N = n())
priori.seq["Summary"]$Summary %>% group_by(n.BIN.final.NWA) %>% summarise(N = n())


postseq.seq["Summary"]$Summary %>% filter(n.DNA.final != 0) %>% nrow()


priori.seq["Summary"]$Summary %>% pull(Alig.length) %>% table()


# First validation: BIN ---------------------------------------------------

postseq.bin <- postseq.seq$Meta.seq %>% select(Species.search, bin_uri) %>% 
  unique() %>%
  mutate(BIN.species = NA,
         BIN.species.N = NA,
         BIN.species.freq = NA) 


for(x in 1:nrow(postseq.bin)){
  #for(x in 1:10){
  
  print(postseq.bin$Species.search[x])
  
  BIN.int <- bold_seqspec(bin = postseq.bin$bin_uri[x])
  
  BIN.sp <- BIN.int %>% group_by(species_name) %>% 
    summarise(N = n()) %>% 
    mutate(freq = N / sum(N)) %>% 
    filter(species_name != "") %>% 
    arrange(desc(N))
  
  postseq.bin$BIN.species[x]      <- BIN.sp$species_name[1]
  postseq.bin$BIN.species.N[x]      <- BIN.sp$N[1]
  postseq.bin$BIN.species.freq[x] <- BIN.sp$freq[1]  
  
} 

#save(list = c("postseq.bin"), file = file.path("01_Raw_data/02_Sequences/02_Post-sequencing/postseq_bin_2021-01-14.Rdata"))
     
postseq.bin

postseq.bin %>% View()
postseq.bin %>% names()

postseq.bin %>% mutate(check = ifelse(Species.search == BIN.species, 0,1)) %>% summarise(sum = sum(check))
postseq.bin %>% mutate(check = ifelse(Species.search == BIN.species, 0,1)) %>% 
  left_join(postseqSP.join %>% select(Species, Categorie, TaxonomicProblems, MorphoIDProblems, Comments),
            by = c("Species.search" = "Species")) %>% 
  write_csv2("01_Raw_data/02_Sequences/02_Post-sequencing/postseq_bin_summary_2021-01-14.csv")
  #View()

# Second round of BOLD search ----------------------------------------------

# Is there some species with no info, but maybe with older name?

postseq.seq.list <- postseq.seq$Meta.seq %>% pull(Species.search) %>% unique()

postseq.wrongtaxo.list <- postseqSP.join %>% filter(!is.na(InvalidSpName),
                                Species %nin% postseq.seq.list) %>% select(RigthSpecies = Species, InvalidSpName)
postseq.wrongtaxo.list

postseq.wrongtaxo.res <- BOLD.sratch(postseq.wrongtaxo.list$InvalidSpName, taxo.level = "Species", min.l.pb = 650)
names(postseq.wrongtaxo.res)

postseq.wrongtaxo.res$Summary
postseq.wrongtaxo.res$Meta.filt <- postseq.wrongtaxo.res$Meta.filt %>% left_join(postseq.wrongtaxo.list, by = c("Species.search" = "InvalidSpName")) 

postseq.wrongtaxo.seq <- BOLD.seq(postseq.wrongtaxo.res$Meta.filt, min.l.pb = 650)

#save(list = c("postseq.wrongtaxo.res", "postseq.wrongtaxo.seq"), file = file.path("01_Raw_data/02_Sequences/02_Post-sequencing/postseq.wrongtaxo_raw_2021-01-16.Rdata"))


postseq.wrongtaxo.bin <- postseq.wrongtaxo.seq$Meta.seq %>% select(Species.search, RigthSpecies, bin_uri) %>% 
  unique() %>%
  mutate(BIN.species = NA,
         BIN.species.N = NA,
         BIN.species.freq = NA) 


for(x in 1:nrow(postseq.wrongtaxo.bin)){
  #for(x in 1:10){
  
  print(postseq.wrongtaxo.bin$Species.search[x])
  
  BIN.int <- bold_seqspec(bin = postseq.wrongtaxo.bin$bin_uri[x])
  
  BIN.sp <- BIN.int %>% group_by(species_name) %>% 
    summarise(N = n()) %>% 
    mutate(freq = N / sum(N)) %>% 
    filter(species_name != "") %>% 
    arrange(desc(N))
  
  postseq.wrongtaxo.bin$BIN.species[x]      <- BIN.sp$species_name[1]
  postseq.wrongtaxo.bin$BIN.species.N[x]      <- BIN.sp$N[1]
  postseq.wrongtaxo.bin$BIN.species.freq[x] <- BIN.sp$freq[1]  
  
} 

#save(list = c("postseq.wrongtaxo.bin"), file = file.path("01_Raw_data/02_Sequences/02_Post-sequencing/postseq.wrongtaxo_bin_2021-01-19.Rdata"))



postseq.wrongtaxo.bin

postseq.wrongtaxo.bin %>% View()
postseq.wrongtaxo.bin %>% names()

postseq.wrongtaxo.bin %>% mutate(check = ifelse(Species.search == BIN.species, 0,1)) %>% summarise(sum = sum(check))
postseq.wrongtaxo.bin %>% mutate(check = ifelse(Species.search == BIN.species, 0,1)) %>% 
  left_join(postseqSP.join %>% select(Species, Categorie, TaxonomicProblems, MorphoIDProblems, Comments),
            by = c("RigthSpecies" = "Species")) %>% 
  write_csv2("01_Raw_data/02_Sequences/02_Post-sequencing/postseq.wrongtaxo_bin_summary_2020-01-19.csv")
#View()


# Last round of BOLD -----------------------------------------------------

# genus, NWA only, with BIN not previously observed

postseq.genus.res <- BOLD.sratch(postseq.genus.join$Genus, taxo.level = "Genus", min.l.pb = 650)

names(postseq.genus.res)

postseq.genus.res$Summary


#load(file.path("01_Raw_data/02_Sequences/01_Pre-sequencing/priori_bin_2020-09-30.Rdata"))
#load(file.path("01_Raw_data/02_Sequences/01_Pre-sequencing/wrongtaxo_bin_2020-09-30.Rdata"))
#load(file.path("01_Raw_data/02_Sequences/01_Pre-sequencing/add_bin_2020-09-30.Rdata"))
#load(file.path("01_Raw_data/02_Sequences/01_Pre-sequencing/genus_bin_2020-09-30.Rdata"))
#load(file.path("01_Raw_data/02_Sequences/01_Pre-sequencing/add.postgenus_bin_2020-09-30.Rdata"))

ls(pattern = "bin")

# What we already have
SP.BIN <- c(postseq.bin$bin_uri,postseq.wrongtaxo.bin$bin_uri,
            priori.bin$bin_uri, add.bin$bin_uri, wrongtaxo.bin$bin_uri, genus.bin$bin_uri, add.postgenus.bin$bin_uri 
            ) %>% unique()
SP.BIN

postseq.genus.res$Meta.filt %>% nrow()

postseq.genus.res$Meta.filt.NWA.NewBIN <- postseq.genus.res$Meta.filt %>% filter(bin_uri %nin% SP.BIN,
                                       NWA_code == 1) 
postseq.genus.res$Meta.filt.NWA.NewBIN %>% nrow()                                      

postseq.genus.res$Meta.filt.NWA.NewBIN %>% select(bin_uri, species_name, Species.search) %>% 
  unique() %>% View()
 # write_csv2("clipboard")



BIN.genus.filt <- postseq.genus.join %>% filter(!is.na(BIN)) %>% pull(BIN) %>% paste(collapse = ";") %>% str_split(";") %>% `[[`(1)
BIN.genus.filt

postseq.genus.res$Meta.filt.NWA.NewBIN.final <- postseq.genus.res$Meta.filt.NWA.NewBIN %>% filter(bin_uri %in% BIN.genus.filt) 

postseq.genus.seq <- BOLD.seq(postseq.genus.res$Meta.filt.NWA.NewBIN.final, min.l.pb = 650)

#save(list = c("postseq.genus.res", "postseq.genus.seq"), file = file.path("01_Raw_data/02_Sequences/02_Post-sequencing/postseq.genus_raw_2021-01-19.Rdata"))


postseq.genus.bin <- postseq.genus.seq$Meta.seq %>% select(Species.search, bin_uri) %>% 
  unique() %>%
  mutate(BIN.species = NA,
         BIN.species.N = NA,
         BIN.species.freq = NA) 


for(x in 1:nrow(postseq.genus.bin)){
  #for(x in 1:10){
  
  print(postseq.genus.bin$Species.search[x])
  
  BIN.int <- bold_seqspec(bin = postseq.genus.bin$bin_uri[x])
  
  BIN.sp <- BIN.int %>% group_by(species_name) %>% 
    summarise(N = n()) %>% 
    mutate(freq = N / sum(N)) %>% 
    filter(species_name != "") %>% 
    arrange(desc(N))
  
  postseq.genus.bin$BIN.species[x]      <- BIN.sp$species_name[1]
  postseq.genus.bin$BIN.species.N[x]      <- BIN.sp$N[1]
  postseq.genus.bin$BIN.species.freq[x] <- BIN.sp$freq[1]  
  
} 

#save(list = c("postseq.genus.bin"), file = file.path("01_Raw_data/02_Sequences/02_Post-sequencing/postseq.genus_bin_2021-01-19.Rdata"))

postseq.genus.bin

postseq.genus.bin %>% View()
postseq.genus.bin %>% names()

postseq.genus.bin %>% mutate(check = ifelse(Species.search == BIN.species, 0,1)) %>% summarise(sum = sum(check))
postseq.genus.bin %>% mutate(check = ifelse(Species.search == BIN.species, 0,1)) %>% 
  left_join(postseqSP.join %>% select(Species, Categorie, TaxonomicProblems, MorphoIDProblems, Comments),
            by = c("Species.search" = "Species")) %>% 
  write_csv2("01_Raw_data/02_Sequences/02_Post-sequencing/postseq.genus_bin_summary_2021-01-19.csv")
#View()


# Check also for the addition

postseq.add.postgenus.res <- BOLD.sratch(c(postseq.add.postgenus.list), taxo.level = "Species", min.l.pb = 650)

postseq.add.postgenus.seq <- BOLD.seq(postseq.add.postgenus.res$Meta.filt, min.l.pb = 650)

names(postseq.add.postgenus.res)
names(postseq.add.postgenus.seq)

#save(list = c("postseq.add.postgenus.res", "postseq.add.postgenus.seq"), file = file.path("01_Raw_data/02_Sequences/02_Post-sequencing/postseq.add.postgenus_raw_2020-01-20.Rdata"))

postseq.add.postgenus.bin <- postseq.add.postgenus.seq$Meta.seq %>% select(Species.search, bin_uri) %>% 
  unique() %>%
  mutate(BIN.species = NA,
         BIN.species.N = NA,
         BIN.species.freq = NA) 


for(x in 1:nrow(postseq.add.postgenus.bin)){
  #for(x in 1:10){
  
  print(postseq.add.postgenus.bin$Species.search[x])
  
  BIN.int <- bold_seqspec(bin = postseq.add.postgenus.bin$bin_uri[x])
  
  BIN.sp <- BIN.int %>% group_by(species_name) %>% 
    summarise(N = n()) %>% 
    mutate(freq = N / sum(N)) %>% 
    filter(species_name != "") %>% 
    arrange(desc(N))
  
  postseq.add.postgenus.bin$BIN.species[x]      <- BIN.sp$species_name[1]
  postseq.add.postgenus.bin$BIN.species.N[x]      <- BIN.sp$N[1]
  postseq.add.postgenus.bin$BIN.species.freq[x] <- BIN.sp$freq[1]  
  
} 

#save(list = c("postseq.add.postgenus.bin"), file = file.path("01_Raw_data/02_Sequences/02_Post-sequencing/postseq.add.postgenus_bin_2020-01-20.Rdata"))

postseq.add.postgenus.bin

postseq.add.postgenus.bin %>% View()
postseq.add.postgenus.bin %>% names()

postseq.add.postgenus.bin %>% mutate(check = ifelse(Species.search == BIN.species, 0,1)) %>% summarise(sum = sum(check))
postseq.add.postgenus.bin %>% mutate(check = ifelse(Species.search == BIN.species, 0,1)) %>% 
  left_join(postseq.addSP.join %>% select(Species, Categorie, TaxonomicProblems, MorphoIDProblems, Comments),
            by = c("Species.search" = "Species")) %>% 
  write_csv2("01_Raw_data/02_Sequences/02_Post-sequencing/postseq.add.postgenus_bin_summary_2021-01-20.csv")

# Then select also the genus we want to keep: 


# Create_fasta_file -------------------------------------------------------

# Change the name if necessary

change.df <- postseq.wrongtaxo.seq$Meta.seq %>% select(Species.search, RigthSpecies) %>% unique()
change.df
DNA.wrong.newname <- names(postseq.wrongtaxo.seq$DNA) %>% stringi::stri_replace_all_fixed(pattern = change.df$Species.search, replacement = as.character(change.df$RigthSpecies), vectorize_all = F)

names(postseq.wrongtaxo.seq$DNA) <- DNA.wrong.newname 

# One manual change

names(postseq.genus.seq$DNA) <- names(postseq.genus.seq$DNA) %>% str_replace("_Brada_", "_Bradabissa_")


# Join everything 

total.DNA <- c(postseq.seq$DNA, 
               postseq.wrongtaxo.seq$DNA,
               #add.seq$DNA,
               postseq.add.postgenus.seq$DNA,
               postseq.genus.seq$DNA)

# Remove what is black listed

BlackList <- postseq.BlackList %>% mutate(pattern = paste(Species.search, bin_uri, sep = "_"))

paste(BlackList$pattern, collapse = "|")

total.DNA.filt <- total.DNA[names(total.DNA) %>% str_detect(paste(BlackList$pattern, collapse = "|"), negate = T)]   

length(total.DNA)
length(total.DNA.filt)


#writeXStringSet(total.DNA.filt, filepath = "01_Raw_data/02_Sequences/02_Post-sequencing/postseq.DNA.total.2021-01-20.fasta")

total.DNA.filt <- readDNAStringSet("01_Raw_data/02_Sequences/02_Post-sequencing/postseq.DNA.total.2021-01-20.fasta")

DNA.alig.total <- msa(total.DNA.filt, method = "ClustalW")
DNA.alig.total

#writeXStringSet(DNAStringSet(DNA.alig.total), filepath = "01_Raw_data/02_Sequences/02_Post-sequencing/postseq.DNA.align.total.2021-01-20.fasta")

DNA.alig.total <- readDNAStringSet("01_Raw_data/02_Sequences/02_Post-sequencing/postseq.DNA.align.total.2021-01-20.fasta")
DNA.alig.total

# Try to get the start value of the alignment -  which is 120 - 804

SEQ.start <- data.frame(
  id = names(DNA.alig.total),
  A.pos = DNA.alig.total %>% as.data.frame() %>% pull(x) %>%  str_locate(c("A")) %>% as_tibble() %>% pull("start"),
  C.pos = DNA.alig.total %>% as.data.frame() %>% pull(x) %>%  str_locate(c("C")) %>% as_tibble() %>% pull("start"),
  T.pos = DNA.alig.total %>% as.data.frame() %>% pull(x) %>%  str_locate(c("T")) %>% as_tibble() %>% pull("start"),
  G.pos = DNA.alig.total %>% as.data.frame() %>% pull(x) %>%  str_locate(c("G")) %>% as_tibble() %>% pull("start")
)

SEQ.start <- SEQ.start %>% group_by(id) %>% mutate(Start.pos = min(A.pos, C.pos, T.pos, G.pos))

SEQ.stop <- data.frame(
  id = names(DNA.alig.total),
  A.pos = DNA.alig.total %>% as.data.frame() %>% pull(x) %>% stringi::stri_locate_last_coll(c("A")) %>% as_tibble() %>% pull("start"),
  C.pos = DNA.alig.total %>% as.data.frame() %>% pull(x) %>% stringi::stri_locate_last_coll(c("C")) %>% as_tibble() %>% pull("start"),
  T.pos = DNA.alig.total %>% as.data.frame() %>% pull(x) %>% stringi::stri_locate_last_coll(c("T")) %>% as_tibble() %>% pull("start"),
  G.pos = DNA.alig.total %>% as.data.frame() %>% pull(x) %>% stringi::stri_locate_last_coll(c("G")) %>% as_tibble() %>% pull("start")
)

SEQ.stop <- SEQ.stop %>% group_by(id) %>% mutate(Stop.pos = max(A.pos, C.pos, T.pos, G.pos))

# Probably the start and stop position
getmode(SEQ.start$Start.pos)
getmode(SEQ.stop$Stop.pos)


DNA.total <- subseq(DNA.alig.total, start = getmode(SEQ.start$Start.pos), end = getmode(SEQ.stop$Stop.pos))

DNA.total

writeXStringSet(DNA.total, filepath = "01_Raw_data/02_Sequences/02_Post-sequencing/postseq.DNA.align.cut.total.2021-01-20.fasta")


DNA.total.aligOK <- readDNAStringSet(filepath = "01_Raw_data/02_Sequences/02_Post-sequencing/postseq.DNA.align.cut.total.2021-01-20.fasta")

length(DNA.total)
length(DNA.total.aligOK)

# By species - BIN, remove sequences with "N" or "-" 

BIN.SP <- paste(sapply(names(DNA.total.aligOK) %>% str_split("_"), "[",2),
                sapply(names(DNA.total.aligOK) %>% str_split("_"), "[",3),
                sep = "_") %>% unique()

length(BIN.SP)

x <- "Arctica islandica_BOLD:AAJ7909"                

actu.width <- DNA.total.aligOK@ranges@width %>% unique()
actu.width

DNA.BIN.SP.filt <- DNAStringSet()

for(x in BIN.SP){
  print(x)
  DNA.int1 <- DNA.total.aligOK[str_detect(names(DNA.total.aligOK),x)]
  # Remove duplicates - still some since we recut everything
  DNA.int2 <- DNA.int1[!duplicated(DNA.int1)]
  # keep sequences with least missing nucleotides
  N.vec <- DNA.int2 %>% as.data.frame() %>% pull(x) %>% str_count("N")
  N.min <- N.vec %>% min()
  DNA.int3 <- DNA.int2[N.vec == N.min]
  # Remove sequences with missing nuc at the beginning, if enough sequences available  
  Mis.Start.vec <- DNA.int3 %>% subseq(start = 1, end = 10) %>%  as.data.frame() %>% pull(x) %>% str_count("-")
  Mis.Start.min <- Mis.Start.vec %>% min()
  DNA.int4 <- DNA.int3[Mis.Start.vec == Mis.Start.min]  
  # Remove sequences with missing nuc at the end, if enough sequences available  
  Mis.End.vec <- DNA.int4 %>% subseq(start = actu.width - 9, end = actu.width) %>%  as.data.frame() %>% pull(x) %>% str_count("-")
  Mis.End.min <- Mis.End.vec %>% min()
  DNA.int5 <- DNA.int4[Mis.End.vec == Mis.End.min]   
  
  DNA.BIN.SP.filt <- c(DNA.BIN.SP.filt,DNA.int5)
      
}

DNA.BIN.SP.filt

# Count the number of missing space
Mis.tot.vec <-  DNA.BIN.SP.filt %>%  as.data.frame() %>% pull(x) %>% str_count("-") 

Mis.tot.vec %>% table()
(actu.width - Mis.tot.vec) %>% table()

# remove when there is more than 50 missing values

DNA.BIN.SP.filt.final <- DNA.BIN.SP.filt[Mis.tot.vec < 50]

length(DNA.BIN.SP.filt)
length(DNA.BIN.SP.filt.final)


#writeXStringSet(DNA.BIN.SP.filt.final, filepath = "01_Raw_data/02_Sequences/02_Post-sequencing/postseq.DNA.align.cut.BINfilt.2021-01-20.fasta")

DNA.BIN.SP.filt.final <- readDNAStringSet("01_Raw_data/02_Sequences/02_Post-sequencing/postseq.DNA.align.cut.BINfilt.2021-01-20.fasta")


# Then, check within and without BIN distance

# Stats -------------------------------------------------------------------

DNA.total.df <- data.frame(
  id = names(DNA.BIN.SP.filt.final),
  seq = DNA.BIN.SP.filt.final %>% as.data.frame() %>% pull(x)
)

DNA.total.df$Species <- sapply(str_split(names(DNA.BIN.SP.filt.final), "_"),`[`,2) 
DNA.total.df$ID.BOLD <- sapply(str_split(names(DNA.BIN.SP.filt.final), "_"),`[`,1)
DNA.total.df$bin_uri <- sapply(str_split(names(DNA.BIN.SP.filt.final), "_"),`[`,3)
DNA.total.df$ID.SEQ <- sapply(str_split(names(DNA.BIN.SP.filt.final), "_"),`[`,4)
DNA.total.df$N.SEQ <- sapply(str_split(names(DNA.BIN.SP.filt.final), "_"),`[`,5) %>% str_remove("n")

head(DNA.total.df)



DNA.total.df %>% names()

DNA.total.df %>% nrow()
DNA.total.df %>% pull(bin_uri) %>% unique() %>% length()

DNA.total.df <- DNA.total.df %>% left_join(Meta.seq.total %>% select(processid, NWA_code),
                                           by = c("ID.BOLD" = "processid")) %>% 
                                 mutate(NWA_code = ifelse(NWA_code == 1, "NWA", "OTHER"),
                                        Taxo.level = ifelse(str_detect(Species, " "), "Species", "Genus")) %>% 
                                 left_join(cat.list, by = c("Species" = "Name"))

DNA.total.df %>% head()
# Number by cat

DNA.total.df %>% filter(Taxo.level == "Genus") %>% View()

DNA.total.df %>% group_by(Categorie, Taxo.level) %>% 
  summarise(N.SP = length(unique(Species)),
            N.Family = length(unique(family)),
            N.BIN = length(unique(bin_uri)))

DNA.total.df %>% group_by(Categorie, NWA_code) %>% 
  summarise(N.SP = length(unique(Species)),
            N.Family = length(unique(family)),
            N.BIN = length(unique(bin_uri)))

# The number of BIN per species
SP.BIN.tbl <- DNA.total.df %>% group_by(Categorie, Species) %>% 
                 summarise(N.BIN = length(unique(bin_uri)),
                           BIN = paste(unique(bin_uri), collapse = ", ")) %>% 
                 arrange(desc(N.BIN)) 

SP.BIN.tbl %>% filter(N.BIN > 1) %>% arrange (Categorie, Species, N.BIN) %>% View()
  write_csv("clipboard")


SP.BIN.tbl %>% group_by(N.BIN) %>% summarise(N = n())

SP.BIN.tbl %>% group_by(Categorie, N.BIN) %>% summarise(N = n())

# The number of Species by geocode

DNA.total.df %>% group_by(NWA_code, Taxo.level) %>% summarise(N = length(unique(Species))) 

DNA.total.df %>% group_by(Categorie, NWA_code, Taxo.level) %>% summarise(N = length(unique(Species)))

# The number of species by BIN

BIN.SP.tbl <- DNA.total.df %>% group_by(Categorie, bin_uri) %>% 
  summarise(N.SP = length(unique(Species)),
            Species = paste(unique(Species), collapse = ", ")) %>% 
  arrange(desc(N.SP)) 

BIN.SP.tbl

# Table for supp material
BIN.SP.tbl %>% filter(N.SP > 1) %>% select(-N.SP) #%>% write_csv("clipboard")

BIN.SP.tbl %>% group_by(Categorie, N.SP) %>% summarise(N = n())


# Compute distance between BIN
dist.res <- data.frame(Species = character(), 
                       BIN1 = character(), 
                       BIN2 = character(), 
                       dist.comp = character(), 
                       dist = numeric(),
                       stringsAsFactors = F)

# Check which SP are associated with each BIN
#sp.bin <- data.frame(Species = character(),
#                     bin = character(),
#                     species_name = character(),
#                     N = numeric(),
#                     stringsAsFactors = F)

for(x in DNA.total.df$Species %>% unique()){

print(x)  
  
DNA.sub <- DNA.total.df %>% filter(Species == x)

# compute stuff by bin, to keep at least 1 bin
vec.id <- vector()

for(y in DNA.sub$bin_uri %>% unique()){
  # Filter by bin
  DNA.sub.bin <- DNA.sub %>% filter(bin_uri == y)  

  #thr.start <- min(DNA.sub.bin$Start.pos)
  #DNA.sub.bin <- DNA.sub.bin %>% filter(Start.pos == thr.start)
 #
  #thr.stop <- max(DNA.sub.bin$Stop.pos)
  #DNA.sub.bin <- DNA.sub.bin %>% filter(Stop.pos == thr.stop)  
  #
  #thr.ambi <- min(DNA.sub.bin$N_nuc_ambi)
  #DNA.sub.bin <- DNA.sub.bin %>% filter(N_nuc_ambi == thr.ambi)  
    
  vec.id <- c(vec.id, DNA.sub.bin %>% pull(ID.BOLD))
  
  # Check species by bin
 # sp.bin.int <- bold_seqspec(bin = y)
  
  #sp.bin.int <- sp.bin.int %>% mutate(species_name = ifelse(species_name == "", "Unknown", species_name)) %>% 
  #  group_by(species_name) %>% summarise(N = n()) %>% 
  #  group_by() %>% mutate(bin = y, Species = x)
  
#  sp.bin <- bind_rows(sp.bin, sp.bin.int)
  
}

DNA.sub.filt <- DNA.sub %>% filter(ID.BOLD %in% vec.id)

if(nrow(DNA.sub.filt) > 1) {

  DNA.int <- DNAStringSet(DNA.sub.filt$seq)
  names(DNA.int) <- DNA.sub.filt$bin_uri# paste(DNA.sub.filt$bin_uri, DNA.sub.filt$GSL, sep = "_") 

  mat.dist <- dist.dna(as.DNAbin(DNA.int), model = "K80", as.matrix = T) 

  mat.dist[upper.tri(mat.dist, diag = T)] <- NA

  mat.df <- mat.dist %>% as.data.frame()
  mat.df$BIN1 <- dimnames(mat.dist)[[1]] 

res <- mat.df %>% pivot_longer(cols = names(mat.df) %>% str_subset("BIN1", negate = T),
                      names_to = "BIN2") %>% 
         filter(!is.na(value)) %>% 
         group_by(BIN1, BIN2) %>% 
         summarise(dist = mean(value, na.rm = T)) %>% 
         mutate(Species = x,
                dist.comp = ifelse(BIN1 == BIN2, "intra", "inter")) %>% 
         select(Species, BIN1, BIN2, dist.comp, dist)


dist.res <- bind_rows(dist.res, res)

}

}

dist.res %>% View()
#sp.bin %>% View()

dist.res %>% names()
dist.res %>% head()

# Stats for the manuscript
dist.res %>% subset(str_detect(Species, " ", negate = F)) %>% 
              group_by(dist.comp) %>% 
             summarise(mean = mean(dist),
                       min = min(dist),
                       max = max(dist),
                       N = n()) 
           

