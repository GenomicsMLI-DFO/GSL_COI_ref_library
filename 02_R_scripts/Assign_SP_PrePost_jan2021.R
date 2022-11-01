
# Info --------------------------------------------------------------------

# Assign species to haplotype tables
# 
# Audrey Bourret
# 2020-10-20
#

# Library -----------------------------------------------------------------

library(dada2); packageVersion("dada2")
#library(Biostrings)

library(DECIPHER); packageVersion("DECIPHER")

library(tidyverse)
library(readxl)


library(ggpubr)
library(ggtext)

# ESV Data --------------------------------------------------------------------

# Seqtab
load("./03_Results/01_ASV_creation/01_data/ALLtable.FINAL.data")

ASVtab.COI$SEQ <- row.names(ASVtab.COI)
ASVtab.COI %>% names()
ASVtab.COI.cor %>% names()


# List of samples with less than 200 reads
# (from Compute_Nreads.R)

BAD.samples <- c("PPO-025","PPO-026","PPO-027","PPO-028","PPO-112","PPO-121")

# Check for low frequency reads

ASVtab.COI %>% pivot_longer(names(.) %>% str_subset("PPO"), 
                            names_to = "ID_labo", 
                            values_to = "Nreads") %>% 
               group_by(ID) %>% 
               summarise(Nreads = sum(Nreads)) %>% 
               arrange(Nreads) %>% 
               ggplot(aes(x = Nreads)) + geom_histogram() +
               scale_x_continuous(limits = c(0,1000))


ASV.stats <- ASVtab.COI %>% pivot_longer(names(.) %>% str_subset("PPO"), 
                            names_to = "ID_labo", 
                            values_to = "Nreads") %>% 
              # Enlever les 6 échantillons louches
              filter(ID_labo %nin% BAD.samples) %>%   
              group_by(ID) %>% 
              summarise(Nreads = sum(Nreads)) %>% 
              arrange(Nreads) 

ASV.stats %>% filter(Nreads <=10) %>% pull(Nreads) %>% table()

# The number of sigleton
ASV.stats %>% filter(Nreads <=1) %>% nrow()

ASV.stats %>% nrow()
ASV.stats %>% tail(1000)

ASV.singleton <- ASV.stats %>% filter(Nreads <=1) %>% pull(ID)
ASV.singleton

# List of ASV that need to be analyze
ASV.sample <- ASVtab.COI.cor %>% pivot_longer(names(.) %>% str_subset("PPO"), 
                                          names_to = "ID_labo", 
                                          values_to = "Nreads") %>% 
              # Enlever les 6 échantillons louches
              filter(ID_labo %nin% BAD.samples,
                     ID %nin% ASV.singleton) %>% 
              group_by(ID) %>% 
              summarise(Nreads = sum(Nreads)) %>% 
              arrange(Nreads) %>% 
              filter(Nreads > 0) %>% pull(ID)


ASV.sample %>% length()
#

# Table of all ASV

ASV <- DNAStringSet(ASVtab.COI$SEQ)
names(ASV) <- ASVtab.COI$ID

ASV[1]

str(ASV)


# Peut-être réviser ça

ASV.meta <- data.frame(ID = names(ASV), width = ASV@ranges@width)

ASV.meta %>% filter(ID %in% ASV.sample) %>% summary()

# Primer
mICOIintF <- "GGWACWGGWTGAACWGTWTAYCCYCC"

# Taxonomy data ----------------------------------------------------------------

# Use the final file

Taxa.Seq.final <- read_csv2("03_Results/Taxa.Seq.final.csv")

# Add a "sp" to those assigned at the genus level

Taxa.Seq.final <-  Taxa.Seq.final %>% mutate(Species = ifelse(Level == "Genus", paste(Name, "sp"), Name))
Taxa.Seq.final %>% names()

# Check if there is missing taxo

Taxa.Seq.final %>% filter(is.na(kingdom))
Taxa.Seq.final %>% filter(is.na(phylum))
Taxa.Seq.final%>% filter(is.na(class))
Taxa.Seq.final %>% filter(is.na(order))

#Taxa.Seq.final %>% filter(is.na(order)) %>% View()

Taxa.Seq.final <- Taxa.Seq.final %>% mutate(order = ifelse(is.na(order), paste(family, "order_incertae_sedis", sep = "_"), as.character(order)))

Taxa.Seq.final %>% filter(is.na(family))
Taxa.Seq.final %>% filter(is.na(genus))
Taxa.Seq.final %>% filter(is.na(Species))


# Load taxonomy data

ncbi.tax <- readr::read_tsv(file.path("rankedlineage.dmp"), 
                     col_names = c("id", "name", "species", "genus", "family", "order", "class","phylum", "kingdom", "superkingdom"), 
                     col_types=("i-c-c-c-c-c-c-c-c-c-"))

ncbi.tax %>% pull(species) %>% unique() %>% head()



# Data Sample -------------------------------------------------------------

# dataseq <- read_excel(file.path(get.value("info.path"), "Samples_PPO_BA.xlsx"), 
#                       sheet = "Samples", na = "NA")
# dataseq
# 
# 
# datalab <- read_excel(file.path("../../Metabarcodage/02_Metabar_Marin_2019/MetaBAR_2019/01_Raw_data/00_FileInfos/Fichier_projet_ADNe_PPO.xlsx"),
#                       sheet = "Labo")
# 
# datalab$Site <- sapply(str_split(datalab$Station, "-"), `[`,1)
# datalab$Profondeur <- sapply(str_split(datalab$Station, "-"), `[`,2)
# 
# dataseq <- dataseq %>% left_join(datalab %>% select(-c(Type_echantillon, Neg_prelevement,Neg_filtration,Neg_extraction)), by = "ID_labo")
# 
# 
# ASV.read.byProjet <- ASVtab.COI.cor %>% pivot_longer(cols = names(ASVtab.COI) %>% str_subset("ID|SEQ", negate = TRUE),
#                                                      names_to = "ID_labo",
#                                                      values_to = "Nreads") %>% 
#   left_join(dataseq) %>% 
#   filter(Type_echantillon == "Echantillon") %>% 
#   group_by(ID, ID_projet) %>% summarise(Nreads = sum(Nreads)) %>%
#   pivot_wider(id_cols = ID, names_from = ID_projet, values_from = Nreads) %>% 
#   mutate(Nreads.tot = BA + PPO)
# 
# write_csv(ASV.read.byProjet, file = "03_Results/ESV_reads.csv")

# A priori database shorten -------------------------------------------

REF.STEP1 <- readDNAStringSet(file.path("./01_Raw_data/02_Sequences/01_Pre-sequencing","DNA.align.cut.BINfilt.2020-10-06.fasta"))

res.F <- vmatchPattern(DNAString(mICOIintF), REF.STEP1, max.mismatch = 11)
print(res.F@ends %>% unlist() %>% table())

REF.STEP1.short <- subseq(REF.STEP1, start = 355)

SPbin.list <- paste(sapply(str_split(names(REF.STEP1.short), "_"), `[`,2),
                    sapply(str_split(names(REF.STEP1.short), "_"), `[`,3),
                           sep = "_")%>% unique() 

REF.STEP1.short.nodup <- DNAStringSet()


for(x in SPbin.list){
  print(x)
  
  # Remove duplicates for this small part
  DNA <- REF.STEP1.short[str_detect(names(REF.STEP1.short), x)]
  DNA <- DNA[DNA %>% duplicated() == F]
  
  # Remove sequence that are of lesser quality
  
  quality <- data.frame(ID = names(DNA),
                        SEQ = DNA %>% as.data.frame() %>% pull(x)) %>% 
             mutate(N_nuc_A = str_count(SEQ, "A"),
                    N_nuc_C = str_count(SEQ, "C"),
                    N_nuc_T = str_count(SEQ, "T"),
                    N_nuc_G = str_count(SEQ, "G"),
                    N_nuc = N_nuc_A + N_nuc_C + N_nuc_T + N_nuc_G)
  
  max.trh <- quality$N_nuc %>% max()
  quality <- quality %>% filter(N_nuc == max.trh)
  
  DNA <- DNAStringSet(quality$SEQ) 
  names(DNA) <- quality$ID
  
  # Remove gaps
  
  DNA <- RemoveGaps(DNA,
             removeGaps = "all",
             processors = 1)
  
  REF.STEP1.short.nodup <- c(REF.STEP1.short.nodup, DNA)
  
  }

length(REF.STEP1.short)
length(REF.STEP1.short.nodup)

table(REF.STEP1.short.nodup@ranges@width)

REF.STEP1.short.nodup # 

taxid <- data.frame(SEQ = names(REF.STEP1.short.nodup),
                    ID = sapply(str_split(names(REF.STEP1.short.nodup), "_"),`[`,1),
                    SP = sapply(str_split(names(REF.STEP1.short.nodup), "_"),`[`,2),
                    stringsAsFactors = F)

taxid <- taxid %>% mutate(Species = ifelse(str_detect(SP, " "), SP, paste(SP, "sp"))) %>% 
  left_join(Taxa.Seq.final %>% unique(), by = c("Species")) %>% #head()
  mutate(Taxo = paste("Root", kingdom, phylum, class, order, family, genus, Species, sep = "; "),
         NewSEQ = paste(ID, Taxo))

names(REF.STEP1.short.nodup) <- taxid$NewSEQ


REF.STEP1.short.nodup


#writeXStringSet(REF.STEP1.short.nodup, file.path("./01_Raw_data/02_Sequences/03_Final-reference","REF.priori.Cut.withTAXO.2021-01-21.fasta"))


# A posteriori database shorten -------------------------------------------

#REF.STEP2 <- readDNAStringSet(file.path("./01_Raw_data/02_Sequences/02_Post-sequencing/","postseq.DNA.align.cut.BINfilt.2021-01-20.fasta"))
REF.STEP2 <- readDNAStringSet(file.path("./01_Raw_data/02_Sequences/03_Final-reference/","DNA.align.cut.BINfilt.2021-01-20.fasta"))


REF.STEP2

res.F <- vmatchPattern(DNAString(mICOIintF), REF.STEP2, max.mismatch = 11)
print(res.F@ends %>% unlist() %>% table())

REF.STEP2.short <- subseq(REF.STEP2, start = 355)

SPbin.list <- paste(sapply(str_split(names(REF.STEP2.short), "_"), `[`,2),
                    sapply(str_split(names(REF.STEP2.short), "_"), `[`,3),
                    sep = "_")%>% unique() 

REF.STEP2.short.nodup <- DNAStringSet()


for(x in SPbin.list){
  print(x)
  
  # Remove duplicates for this small part
  DNA <- REF.STEP2.short[str_detect(names(REF.STEP2.short), x)]
  DNA <- DNA[DNA %>% duplicated() == F]
  
  # Remove sequence that are of lesser quality
  
  quality <- data.frame(ID = names(DNA),
                        SEQ = DNA %>% as.data.frame() %>% pull(x)) %>% 
    mutate(N_nuc_A = str_count(SEQ, "A"),
           N_nuc_C = str_count(SEQ, "C"),
           N_nuc_T = str_count(SEQ, "T"),
           N_nuc_G = str_count(SEQ, "G"),
           N_nuc = N_nuc_A + N_nuc_C + N_nuc_T + N_nuc_G)
  
  max.trh <- quality$N_nuc %>% max()
  quality <- quality %>% filter(N_nuc == max.trh)
  
  DNA <- DNAStringSet(quality$SEQ) 
  names(DNA) <- quality$ID
  
  # Remove gaps
  
  DNA <- RemoveGaps(DNA,
                    removeGaps = "all",
                    processors = 1)
  
  REF.STEP2.short.nodup <- c(REF.STEP2.short.nodup, DNA)
  
}

length(REF.STEP2.short)
length(REF.STEP2.short.nodup)

table(REF.STEP2.short.nodup@ranges@width)

REF.STEP2.short.nodup # 

taxid <- data.frame(SEQ = names(REF.STEP2.short.nodup),
                    ID = sapply(str_split(names(REF.STEP2.short.nodup), "_"),`[`,1),
                    SP = sapply(str_split(names(REF.STEP2.short.nodup), "_"),`[`,2),
                    stringsAsFactors = F)

taxid <- taxid %>% mutate(Species = ifelse(str_detect(SP, " "), SP, paste(SP, "sp"))) %>% 
  left_join(Taxa.Seq.final %>% unique(), by = c("Species")) %>% #head()
  mutate(Taxo = paste("Root", kingdom, phylum, class, order, family, genus, Species, sep = "; "),
         NewSEQ = paste(ID, Taxo))


names(REF.STEP2.short.nodup) <- taxid$NewSEQ


REF.STEP2.short.nodup


#writeXStringSet(REF.STEP2.short.nodup, file.path("./01_Raw_data/02_Sequences/03_Final-reference","REF.posteriori.Cut.withTAXO.2021-01-21.fasta"))


# Rename full dataset too

taxid.ori <- data.frame(SEQ = names(REF.STEP2),
                    ID = sapply(str_split(names(REF.STEP2), "_"),`[`,1),
                    SP = sapply(str_split(names(REF.STEP2), "_"),`[`,2),
                    stringsAsFactors = F)

taxid.ori <- taxid.ori %>% mutate(Species = ifelse(str_detect(SP, " "), SP, paste(SP, "sp"))) %>% 
  left_join(Taxa.Seq.final %>% unique(), by = c("Species")) %>% #head()
  mutate(Taxo = paste("Root", kingdom, phylum, class, order, family, genus, Species, sep = "; "),
         NewSEQ = paste(ID, Taxo))

names(REF.STEP2) <- taxid.ori$NewSEQ

#writeXStringSet(REF.STEP2, file.path("./01_Raw_data/02_Sequences/03_Final-reference","GSL-rl_COI_Folmer650pb_450taxa_1304seq.fasta"))


# Data - Reference DB -----------------------------------------------------

DNA.priori.REF  <- readDNAStringSet(file.path("./01_Raw_data/02_Sequences/03_Final-reference","REF.priori.Cut.withTAXO.2021-01-21.fasta"))
DNA.postseq.REF <- readDNAStringSet(file.path("./01_Raw_data/02_Sequences/03_Final-reference","REF.posteriori.Cut.withTAXO.2021-01-21.fasta"))


# Peut-être réviser ça

names(DNA.postseq.REF)[1:10]

postseq.meta <- data.frame(ID = sapply(str_split(names(DNA.postseq.REF), " "), `[`, 1),
                           Name = names(DNA.postseq.REF),
                           width = DNA.postseq.REF@ranges@width)


# TaxID : Training the classifier -------------------------------------------------

# train the classifier

set.seed(111) # to get the same result everytimes

trainingSet.priori <- LearnTaxa(DNA.priori.REF,
                         names(DNA.priori.REF))

trainingSet.priori 

trainingSet.postseq <- LearnTaxa(DNA.postseq.REF,
                                names(DNA.postseq.REF))

trainingSet.postseq

set.seed(NULL)

# If there is problematic sequences, we should do something

plot(trainingSet.priori)
plot(trainingSet.postseq)

# TaxID : Assigning taxo ----------------------------------------------------------

# Very high confidence
ids.priori.60 <- DECIPHER::IdTaxa(ASV,
              trainingSet.priori,
              type="extended",
              strand="top",
              bootstraps = 100,
              threshold=60,
              processors=1)

# Very high confidence
ids.postseq.60 <- DECIPHER::IdTaxa(ASV,
                                  trainingSet.postseq,
                                  type="extended",
                                  strand="top",
                                  bootstraps = 100,
                                  threshold=60,
                                  processors=1)


# High confidence
ids.priori.50 <- DECIPHER::IdTaxa(ASV,
                                  trainingSet.priori,
                                  type="extended",
                                  strand="top",
                                  bootstraps = 100,
                                  threshold=50,
                                  processors=1)

# High confidence
ids.postseq.50 <- DECIPHER::IdTaxa(ASV,
                                  trainingSet.postseq,
                                  type="extended",
                                  strand="top",
                                  bootstraps = 100,
                                  threshold=50,
                                  processors=1)


#Moderate confidence
ids.priori.40 <- DECIPHER::IdTaxa(ASV,
                                  trainingSet.priori,
                                  type="extended",
                                  strand="top",
                                  bootstraps = 100,
                                  threshold=40,
                                  processors=1)

#Moderate confidence
ids.postseq.40 <- DECIPHER::IdTaxa(ASV,
                                  trainingSet.postseq,
                                  type="extended",
                                  strand="top",
                                  bootstraps = 100,
                                  threshold=40,
                                  processors=1)

add.ranks <- function(ids, ranks = c("root", "kingdom", "phylum", "class", "order", "family", "genus", "species")){
  for(x in 1:length(ids)){
    N <- length(ids[[x]]$taxon)
    ids[[x]]$rank <- ranks[1:N]
  }

  return(ids)
} # End of add.ranks function

ids.priori.60 <- add.ranks(ids.priori.60)
ids.priori.50 <- add.ranks(ids.priori.50)
ids.priori.40 <- add.ranks(ids.priori.40)

ids.postseq.60 <- add.ranks(ids.postseq.60)
ids.postseq.50 <- add.ranks(ids.postseq.50)
ids.postseq.40 <- add.ranks(ids.postseq.40)

#save(file = get.value("IDT.TS.data"), 
#     list = c("trainingSet.priori", "trainingSet.postseq", 
#              "ids.priori.40", "ids.priori.50","ids.priori.60",
#              "ids.postseq.40", "ids.postseq.50","ids.postseq.60")) 


#load(get.value("IDT.TS.data"))


#save(file = "Training_COI_2021-02-08.data" , list = "trainingSet.postseq")

###

ids.to.tbl <- function(ids, ranks = c("root", "kingdom", "phylum", "class", "order", "family", "genus", "species")){
  ids.tbl <- t(sapply(ids, function(x) {
                  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
})) %>% as.data.frame()
  
  names(ids.tbl) <- ranks
  return(ids.tbl)
}

add.ASV <- function(ids.tbl, ASVtab){
  ids.tbl$ID <- ASVtab$ID
  ids.tbl$SEQ <- ASVtab$SEQ
  
  return(ids.tbl)
}

add.taxon <- function(ids.tbl){
  ids.tbl <- ids.tbl %>% mutate(Taxon = ifelse(!is.na(species), species %>% as.character(),
                                        ifelse(!is.na(genus), genus %>% as.character(),
                                        ifelse(!is.na(family), family %>% as.character(),
                                        ifelse(!is.na(order), order %>% as.character(),
                                        ifelse(!is.na(class), class %>% as.character(),
                                        ifelse(!is.na(phylum), phylum %>% as.character(),
                                        ifelse(!is.na(kingdom), kingdom %>% as.character(), "Unassigned"))))))),
                                Levels = ifelse(!is.na(species), "species",
                                         ifelse(!is.na(genus), "genus",
                                         ifelse(!is.na(family), "family",
                                         ifelse(!is.na(order), "order",
                                         ifelse(!is.na(class), "class",
                                         ifelse(!is.na(phylum), "phylum",
                                         ifelse(!is.na(kingdom), "kingdom", NA))))))),
                                Levels = factor(Levels, levels = c("species", "genus", "family", "order", "class", "phylum", "kingdom"))
                                
                                )

  return(ids.tbl)
}

ids.tbl.priori.60 <- ids.to.tbl(ids.priori.60) %>% add.ASV(ASVtab.COI) %>% add.taxon() %>% mutate(threshold = 60, DB = "priori")
ids.tbl.priori.50 <- ids.to.tbl(ids.priori.50) %>% add.ASV(ASVtab.COI) %>% add.taxon() %>% mutate(threshold = 50, DB = "priori")
ids.tbl.priori.40 <- ids.to.tbl(ids.priori.40) %>% add.ASV(ASVtab.COI) %>% add.taxon() %>% mutate(threshold = 40, DB = "priori")

ids.tbl.postseq.60 <- ids.to.tbl(ids.postseq.60) %>% add.ASV(ASVtab.COI) %>% add.taxon() %>% mutate(threshold = 60, DB = "postseq")
ids.tbl.postseq.50 <- ids.to.tbl(ids.postseq.50) %>% add.ASV(ASVtab.COI) %>% add.taxon() %>% mutate(threshold = 50, DB = "postseq")
ids.tbl.postseq.40 <- ids.to.tbl(ids.postseq.40) %>% add.ASV(ASVtab.COI) %>% add.taxon() %>% mutate(threshold = 40, DB = "postseq")


ids.tbl.all <- bind_rows(ids.tbl.priori.60, ids.tbl.priori.50, ids.tbl.priori.40,
                         ids.tbl.postseq.60, ids.tbl.postseq.50, ids.tbl.postseq.40) %>% 
               mutate(method = "IDtaxa")

#save(file = get.value("IDT.TS.data"), 
#     list = c("trainingSet.priori", "trainingSet.postseq", 
#              "ids.priori.40", "ids.priori.50","ids.priori.60",
#              "ids.postseq.40", "ids.postseq.50","ids.postseq.60", "ids.tbl.all")) 


# DATA - IDTAXA results ---------------------------------------------------

load(get.value("IDT.TS.data"))


# Stats

ids.tbl.all %>% filter(ID %in% ASV.sample) %>% 
                   group_by(DB,Levels) %>% 
                   summarise(N40 = length(Taxon[threshold == 40]),
                             N50 = length(Taxon[threshold == 50]),
                             N60 = length(Taxon[threshold == 60]))

ids.tbl.all %>% filter(ID %in% ASV.sample,
                          !is.na(Levels)) %>% 
  group_by(DB, Levels, Taxon) %>% 
  summarise(N40 = length(Taxon[threshold == 40]),
            N50 = length(Taxon[threshold == 50]),
            N60 = length(Taxon[threshold == 60])) %>% View()


ids.tbl.all %>% left_join(ASV.read.byProjet) %>% filter(Taxon != "Unassigned", threshold == 60) %>% 
  group_by(phylum, Taxon, DB) %>% summarise(PPO = sum(PPO), BA = sum(BA)) %>% 
  mutate(Ntot = PPO + BA,
         DB = factor(DB, levels =c("priori", "postseq"))) %>% 
  filter(Ntot > 0) %>% arrange(desc(Ntot)) %>% 
  select(-Ntot) %>%  
  pivot_longer(cols = c(PPO, BA), names_to = "ID_projet", values_to = "Nreads") %>% 
  ggplot(aes(x = DB, y = Taxon, fill = Nreads)) + 
  geom_bin2d(col = "gray") +
  scale_fill_distiller(palette = "Spectral", trans = "log10", na.value = "white") +
  guides(fill = guide_colourbar(title = "N reads", title.hjust = 0)) + 
  theme_bw()+
  facet_grid(phylum ~ ID_projet, scale = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0)) #+ coord_flip()

# Phyloseq test -----------------------------------------------------------

library(phyloseq)


ids.work <- ids.tbl %>% left_join(ASV.read.byProjet) %>% filter(Taxon != "Unassigned") %>% 
  mutate(Ntot = PPO + BA) %>% 
  filter(Ntot > 0) %>% arrange(desc(Ntot)) %>% 
  select(-Ntot) 

otumat <- ids.work %>% select(PPO, BA) %>% as.matrix() 

ids.work <- ids.tbl %>% left_join(ASVtab.COI.cor) %>% filter(Taxon != "Unassigned") 

otumat <-  ids.work %>% select(names(ASVtab.COI.cor) %>% str_subset("ID|SEQ", negate = T)) %>% as.matrix() 

row.names(otumat) <- ids.work %>% pull(ID)

taxmat <- ids.work %>% select(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% as.matrix()
row.names(taxmat) <- ids.work %>% pull(ID)


# transforms the matrix in the phyloseq format

OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
OTU
TAX

datmat <- dataseq %>% filter(ID_labo %in% colnames(otumat)) %>% as.data.frame()
row.names(datmat) <- dataseq %>% filter(ID_labo %in% colnames(otumat)) %>% pull(ID_labo)

DAT <- sample_data(datmat)

physeq = phyloseq(OTU, TAX, phy_tree(fitGTR$tree), DAT)
physeq

subset_samples(physeq, Type_echantillon == "Echantillon" & ID_projet == "PPO")

plot_bar(physeq, fill = "Phylum")

plot_bar(subset_samples(physeq, Type_echantillon == "Echantillon" & ID_projet == "PPO"), 
         x= "Profondeur", fill = "Phylum") + facet_grid(Site~., scale = "free")

plot_tree(physeq, color = "Phylum", shape = "Sample") + coord_polar(theta="y")

phy_tree(physeq)


data(GlobalPatterns)
str(GlobalPatterns)

# Create a tree
library(msa)

seqs <- ids.work %>% pull(SEQ)
names(seqs) <- ids.work %>% pull(ID) # This propagates to the tip labels of the tree
mult <- msa(seqs, method="ClustalW", type="dna", order="input")

library("phangorn")
phang.align <- as.phyDat(mult, type="DNA", names=ids.work %>% pull(SEQ))
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)



# BLAST function ---------------------------------------------------------------

# Get the right emplacement for blast modules
if(get_os() == "linux"){
  
  makeblastdb <- "makeblastdb"
  blastn <- "blastn"
  
} else{
  
  makeblastdb <- get.value("makeblastdb")
  blastn <- get.value("blastn")
  
}

# Performed LCA at a given threshold, and return everything
BLAST_LCA <- function(RES, threshold = 0.97){
  DF <- data.frame()
  
  RES.OK <- RES %>% filter(AlignmentLength >= .95 * width,
                           Identity >= threshold) %>% 
    mutate(Taxon = NA,
           Levels = NA)
  
  # if from ncbi
  if(str_count(names(RES.OK), "SciName") %>% sum() == 1){
    RES.OK <- RES.OK %>% filter(str_detect(SciName, "environmental sample|uncultured|predicted", negate = T))
    RES.OK$species <- paste(sapply(str_split(RES.OK$SciName, " "),`[`,1),
                            sapply(str_split(RES.OK$SciName, " "),`[`,2))
    RES.OK$genus <- sapply(str_split(RES.OK$specie, " "),`[`,1)
    RES.OK <- RES.OK %>% mutate(species = ifelse(str_detect(species, " sp[.]| cf[.]| aff."), NA, species))
  }
  
  ASV <- RES.OK %>% pull(QueryAccVer) %>% unique()
  for(x in seq_along(ASV)){
    RES.INT <- RES.OK %>% filter(QueryAccVer == ASV[x])
    
    # loops around ranks
    for(y in c("species", "genus", "family", "order", "class", "phylum", "kingdom")) {
      
      N.LCA <- RES.INT %>% filter(!is.na(y)) %>%  pull(y) %>% unique() %>% str_subset("NA", negate = T) %>% length()
      
      if(N.LCA==1){
        RES.INT$Taxon <-  RES.INT %>% filter(!is.na(y)) %>% pull(y) %>% unique()  %>% str_subset("NA", negate = T)
        RES.INT$Levels <-  y
        break;
      }else{
        RES.INT[,y] <- NA
      }
    } # END of the loop around columns
    #RES.INT <- RES.INT %>% select(QueryAccVer, kingdom, phylum, class, order, family, genus, species, Taxon, Levels) %>% 
    #                     distinct(.keep_all = T)
    
    if(nrow(RES.INT[which(!is.na(RES.INT[,y])),])>0){
      DF <- bind_rows(DF, RES.INT[which(!is.na(RES.INT[,y])),])      
    }

    
  }
  return(DF)  
  
}   

BLAST_TOPHIT <- function(RES, threshold = 0.95){
  DF <- data.frame()
  
  RES.OK <- RES %>% filter(AlignmentLength >= .95 * width,
                           Identity >= threshold) %>% 
    mutate(Taxon = NA,
           Levels = NA)
  
  # if from ncbi
  if(str_count(names(RES.OK), "SciName") %>% sum() == 1){
    RES.OK <- RES.OK %>% filter(str_detect(SciName, "environmental sample|uncultured|predicted", negate = T))
    RES.OK$species <- paste(sapply(str_split(RES.OK$SciName, " "),`[`,1),
                            sapply(str_split(RES.OK$SciName, " "),`[`,2))
    RES.OK$genus <- sapply(str_split(RES.OK$specie, " "),`[`,1)
    RES.OK <- RES.OK %>% mutate(species = ifelse(str_detect(species, " sp[.]| cf[.]| aff."), NA, species))
  }
  
  ASV <- RES.OK %>% pull(QueryAccVer) %>% unique()
  for(x in seq_along(ASV)){
    RES.INT <- RES.OK %>% filter(QueryAccVer == ASV[x])
    
    evalue.min <- min(RES.INT$evalue)
    RES.INT <- RES.INT %>% filter(evalue == evalue.min)
    
    #identity.max <- max(RES.INT$Identity)
    #RES.INT <- RES.INT %>% filter(Identity == identity.max)
    
    # loops around ranks
    for(y in c("species", "genus", "family", "order", "class", "phylum", "kingdom")) {
      
      N.LCA <- RES.INT %>% filter(!is.na(y)) %>%  pull(y) %>% unique() %>% str_subset("NA", negate = T) %>% length()
      
      if(N.LCA==1){
        RES.INT$Taxon <-  RES.INT %>% filter(!is.na(y)) %>% pull(y) %>% unique()  %>% str_subset("NA", negate = T)
        RES.INT$Levels <-  y
        break;
      }else{
        RES.INT[,y] <- NA
      }
    } # END of the loop around columns
    #RES.INT <- RES.INT %>% select(QueryAccVer, kingdom, phylum, class, order, family, genus, species, Taxon, Levels) %>% 
    #                     distinct(.keep_all = T)
    
    if(nrow(RES.INT[which(!is.na(RES.INT[,y])),])>0){
      DF <- bind_rows(DF, RES.INT[which(!is.na(RES.INT[,y])),])      
    }
    
  }
  return(DF)  
  
}   

sum.BLAST <- function(DF){
  RES <- DF %>% select(QueryAccVer, Taxon, Levels, species, genus, family, order, class, phylum, kingdom) %>% unique() 
  
  return(RES)
}


makeblastdb
blastn



# Blast - ASV on priori ---------------------------------------------------

# Create DB with makeblastdb

cmd <- paste("-in", file.path("./01_Raw_data/02_Sequences/01_Pre-sequencing","DNA.Cut.withTAXO.2020-10-21.fasta"),
             "-dbtype", "nucl",
             "-parse_seqids", 
             sep = " ") # forward adapter

#system2(makeblastdb, cmd, stdout=T, stderr=T) 


# Blast


cmd1 <- paste("-db", file.path("./01_Raw_data/02_Sequences/01_Pre-sequencing","DNA.Cut.withTAXO.2020-10-21.fasta"), # Already set in the environment to use taxid
              #"-db", file.path(get.value("ref.path"),"Blast","nt"),
              "-query",  file.path("./01_Raw_data/03_eDNA/03g_ASV","all.COI_ASV.fasta"),
              
              #"-taxids",file.path("C:\\Users\\bourreta\\Documents\\Projets\\Banque_REF_NCBI.git\\01_Data\\01_blastdb","taxdb"), 
              #"-query",  "ASV.12S.db.fasta",
              "-outfmt", 7, #"\"7 qseqid sacc staxid ssciname sskingdom pident length mismatch gapopen qstart qend sstart send evalue bitscore\"",
              "-out", file.path("./03_Results/02_TaxoAssign/01_Priori","Blast.Priori.COI_ASV_2020-10-28.out"), 
              "-perc_identity", 95,
              "-num_threads", 1,
              #"-max_target_seqs", 10, 
              sep = " ")# forward adapter

#A<-system2(blastn, cmd1, stdout=T, stderr=T)

RES.COI.priori <- read.table(file.path("./03_Results/02_TaxoAssign/01_Priori","Blast.Priori.COI_ASV_2020-10-28.out"),
                            sep="\t")


RES.COI.priori %>% head()
names(RES.COI.priori) <- c("QueryAccVer", "SubjectAccVer", "Identity", "AlignmentLength", "mismatches", "gap opens", "q. start", "q. end", "s. start", "s. end", "evalue", "bit score")

RES.COI.priori <- RES.COI.priori %>% 
  left_join(taxid %>% select(ID, species = Species, genus, family, order, class, phylum, kingdom ), 
            by = c("SubjectAccVer" = "ID")) %>% 
  left_join(ASV.meta, by = c("QueryAccVer" = "ID")) %>% 
  #left_join(ids.tbl %>% select(ID, Taxon),
  #          by = c("QueryAccVer" = "ID")) %>% 
  filter(QueryAccVer %in% ASV.sample)


RES.COI.priori %>% filter(AlignmentLength >= .95 * width) %>% pull(QueryAccVer) %>% unique() %>% length()


# FILT blast RES

COI.priori.BlastLCA.95 <- BLAST_LCA(RES.COI.priori, threshold = 95) 
COI.priori.BlastLCA.97 <- BLAST_LCA(RES.COI.priori, threshold = 97) 
COI.priori.BlastLCA.99 <- BLAST_LCA(RES.COI.priori, threshold = 99) 

COI.priori.BlastTOP.95 <- BLAST_TOPHIT(RES.COI.priori, threshold = 95) 
COI.priori.BlastTOP.97 <- BLAST_TOPHIT(RES.COI.priori, threshold = 97) 
COI.priori.BlastTOP.99 <- BLAST_TOPHIT(RES.COI.priori, threshold = 99) 


#save(file = file.path("./03_Results/02_TaxoAssign/01_Priori", "BLAST.priori.data"), 
#     list = c("RES.COI.priori",
#              "COI.priori.BlastLCA.95", "COI.priori.BlastLCA.97","COI.priori.BlastLCA.99",
#              "COI.priori.BlastTOP.95", "COI.priori.BlastTOP.97","COI.priori.BlastTOP.99")) 

#load(file.path("./03_Results/02_TaxoAssign/01_Priori", "BLAST.priori.data"))


#DF <- COI.priori.BlastLCA.95 

sum.BLAST <- function(DF){
  RES <- DF %>% select(QueryAccVer, Taxon, Levels, species, genus, family, order, class, phylum, kingdom) %>% unique() 
  
  return(RES)
}


COI.priori.BlastLCA.97 %>% sum.BLAST() %>% group_by(Levels) %>% summarise(N = n())
COI.priori.BlastLCA.99 %>% sum.BLAST() %>% group_by(Levels) %>% summarise(N = n())


# Blast - ASV on postseq ---------------------------------------------------

# Create DB with makeblastdb

cmd <- paste("-in", file.path("./01_Raw_data/02_Sequences/03_Final-reference","REF.posteriori.Cut.withTAXO.2021-01-21.fasta"),
             "-dbtype", "nucl",
             "-parse_seqids", 
             sep = " ") # forward adapter

system2(makeblastdb, cmd, stdout=T, stderr=T) 


# Blast


cmd1 <- paste("-db",file.path("./01_Raw_data/02_Sequences/03_Final-reference","REF.posteriori.Cut.withTAXO.2021-01-21.fasta"), # Already set in the environment to use taxid
              #"-db", file.path(get.value("ref.path"),"Blast","nt"),
              "-query",  file.path("./01_Raw_data/03_eDNA/03g_ASV","all.COI_ASV.fasta"),
              
              #"-taxids",file.path("C:\\Users\\bourreta\\Documents\\Projets\\Banque_REF_NCBI.git\\01_Data\\01_blastdb","taxdb"), 
              #"-query",  "ASV.12S.db.fasta",
              "-outfmt", 7, #"\"7 qseqid sacc staxid ssciname sskingdom pident length mismatch gapopen qstart qend sstart send evalue bitscore\"",
              "-out", file.path("./03_Results/02_TaxoAssign/02_Post-sequencing","Blast.postseq.COI_ASV_2021-01-21.out"), 
              "-perc_identity", 95,
              "-num_threads", 1,
              #"-max_target_seqs", 10, 
              sep = " ")# forward adapter

#A<-system2(blastn, cmd1, stdout=T, stderr=T)

RES.COI.postseq <- read.table(file.path("./03_Results/02_TaxoAssign/02_Post-sequencing","Blast.postseq.COI_ASV_2021-01-21.out"),
                             sep="\t")


RES.COI.postseq %>% head()
names(RES.COI.postseq) <- c("QueryAccVer", "SubjectAccVer", "Identity", "AlignmentLength", "mismatches", "gap opens", "q. start", "q. end", "s. start", "s. end", "evalue", "bit score")

RES.COI.postseq <- RES.COI.postseq %>% 
  left_join(taxid %>% select(ID, species = Species, genus, family, order, class, phylum, kingdom ), 
            by = c("SubjectAccVer" = "ID")) %>% 
  left_join(ASV.meta, by = c("QueryAccVer" = "ID")) %>% 
  #left_join(ids.tbl %>% select(ID, Taxon),
  #          by = c("QueryAccVer" = "ID")) %>% 
  filter(QueryAccVer %in% ASV.sample)


RES.COI.postseq %>% filter(AlignmentLength >= .95 * width) %>% pull(QueryAccVer) %>% unique() %>% length()


# FILT blast RES

COI.postseq.BlastLCA.95 <- BLAST_LCA(RES.COI.postseq, threshold = 95) 
COI.postseq.BlastLCA.97 <- BLAST_LCA(RES.COI.postseq, threshold = 97) 
COI.postseq.BlastLCA.99 <- BLAST_LCA(RES.COI.postseq, threshold = 99) 

COI.postseq.BlastTOP.95 <- BLAST_TOPHIT(RES.COI.postseq, threshold = 95) 
COI.postseq.BlastTOP.97 <- BLAST_TOPHIT(RES.COI.postseq, threshold = 97) 
COI.postseq.BlastTOP.99 <- BLAST_TOPHIT(RES.COI.postseq, threshold = 99) 


#save(file = file.path("./03_Results/02_TaxoAssign/02_Post-sequencing", "BLAST.postseq.data"), 
#     list = c("RES.COI.postseq",
#              "COI.postseq.BlastLCA.95", "COI.postseq.BlastLCA.97","COI.postseq.BlastLCA.99",
#              "COI.postseq.BlastTOP.95", "COI.postseq.BlastTOP.97","COI.postseq.BlastTOP.99")) 

#load(file.path("./03_Results/02_TaxoAssign/02_Post-sequencing", "BLAST.postseq.data"))


#DF <- COI.priori.BlastLCA.95 

sum.BLAST <- function(DF){
  RES <- DF %>% select(QueryAccVer, Taxon, Levels, species, genus, family, order, class, phylum, kingdom) %>% unique() 
  
  return(RES)
}


COI.postseq.BlastLCA.97 %>% sum.BLAST() %>% group_by(Levels) %>% summarise(N = n())
COI.postseq.BlastLCA.99 %>% sum.BLAST() %>% group_by(Levels) %>% summarise(N = n())


# Blast - ASV on NCBI nt -----------------------------------------------------

cmd1 <- paste("-db", file.path("nt"), # Already set in the environment to use taxid
              #"-db", file.path(get.value("ref.path"),"Blast","nt"),
              "-query",  file.path("./01_Raw_data/03_eDNA/03g_ASV","all.COI_ASV.fasta"),
              
              #"-taxids",file.path("C:\\Users\\bourreta\\Documents\\Projets\\Banque_REF_NCBI.git\\01_Data\\01_blastdb","taxdb"), 
              #"-query",  "ASV.12S.db.fasta",
              "-outfmt", "\"7 qseqid sacc staxid ssciname sskingdom pident length mismatch gapopen qstart qend sstart send evalue bitscore\"",
              "-out", file.path("./03_Results/02_TaxoAssign/03_NCBI_nt","Blast.all.COI_ASV_2020-10-26.out"), 
              "-perc_identity", 95,
              "-num_threads", 1,
              #"-max_target_seqs", 10, 
              sep = " ")# forward adapter

#A<-system2(blastn, cmd1, stdout=T, stderr=T)


RES.COI.ncbi <- read.table(file.path("./03_Results/02_TaxoAssign/03_NCBI_nt","Blast.all.COI_ASV_2020-10-26.out"),
                      sep="\t")


names(RES.COI.ncbi) <- c("QueryAccVer", "SubjectAccVer", "TaxoId","SciName", "SKindom", "Identity", "AlignmentLength", "mismatches", "gap opens", "q. start", "q. end", "s. start", "s. end", "evalue", "bit score")

RES.COI.ncbi <- RES.COI.ncbi %>% left_join(ncbi.tax, by = c("TaxoId" = "id")) %>% 
                       left_join(ASV.meta, by = c("QueryAccVer" = "ID")) %>% 
                       #left_join(ids.tbl %>% select(ID, Taxon),
                        #         by = c("QueryAccVer" = "ID")) %>% 
                       # Juste analyser les ASV qui nous intéresse
                        filter(QueryAccVer %in% ASV.sample) #%>% 
                        #mutate(species = SciName)

RES.COI.ncbi %>% View()

COI.ncbi.BlastLCA.95 <- BLAST_LCA(RES.COI.ncbi, threshold = 95) 
COI.ncbi.BlastLCA.97 <- BLAST_LCA(RES.COI.ncbi, threshold = 97) 
COI.ncbi.BlastLCA.99 <- BLAST_LCA(RES.COI.ncbi, threshold = 99) 

COI.ncbi.BlastTOP.95 <- BLAST_TOPHIT(RES.COI.ncbi, threshold = 95) 
COI.ncbi.BlastTOP.97 <- BLAST_TOPHIT(RES.COI.ncbi, threshold = 97) 
COI.ncbi.BlastTOP.99 <- BLAST_TOPHIT(RES.COI.ncbi, threshold = 99) 

#save(file = file.path("./03_Results/02_TaxoAssign/03_NCBI_nt", "BLAST.ncbi.data"), 
#     list = c("RES.COI.ncbi",
#              "COI.ncbi.BlastLCA.95", "COI.ncbi.BlastLCA.97","COI.ncbi.BlastLCA.99",
#              "COI.ncbi.BlastTOP.95", "COI.ncbi.BlastTOP.97","COI.ncbi.BlastTOP.99")) 

#load(file.path("./03_Results/02_TaxoAssign/03_NCBI_nt", "BLAST.ncbi.data"))



TEST %>% View()

RES.COI.ncbi %>% View()
RES.COI.ncbi %>% pull(QueryAccVer) %>% unique() %>% length()


RES.COI.ncbi %>% filter(AlignmentLength >= .95*width) %>% View()
RES.COI.ncbi %>% filter(AlignmetLength >= .95 * width) %>% pull(QueryAccVer) %>% unique() %>% length()

RES.COI.ncbi %>% names()

# Quels superkingdom

RES.COI.ncbi %>% filter(AlignmentLength >= .95 * width) %>%  pull(superkingdom) %>% table()


RES.COI.ncbi %>% filter(AlignmentLength >= .95 * width) %>% 
            group_by(superkingdom,kingdom,phylum,class) %>% 
            summarise(N = length(unique(QueryAccVer))) %>% View() #write_csv("clipboard")

RES.COI.ncbi %>% filter(AlignmentLength >= .95 * width,
                   kingdom == "Metazoa") %>% pull(QueryAccVer) %>% unique() %>% length()

RES.COI.ncbi %>% filter(AlignmentLength >= .95 * width,
                   kingdom == "Metazoa",
                   Taxon == "Unassigned") %>% View()


RES.COI.ncbi %>% filter(AlignmentLength >= .95 * width,
                        kingdom == "Metazoa",
                        Taxon == "Unassigned",
                        Identity >= 97) %>% pull(SciName) %>% unique()

# Select species for revision ---------------------------------------------

REV.ncbi <- RES.COI.ncbi %>% filter(AlignmentLength >= .95 * width,
                        kingdom == "Metazoa",
                        str_detect(SciName, "environmental sample|uncultured|predicted", negate = T)) %>% 
                 select(QueryAccVer, SciName, phylum, class, order, family, genus) %>% 
                 mutate(species = ifelse(str_detect(SciName, " cf[.]| aff."),
                                        paste(sapply(str_split(as.character(SciName), " "),`[`,1),
                                              sapply(str_split(as.character(SciName), " "),`[`,2),
                                              sapply(str_split(as.character(SciName), " "),`[`,3)
                                              ),
                         #normal sp name
                         paste(sapply(str_split(as.character(SciName), " "),`[`,1),
                               sapply(str_split(as.character(SciName), " "),`[`,2)
                         )
                        )) %>% select(-c(SciName)) %>% 
                              arrange(phylum, class, order, family, genus) %>%  unique() #%>% View()

REV.priori <- RES.COI.priori %>% filter(AlignmentLength >= .95 * width,
                          kingdom == "Animalia") %>% 
                  select(QueryAccVer, species) %>% 
                  unique() %>% 
                  group_by(QueryAccVer) %>% 
                  summarise(priori.assign = paste(species, collapse = "; "))

REV.ncbi %>% left_join(REV.priori) %>% select(-c(QueryAccVer)) %>% unique() %>% write_csv("./03_Results/02_TaxoAssign/01_Priori/Verif_NCBI_blast_2020-11-06.csv")

REV.ncbi %>% filter(order == "Alcyonacea") %>%  View()

filter(str_detect(SciName, "environmental sample|uncultured|predicted", negate = T))
RES.OK$species <- paste(sapply(str_split(RES.OK$SciName, " "),`[`,1),
                        sapply(str_split(RES.OK$SciName, " "),`[`,2))
RES.OK$genus <- sapply(str_split(RES.OK$specie, " "),`[`,1)
RES.OK <- RES.OK %>% mutate(species = ifelse(str_detect(species, " sp[.]| cf[.]| aff."), NA, species))



RES.COI.ncbi %>% names()

# Blast - postseq on NCBI nt -----------------------------------------------------

cmd1 <- paste("-db", file.path("nt"), # Already set in the environment to use taxid
              #"-db", file.path(get.value("ref.path"),"Blast","nt"),
              "-query",  file.path("./01_Raw_data/02_Sequences/03_Final-reference/","REF.posteriori.Cut.withTAXO.2021-01-21.fasta"),
              
              #"-taxids",file.path("C:\\Users\\bourreta\\Documents\\Projets\\Banque_REF_NCBI.git\\01_Data\\01_blastdb","taxdb"), 
              #"-query",  "ASV.12S.db.fasta",
              "-outfmt", "\"7 qseqid sacc staxid ssciname sskingdom pident length mismatch gapopen qstart qend sstart send evalue bitscore\"",
              "-out", file.path("./03_Results/02_TaxoAssign/03_NCBI_nt","Blast.all.COI_postseqONncbi_2021-03-09.out"), 
              "-perc_identity", 95,
              "-num_threads", 1,
              #"-max_target_seqs", 10, 
              sep = " ")# forward adapter

#A<-system2(blastn, cmd1, stdout=T, stderr=T)

RES.COI.postONncbi <- read.table(file.path("./03_Results/02_TaxoAssign/03_NCBI_nt","Blast.all.COI_postseqONncbi_2021-03-09.out"),
                           sep="\t")


names(RES.COI.postONncbi) <- c("QueryAccVer", "SubjectAccVer", "TaxoId","SciName", "SKindom", "Identity", "AlignmentLength", "mismatches", "gap opens", "q. start", "q. end", "s. start", "s. end", "evalue", "bit score")

RES.COI.postONncbi %>% head()

RES.COI.postONncbi <- RES.COI.postONncbi %>% left_join(ncbi.tax, by = c("TaxoId" = "id")) %>% 
  left_join(postseq.meta, by = c("QueryAccVer" = "ID")) %>% 
  #left_join(ids.tbl %>% select(ID, Taxon),
  #         by = c("QueryAccVer" = "ID")) %>% 
  # Juste analyser les ASV qui nous intéresse
  #filter(QueryAccVer %in% ASV.sample) #%>% 
  mutate(species = SciName,
         Taxon = sapply(str_split(Name, "; "), `[`, 8))

RES.COI.postONncbi %>% View()


COI.postONncbi.BlastLCA.95 <- BLAST_LCA(RES.COI.postONncbi, threshold = 95) 
COI.postONncbi.BlastLCA.97 <- BLAST_LCA(RES.COI.postONncbi, threshold = 97) 
COI.postONncbi.BlastLCA.99 <- BLAST_LCA(RES.COI.postONncbi, threshold = 99) 

COI.postONncbi.BlastTOP.95 <- BLAST_TOPHIT(RES.COI.postONncbi, threshold = 95) 
COI.postONncbi.BlastTOP.97 <- BLAST_TOPHIT(RES.COI.postONncbi, threshold = 97) 
COI.postONncbi.BlastTOP.99 <- BLAST_TOPHIT(RES.COI.postONncbi, threshold = 99) 

#save(file = file.path("./03_Results/02_TaxoAssign/03_NCBI_nt", "BLAST.postONncbi.data"), 
#    list = c("RES.COI.postONncbi",
#              "COI.postONncbi.BlastLCA.95", "COI.postONncbi.BlastLCA.97","COI.postONncbi.BlastLCA.99",
#              "COI.postONncbi.BlastTOP.95", "COI.postONncbi.BlastTOP.97","COI.postONncbi.BlastTOP.99")) 


COI.postONncbi.BlastLCA.95 %>% sum.BLAST() %>% group_by(Levels) %>% summarise(N = n())

