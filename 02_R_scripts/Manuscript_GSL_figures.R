# Info --------------------------------------------------------------------

# Figures, tables and stats related to the manuscript
# 
# Audrey Bourret
# 2021-10-04
#

# Library -----------------------------------------------------------------

library(here)
library(tidyverse)

library(ggpubr)
library(ggtext)

`%nin%` = Negate(`%in%`)


colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

pie(rep(1, 8), col = colorBlindBlack8)

# Real data - Dataset --------------------------------------------------------------------

#load(file.path("01_Raw_data/MS_All.assignments.Rdata"))
#
#write_csv(RES.total, file.path(here::here(), "01_Raw_data", "MS_AllAssignments.csv"))
#

# Assignment results

# N read observed for the PPO and BA at each ASV
ESV.reads <- readr::read_csv(file.path(here::here(), "01_Raw_data", "MS_ESV_reads.csv"))
ESV.reads

# Taxonomic assignment results for each ASV
RES.total <- readr::read_csv(file.path(here::here(), "01_Raw_data", "MS_AllAssignments.csv"))
RES.total

# Minor tidy changes

RES.total <- RES.total %>% filter(kingdom %in% c("Animalia", "Metazoa"),
                                  Taxon %nin% c("Homo sapiens", "Polyommatus nivescens", "Cupido alcetas", "Unassigned"),
                                  class %nin% c("Arachnida", "Insecta")) %>% 
                          mutate(assigner = paste(method, threshold),
         method.long = ifelse(method == "LCA", "BLAST - LCA",
                          ifelse(method == "Top hit", "BLAST - Top Hit", "Assigner - IDtaxa")),
         method.threshold = paste(method, threshold),
         Levels = ifelse(str_detect(Taxon, " sp"), "genus", Levels),
         Levels = factor(Levels, levels = c("species", "genus", "family", "order", "class", "phylum", "kingdom")),
         Levels.red = ifelse(Levels %in% c("species", "genus"), Levels %>% as.character(), "> genus"),
         Levels.red = factor(Levels.red, levels = c("species", "genus", "> genus")),
         Taxon = str_remove(Taxon, " sp"),
         Taxon = ifelse(str_ends(Taxon, "dujardini"), str_replace(Taxon, "dujardini", "dujardinii"), Taxon),
        
         Taxo.group = ifelse(phylum %in% c("Sipuncula", "Nemertea", "Brachiopoda", "Annelida"), "Annelida<br>Brachipoda<br>Nemertea",
                             ifelse(phylum %in% c("Porifera", "Cnidaria"), "Cnidaria<br>Porifera", 
                                    phylum)),
         Taxo.group = factor(Taxo.group, levels = c("Rotifera","Cnidaria<br>Porifera", "Bryozoa", "Annelida<br>Brachipoda<br>Nemertea", "Mollusca", "Echinodermata", "Arthropoda", "Chordata" )),
         Phylum = factor(phylum, levels = c("Rotifera","Cnidaria", "Porifera", "Bryozoa", "Annelida", "Brachipoda", "Nemertea", "Sipuncula", "Mollusca", "Arthropoda", "Echinodermata","Chordata" ))
                          )
# COI ref db --------------------------------------------------------------

# Load ranking category!
Taxa.Seq.final <- readr::read_csv(file.path(here::here(), "01_Raw_data", "GSL_included_taxa.csv"))

Taxa.Seq.final

#View(Taxa.Seq.final )

# Rename ranking 
Validity.df <- data.frame(Validity =  c("Sufficient data", "BIN sharing", "Insufficient data", "No genetic data"),
                          Validity.new = c("Reliable",  "Unreliable - BIN sharing", "Unreliable - gaps","No sequences available"))

# Minor tidy changes
Taxa.Seq.final <- Taxa.Seq.final %>% 
  left_join(Validity.df) %>% 
  mutate(
  Taxo.group = ifelse(phylum %in% c("Sipuncula", "Nemertea", "Brachiopoda", "Annelida"), "Annelida<br>Brachipoda<br>Nemertea",
                      ifelse(phylum %in% c("Porifera", "Cnidaria"), "Cnidaria<br>Porifera", 
                             phylum)),
  Taxo.group = factor(Taxo.group, levels = c("Rotifera","Cnidaria<br>Porifera", "Bryozoa", "Annelida<br>Brachipoda<br>Nemertea", "Mollusca", "Echinodermata", "Arthropoda", "Chordata" )),
  
  Validity = factor(Validity, levels = c("Sufficient data", "BIN sharing", "Insufficient data", "No genetic data")),
  Validity.new = factor(Validity.new, levels = c("Reliable",  "Unreliable - BIN sharing","Unreliable - gaps", "No sequences available")),
  Phylum = factor(phylum, levels = c("Rotifera","Cnidaria", "Porifera", "Bryozoa", "Annelida", "Brachiopoda", "Sipuncula", "Nemertea", "Mollusca", "Echinodermata", "Arthropoda", "Chordata" ))
)



# Fig 2. GSL-rl Overview --------------------------------------------------


gg.overview.GSLrl <- Taxa.Seq.final %>% filter(WithinNWA == "Yes", Level == "Species") %>% 
                         group_by(Phylum, Validity.new) %>% 
  summarise(N = n()) %>% 
  mutate(SUM = sum(N),
       freq = N / sum(N)) %>% 
  ggplot(aes(x = 1, y = freq, fill = Validity.new)) +
  geom_bar(width = , stat = "identity", color = "gray10", cex = 0.2) +
  coord_polar("y", start = 0) + 
  scale_fill_manual(name = "Species detection category", 
                    values = c(colorBlindBlack8[4], colorBlindBlack8[5], colorBlindBlack8[2],"gray"))+
  
  geom_text(aes(y = 0.1, label = paste0("n=",SUM)), vjust = 4, col = "black", cex = 3) +
  facet_wrap(~Phylum, nrow = 2) + theme_void() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10, unit = "pt"))

gg.overview.GSLrl

#ggsave(filename = "03_Results/fig_GSLrl_coverage.png",
#        plot = gg.overview.GSLrl,
#       width = 7, height = 3, units = "in")


# Real data - SP level Fig 4 ----------------------------------------------------------------

Clef.taxa <- Taxa.Seq.final %>% select(Taxon = "Name", Validity.new) %>% mutate(Validity.new = as.character(Validity.new))

# NEW FIG
SP.order <- RES.total %>% arrange(phylum, class, order, family, genus, Taxon) %>% pull(Taxon) %>% unique()


# FIg 4A

RES.total <- RES.total %>% 
left_join(Clef.taxa) %>% 
  mutate(Validity.new =  ifelse(Taxon %in% c("Leptodiaptomus minutus", "Caprella unica", "Alcyonidium mamillatum",  "Parasmittina jeffreysi", "Catablema multicirratum", "Serripes laperousii", "Halichondria coerulea", "Polyarthra dolichoptera"), "Unlikely",
                            ifelse(DB == "NCBI_nt", "Likely",
                                   Validity.new)),
         Validity.new = factor(Validity.new, levels = c("Reliable",  "Unreliable - BIN sharing","Unreliable - gaps",   "Likely", "Unlikely" )),
         Taxon = factor(Taxon, levels = SP.order))

figSPa <- RES.total %>%  filter(Levels == "species") %>% 
  left_join(ESV.reads %>% select(QueryAccVer = ID, Nreads.tot)) %>% #View()
  group_by(Phylum, Taxon,method.threshold, Validity.new, DB) %>% summarise(Ntot = sum(Nreads.tot)) %>% 
  filter(#met.com %in% c("IDtaxa 60", "LCA 97", "Top hit 97"),
    Ntot > 0) %>% 
  mutate(Taxon = factor(Taxon, levels = SP.order),
         DB = DB %>% str_replace("_", "-")) %>%
  ggplot(aes(x = method.threshold, y = Taxon, fill = Validity.new)) + 
  geom_bin2d(col = "gray") +
  labs(x="", y="") +
  #  scale_fill_manual(name = "Species detection category", 
  #                    values = c("cornflowerblue", "darkslategray1", "darkgoldenrod1" , "deepskyblue1", "brown1"))+
  #scale_fill_manual(name = "Species detection category", 
  #                  values = c("cornflowerblue", "darkorchid", "darkgoldenrod1" , "cornflowerblue", "brown1"))+
  scale_fill_manual(name = "Species detection category", 
                    values = c(colorBlindBlack8[4], colorBlindBlack8[5], colorBlindBlack8[2], colorBlindBlack8[3], colorBlindBlack8[7]))+
  scale_y_discrete(position = "left") +
  guides(fill=guide_legend(ncol=2,byrow=F)) +
  theme_bw()+
  facet_grid(Phylum ~ DB, scale = "free", space = "free") + #, switch = "y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        axis.text.y = element_text(face = "italic", size = 8, vjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0, size = 8), 
        panel.spacing.y=unit(0.05, "lines"),
        legend.position = "left",
        legend.text = element_text(size = 8),
        panel.grid = element_blank()
  ) #+ coord_flip()

figSPa


# Fig 4 b
figSPb <- RES.total %>%   filter(Levels == "species") %>% 
  left_join(ESV.reads %>% select(QueryAccVer = ID, Nreads.tot)) %>% #View()  group_by(Phylum, Taxon,met.com, Validity, DB) %>% 
  group_by(method.threshold, Validity.new, DB, Taxon) %>%  
  summarise(Ntot = sum(Nreads.tot)) %>% 
  filter(Ntot > 0) %>% 
   group_by(method.threshold, Validity.new, DB) %>% 
  summarise(Ndetect = length(unique(Taxon))) %>% 
  mutate(DB = DB %>% str_replace("_", "-")) %>% 
  ggplot(aes(x = method.threshold, y = Ndetect, fill = Validity.new)) + 
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  #geom_hline(yintercept = c(68, 72)) +
  labs(y="N species detected", x="") +
  #  scale_fill_manual(name = "Species detection category", 
  #                    values = c("cornflowerblue", "darkslategray1", "darkgoldenrod1" , "deepskyblue1", "brown1"))+
  #scale_fill_manual(name = "Species detection category", 
  #                  values = c("cornflowerblue", "darkorchid", "darkgoldenrod1" , "cornflowerblue", "brown1"))+
  scale_fill_manual(name = "Species rank", 
                    values = c(colorBlindBlack8[4], colorBlindBlack8[5], colorBlindBlack8[2], colorBlindBlack8[3], colorBlindBlack8[7]))+
  
  guides(fill=guide_legend(ncol=2,byrow=F)) +
  theme_bw()+
  facet_grid(. ~ DB, scale = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        axis.text.y = element_text( size = 8, vjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0, size = 8), 
        panel.spacing.y=unit(0.1, "lines"),
        legend.position = "right",
        legend.text = element_text(size = 8)#,
        #panel.grid = element_blank()
  ) #+ coord_flip()
figSPb


figSPc <- RES.total %>% filter(Levels == "species") %>% 
  left_join(ESV.reads %>% select(QueryAccVer = ID, Nreads.tot)) %>% #View()  group_by(Phylum, Taxon, Validity, DB) %>% 
  
 # filter(assigner %in% c("IDtaxa 40", "Top hit 95"),
  #       !(method.threshold == "Top hit 95" & DB == "GSL_rl")) %>% 
  group_by(Phylum, Taxon, Validity.new, DB) %>% 
  summarise(Ntot = sum(Nreads.tot)) %>% 
  filter(Ntot > 0) %>% 
  mutate(Taxon = factor(Taxon, levels = SP.order),
         Validity.new = as.character(Validity.new)) %>%
  select(-Ntot) %>% 
  pivot_wider(names_from = DB, values_from = Validity.new, values_fill = "No assigments") %>%  
  group_by(`NCBI_nt`,`GSL_rl`) %>% 
  summarise(N = n()) %>% 
  mutate(`NCBI_nt` = factor(`NCBI_nt`, levels = c("Likely", "Unlikely", "No assigments" )),
         `GSL_rl` = factor(`GSL_rl`, levels = c("Reliable",  "Unreliable - BIN sharing","Unreliable - gaps", "No assigments" )))  %>% 
  complete(`NCBI_nt`,`GSL_rl`) %>% 
  
  ggplot(aes(x = `NCBI_nt`, y = `GSL_rl`, fill = N)) + 
  geom_bin2d(col = "black")  +
  geom_text(aes(label = N), size = 3) +
  scale_fill_gradient(low = "gray90", high = "gray50", na.value = "white") +
  labs(x ="NCBI-nt/Top hit 95", y = "GSL-rl/IDtaxa 40") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        axis.text.y = element_text( size = 8, vjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0, size = 8), 
        panel.spacing.y=unit(0.1, "lines"),
        legend.position = "none",
        legend.text = element_text(size = 8)#,
        #panel.grid = element_blank()
  ) #+ coord_flip()
figSPc


figSPc.v2 <- RES.total %>% filter(Levels == "species") %>% 
  left_join(ESV.reads %>% select(QueryAccVer = ID, Nreads.tot)) %>% #View()  group_by(Phylum, Taxon, Validity, DB) %>% 
  
  # filter(assigner %in% c("IDtaxa 40", "Top hit 95"),
  #       !(method.threshold == "Top hit 95" & DB == "GSL_rl")) %>% 
  group_by(Phylum, Taxon, Validity.new, DB) %>% 
  summarise(Ntot = sum(Nreads.tot)) %>% 
  filter(Ntot > 0) %>% 
  mutate(Taxon = factor(Taxon, levels = SP.order),
         Validity.new = as.character(Validity.new)) %>%
  select(-Ntot) %>% 
  pivot_wider(names_from = DB, values_from = Validity.new, values_fill = "No assigments") %>%  
  group_by(`NCBI_nt`,`GSL_rl`) %>% 
  summarise(N = n()) %>% 
  mutate(`NCBI_nt` = factor(`NCBI_nt`, levels = c("Likely", "Unlikely", "No assigments" )),
         `GSL_rl` = factor(`GSL_rl`, levels = c("Reliable",  "Unreliable - BIN sharing","Unreliable - gaps", "No assigments" )))  %>% 
  complete(`NCBI_nt`,`GSL_rl`) %>% 
  
  ggplot(aes(x = `NCBI_nt`, y = `GSL_rl`, fill = N)) + 
  geom_bin2d(col = "black")  +
  geom_text(aes(label = N), size = 3) +
  scale_fill_gradient(low = "gray90", high = "gray50", na.value = "white") +
  scale_y_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 15))+
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 15))+
  labs(x ="NCBI-nt", y = "GSL-rl") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        axis.text.y = element_text( size = 8, vjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0, size = 8), 
        panel.spacing.y=unit(0.1, "lines"),
        legend.position = "none",
        legend.text = element_text(size = 8)#,
        #panel.grid = element_blank()
  ) #+ coord_flip()
figSPc.v2




figSP2 <- ggpubr::ggarrange(figSPa + theme(legend.position = "none"), 
                            ggarrange(figSPb + theme(legend.position = "none"),
                                     
                                     figSPc,
                                     get_legend(figSPb),
                                     labels = c("B", "C", ""),
                                     nrow =3, heights = c(3,3,2)), 
                              
                           # bar plot spaning two columns                       # box plot and scatter plot
                           labels = c("A", ""),
                           ncol = 2, nrow = 1, widths = c(5,3), 
                           common.legend = F)

figSP2

# EN changeant 3c
figSP2.v2 <- ggpubr::ggarrange(figSPa + theme(legend.position = "none"), 
                            ggarrange(figSPb + theme(legend.position = "none"),
                                      
                                      figSPc.v2,
                                      get_legend(figSPb),
                                      labels = c("B", "C", ""),
                                      nrow =3, heights = c(3,3,2)), 
                            
                            # bar plot spaning two columns                       # box plot and scatter plot
                            labels = c("A", ""),
                            ncol = 2, nrow = 1, widths = c(5,3), 
                            common.legend = F)

figSP2.v2

#ggsave(filename = file.path(here::here(), "03_Results", "fig_Assignement_SP_v2.png"), 
#       plot = figSP2.v2, width = 8, height = 9, units = "in", bg = "white")



# Stats - Number of unique species
RES.total %>%   filter( Levels == "species") %>%# head()
  left_join(ESV.reads %>% select(QueryAccVer = ID, Nreads.tot)) %>% 
  group_by(Phylum, Taxon,method.threshold, method, Validity.new, DB) %>% summarise(Ntot = sum(Nreads.tot)) %>% 
  filter(Ntot > 0) %>% 
  mutate(Taxon = factor(Taxon, levels = SP.order)) %>%
  group_by(DB) %>% summarise(N = length(unique(Taxon)))


RES.total %>%   filter(Levels == "species") %>% 
  left_join(ESV.reads %>% select(QueryAccVer = ID, Nreads.tot)) %>% #View()  group_by(Phylum, Taxon,met.com, Validity, DB) %>% 
  group_by(method.threshold, Validity.new, DB, Taxon) %>%  
  summarise(Ntot = sum(Nreads.tot)) %>% 
  filter(Ntot > 0) %>% 
  group_by(DB, method.threshold, Validity.new) %>% 
  summarise(Ndetect = length(unique(Taxon))) %>% 
  mutate(Ntotal = sum(Ndetect),
         freq = Ndetect / Ntotal) %>% 
  filter(Validity.new == "Unlikely")

# NCBI Blast over GSL - Dataset Fig 3---------------------------------------------------------

#load(file.path("01_Raw_data/MS_NCBI.test.assignments.Rdata"))
#
#head(RES.postONncbi)
#
#write_csv(RES.postONncbi, file.path(here::here(), "01_Raw_data", "MS_NCBI.TestAssignments.csv"))

# Load assignment
RES.NCBI.test <- read_csv(file.path(here::here(), "01_Raw_data", "MS_NCBI.TestAssignments.csv"))


# Tidy changes
RES.NCBI.test <-RES.NCBI.test %>% filter(str_ends(Taxon.ori, " sp", negate = T)) %>% 
  mutate(Taxo.group = ifelse(Phylum.ori %in% c("Sipuncula", "Nemertea", "Brachiopoda", "Annelida"), "Annelida<br>Brachiopoda<br>Nemertea",
                             ifelse(Phylum.ori %in% c("Porifera", "Cnidaria"), "Cnidaria<br>Porifera", Phylum.ori)),
         Taxo.group = factor(Taxo.group, levels = c("Cnidaria<br>Porifera", "Annelida<br>Brachiopoda<br>Nemertea", "Mollusca", "Arthropoda",  "Echinodermata","Chordata" )),
         Levels.group = ifelse(Levels %in% c("order", "phylum", "family", "class", "kingdom", "Unassigned"), "> genus or unassigned", Levels),
         Levels.group = factor(Levels.group, levels = c("species", "genus", "> genus or unassigned")),    
         
         Levels.group2 = ifelse(Levels %in% c("order", "phylum", "family", "class", "kingdom"), "> genus", Levels),
         Levels.group2 = factor(Levels.group2, levels = c("species", "genus", "> genus", "Unassigned")),    
         
         #Similar = ifelse(Taxon.ori == Taxon & Levels == "species", "Right species identification",
          #                ifelse(Levels == "species", "Wrong identification at species levels",
          #                       ifelse(Genus.ori == genus & Levels == "genus", "Right genus identification",
          #                              ifelse(Levels == "genus", "Wrong identification at genus levels",
          #                                     #ifelse(Order.ori == order, "Order",
          #                                     "Assignation at higher taxonomic level or unassigned"
          #                              )))),
         #Similar = ifelse(is.na(Similar), "Unassigned", Similar),
        # 
        # Similar = factor(Similar, levels = c("Right species identification",
        #                                      "Right genus identification",
        #                                      "Wrong identification at species levels",
        #                                      "Wrong identification at genus levels",
        #                                      "Assignation at higher taxonomic level or unassigned"
        # ))
        Similar = ifelse(Taxon.ori == Taxon & Levels == "species", "Right identification",
                         ifelse(Levels == "species", "Wrong identification",
                                ifelse(Genus.ori == genus & Levels == "genus", "Right identification",
                                       ifelse(Levels == "genus", "Wrong identification",
                                              #ifelse(Order.ori == order, "Order",
                                              "PROBLEMS"
                                       )))),
        method.graph = ifelse(method == "LCA", "<br>LCA<br>", method)
        ) 



# Interesting stat
RES.NCBI.test %>%  group_by(method, threshold, Levels) %>% summarise(N = n()) %>% 
                   pivot_wider(names_from = Levels, values_from = N)



# Stats for the article

RES.NCBI.test %>%   group_by(method, threshold, Levels.group) %>% summarise(N = n()) %>%
  mutate(freq = N / sum(N)) %>% 
  filter(Levels.group == "species") %>% 
  group_by(method, threshold) %>% 
  summarise(min = min(freq),
            max = max(freq))


RES.NCBI.test %>% group_by(method, threshold, Levels.group) %>% summarise(N = n()) %>%
  mutate(freq = N / sum(N)) %>% 
  filter(Levels.group == "genus") %>% 
  group_by(method, threshold, Levels.group) %>% 
  summarise(mean = mean(freq))



RES.NCBI.test %>% group_by(method, threshold, Taxo.group, Levels.group) %>% summarise(N = n()) %>%
  mutate(freq = N / sum(N)) %>% 
  filter(Levels.group == "species", threshold == 97) %>% 
  arrange(method, freq)

# Range of accuracy (fig3A)
RES.NCBI.test  %>%  group_by(method, threshold, Levels.group, Similar) %>% summarise(N = n()) %>%
  mutate(freq = N / sum(N)) %>%
  filter(Similar == "Right identification", Levels.group %in% c("species","genus")) %>% 
  group_by(Levels.group, method) %>% 
  summarise(min = min(freq),
            max = max(freq))

# Range of accuracy (fig3B)
RES.NCBI.test  %>% filter( Levels.group == "species", threshold == 97) %>% 
  group_by(method, threshold, Levels.group,  Taxo.group, Similar) %>% summarise(N = n()) %>%
  mutate(freq = N / sum(N)) %>%
  filter(Similar == "Right identification", Levels.group %in% c("species","genus"), threshold == 97) %>% 
  group_by(Levels.group, method) %>% 
  summarise(min = min(freq),
            max = max(freq))

RES.NCBI.test %>% group_by(method, threshold,  Levels.group, Similar) %>% summarise(N = n()) %>%
  mutate(freq = N / sum(N)) %>%
  #filter(Similar == "Wrong identification", Levels.group == "species", threshold == 97) %>% 
  arrange(method, freq) #%>% View()

#RES.NCBI.test %>% pull(QueryAccVer) %>% unique() %>% length()


# NEW VERSION

grapha.acc.over <- RES.NCBI.test %>% 
  group_by(method.graph, threshold, Levels.group) %>%  summarise(Nok = length(SeqName[Similar == "Right identification"]),
                                                                 Nnotok = length(SeqName[Similar != "Right identification"]),
                                                                 N = n()) %>%
  mutate(freqOK = Nok / sum(N),
         freqnotOK = Nnotok / sum(N),
         threshold = factor(threshold),
         Levels.group = factor(Levels.group, levels =  c("> genus or unassigned", "genus", "species"))) %>% #View()
  filter(Levels.group %in% c("species")) %>%
  pivot_longer(cols = c(freqOK, freqnotOK), values_to = "freq") %>% 
  # filter(Similar == "Right identification") %>% 
  ggplot(aes(y=freq, x = threshold, fill = name)) +
  #geom_point(position=position_dodge(width=0.3), size = 3) +
  scale_y_continuous(limits = c(0,1)) +
  geom_bar(stat= "identity", col = "black")+
  #geom_bar(stat= "identity",col = "darkgray", position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values = c("gray60", "gray95"), labels = c("Inaccurate", "Accurate") )+
  
  facet_grid(. ~method.graph , scale = "free") + 
  xlab("") +
  ylab(expression(paste("Prop. of assigments"))) +  theme_bw() + 
  theme_bw() + 
  
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom",
        strip.text = element_markdown(),
        legend.title = element_blank())

grapha.acc.over


grapha.acc.sp <- RES.NCBI.test %>% filter(threshold == 97) %>% 
  group_by(method, threshold, Taxo.group, Levels.group) %>% summarise(Nok = length(SeqName[Similar == "Right identification"]),
                                                                      Nnotok = length(SeqName[Similar != "Right identification"]),
                                                                      N = n()) %>%
  mutate(freqOK = Nok / sum(N),
         freqnotOK = Nnotok / sum(N),
         Levels.group = factor(Levels.group, levels =  c("> genus or unassigned", "genus", "species"))) %>% #View()
  filter(Levels.group %in% c("species")) %>%
  pivot_longer(cols = c(freqOK, freqnotOK), values_to = "freq") %>% 
  # filter(Similar == "Right identification") %>% 
  ggplot(aes(y=freq, x = method, fill = name)) +
  #geom_point(position=position_dodge(width=0.3), size = 3) +
  scale_y_continuous(limits = c(0,1)) +
  geom_bar(stat= "identity", col = "black")+
  #geom_bar(stat= "identity",col = "darkgray", position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values = c("gray60", "gray95"), labels = c("Inaccurate", "Accurate") )+
  facet_grid(. ~Taxo.group , scale = "free") + 
  xlab("") +
  ylab(expression(paste("Prop. of assigments"))) +  theme_bw() + 
  theme_bw() + 
  
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom",
        strip.text = element_markdown(),
        legend.title = element_blank())

grapha.acc.sp

# 
# graphb.over <- RES.NCBI.test %>% filter(Levels %in% c("species", "genus")) %>% 
#   group_by(method.graph, threshold, Levels.group, Similar) %>% summarise(N = n()) %>%
#   mutate(freq = N / sum(N),
#          threshold = factor(threshold)) %>% 
#   # Levels.group = factor(Levels.group, levels =  c("> genus or unassigned", "genus", "species"))) %>% #View()
#   filter(Similar != "Wrong identification",
#          Levels.group %in% c("species")) %>% 
#   ggplot(aes(y=freq, x = threshold, fill = Levels.group, group = Levels.group)) +
#   #geom_point(position=position_dodge(width=0.3), size = 3) +
#   #scale_y_continuous(limits = c(,1)) +
#   #geom_point() +
#   ungeviz::geom_hpline()+
#   #geom_bar(stat= "identity", position = "dodge") +
#   scale_y_continuous(limits = c(0.87,1)) +
#    facet_grid(. ~ method.graph , scale = "free") + 
#   xlab("") +
#   ylab(expression(paste("Accuracy"))) +  theme_bw() + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#         legend.position = "bottom",
#         strip.text = element_markdown(),
#         legend.title = element_blank())
# 
# graphb.over

# 


fig2.v1 <- ggarrange(grapha.acc.over, grapha.acc.sp,
          labels = c("A", "B"),
          nrow = 1, ncol = 2, widths = c(2,5), align = "hv",
          common.legend = TRUE, legend = "bottom"
)

fig2.v1



#ggsave(filename = file.path(here::here(), "03_Results", "fig_NCBI.png"), 
#       plot = fig2.v1, width = 9, height = 4, units = "in", bg = "white")


# Stats

 RES.NCBI.test %>% filter(threshold == 99) %>% 
  group_by(method, threshold, Taxo.group, Levels.group) %>% summarise(Nok = length(SeqName[Similar == "Right identification"]),
                                                                      Nnotok = length(SeqName[Similar != "Right identification"]),
                                                                      N = n()) %>%
  mutate(freqOK = Nok / sum(N),
         freqnotOK = Nnotok / sum(N),
         Levels.group = factor(Levels.group, levels =  c("> genus or unassigned", "genus", "species"))) %>% #View()
  filter(Levels.group %in% c("species")) %>%
  pivot_longer(cols = c(freqOK, freqnotOK), values_to = "freq") %>% 
  filter(name == "freqOK")
 
 