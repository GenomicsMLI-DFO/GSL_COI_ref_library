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

ESV.reads <- readr::read_csv(file.path(here::here(), "01_Raw_data", "MS_ESV_reads.csv"))

RES.total <- readr::read_csv(file.path(here::here(), "01_Raw_data", "MS_AllAssignments.csv"))
RES.total

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



# Annotation

Taxa.Seq.final <- readr::read_csv2(file.path(here::here(), "00_GSL_library", "GSL_included_taxa.csv"))

Taxa.Seq.final


Clef.taxa <- Taxa.Seq.final %>% select(Taxon = "Name", Validity)

# NEW FIG
SP.order <- RES.total %>% arrange(phylum, class, order, family, genus, Taxon) %>% pull(Taxon) %>% unique()



RES.total <- RES.total %>% 
left_join(Clef.taxa) %>% 
  mutate(Validity =  ifelse(Taxon %in% c("Leptodiaptomus minutus", "Caprella unica", "Alcyonidium mamillatum",  "Parasmittina jeffreysi", "Catablema multicirratum", "Serripes laperousii", "Halichondria coerulea", "Polyarthra dolichoptera"), "Unlikely",
                            ifelse(DB == "NCBI_nt", "Likely",
                                   Validity)),
         Validity = factor(Validity, levels = c("Sufficient data", "BIN sharing","Insufficient data", "Likely", "Unlikely" )),
         Taxon = factor(Taxon, levels = SP.order))


# Real data - ESV level  -------------------------------------------------------------

# The number of ESV assigned to a taxon, for each method
RES.total %>% group_by(DB, method, threshold) %>% 
              summarise(N=n())


RES.total %>% group_by(DB, method, threshold) %>% 
              summarise(N=n()) %>% ungroup() %>% summarise(mean = mean(N),
                                               sd = sd(N) )

RES.total %>% group_by(Levels) %>% summarise(N = n())
# Figure 2


figESVa <- RES.total %>%  group_by(DB, method.long, assigner, Levels.red) %>% 
  summarise(N=n()) %>% #View()
  filter(method.long == "BLAST - LCA") %>% 
  ggplot(aes(y=N, x = assigner, fill = Levels.red)) +
  geom_bar(stat = "identity", col = "darkgray", position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = c("gray15", "gray60", "gray95"))+
  scale_y_continuous(limits = c(0,165)) +
  facet_grid(. ~ DB, space = "free", scale = "free") + 
  labs(x = "", y = "N ESV") + 
  theme_bw() + theme(legend.position = "bottom",
                     axis.text.x =  element_text(angle = 90, vjust = 0.5), 
                     legend.title = element_blank())

figESVa

figESVb <- RES.total %>% 
  group_by(DB, method.long, assigner, Levels.red) %>% 
  summarise(N=n()) %>% #View()
  filter(method.long == "BLAST - Top Hit") %>% 
  ggplot(aes(y=N, x = assigner, fill = Levels.red)) +
  geom_bar(stat = "identity", col = "darkgray", position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = c("gray15", "gray60", "gray95"))+
  scale_y_continuous(limits = c(0,165)) +
  facet_grid(. ~ DB, space = "free", scale = "free") + 
  labs(x = "", y = "N ESV") + 
  theme_bw() + theme(legend.position = "bottom",
                     axis.text.x =  element_text(angle = 90, vjust = 0.5), 
                     legend.title = element_blank())

figESVb

figESVc <- RES.total %>% 
  group_by(DB, method.long, assigner, Levels.red) %>% 
  summarise(N=n()) %>% #View()
  filter(method.long == "Assigner - IDtaxa") %>% 
  ggplot(aes(y=N, x = assigner, fill = Levels.red)) +
  geom_bar(stat = "identity", col = "darkgray", position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = c("gray15", "gray60", "gray95"))+
  scale_y_continuous(limits = c(0,165)) +
  facet_grid(. ~ DB, space = "free", scale = "free") + 
  labs(x = "", y = "N ESV") + 
  theme_bw() + theme(legend.position = "bottom",
                     axis.text.x =  element_text(angle = 90, vjust = 0.5), 
                     legend.title = element_blank())

figESVc

figESV <- ggpubr::ggarrange(figESVa, figESVb, figESVc, 
          labels = c("A", "B", "C"),
          nrow = 1, ncol = 3, widths = c(5,5,3), align = "hv",
          common.legend = TRUE, legend = "bottom"
          )

figESV

#ggsave(filename = file.path(here::here(), "03_Results", "fig_Assignement_ESV.png"), plot = figESV, width = 6, height = 3, units = "in", bg = "white")



# Real data - SP level ----------------------------------------------------------------


figSPa <- RES.total %>%  filter(Levels == "species") %>% 
  left_join(ESV.reads %>% select(QueryAccVer = ID, Nreads.tot)) %>% #View()
  group_by(Phylum, Taxon,method.threshold, Validity, DB) %>% summarise(Ntot = sum(Nreads.tot)) %>% 
  filter(#met.com %in% c("IDtaxa 60", "LCA 97", "Top hit 97"),
    Ntot > 0) %>% 
  mutate(Taxon = factor(Taxon, levels = SP.order)) %>%
  ggplot(aes(x = method.threshold, y = Taxon, fill = Validity)) + 
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



figSPb <- RES.total %>%   filter(Levels == "species") %>% 
  left_join(ESV.reads %>% select(QueryAccVer = ID, Nreads.tot)) %>% #View()  group_by(Phylum, Taxon,met.com, Validity, DB) %>% 
  group_by(method.threshold, Validity, DB, Taxon) %>%  
  summarise(Ntot = sum(Nreads.tot)) %>% 
  filter(Ntot > 0) %>% 
   group_by(method.threshold, Validity, DB) %>% 
  summarise(Ndetect = length(unique(Taxon))) %>% 
  
  ggplot(aes(x = method.threshold, y = Ndetect, fill = Validity)) + 
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  #geom_hline(yintercept = c(68, 72)) +
  labs(y="N detections", x="") +
  #  scale_fill_manual(name = "Species detection category", 
  #                    values = c("cornflowerblue", "darkslategray1", "darkgoldenrod1" , "deepskyblue1", "brown1"))+
  #scale_fill_manual(name = "Species detection category", 
  #                  values = c("cornflowerblue", "darkorchid", "darkgoldenrod1" , "cornflowerblue", "brown1"))+
  scale_fill_manual(name = "Species detection category", 
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
  group_by(Phylum, Taxon, Validity, DB) %>% 
  summarise(Ntot = sum(Nreads.tot)) %>% 
  filter(Ntot > 0) %>% 
  mutate(Taxon = factor(Taxon, levels = SP.order),
         Validity = as.character(Validity)) %>%
  select(-Ntot) %>% 
  pivot_wider(names_from = DB, values_from = Validity, values_fill = "No assigments") %>%  
  group_by(`NCBI_nt`,`GSL_rl`) %>% 
  summarise(N = n()) %>% 
  mutate(`NCBI_nt` = factor(`NCBI_nt`, levels = c("Likely", "Unlikely", "No assigments" )),
         `GSL_rl` = factor(`GSL_rl`, levels = c("Sufficient data", "BIN sharing","Insufficient data", "No assigments" )))  %>% 
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


figSP <- ggpubr::ggarrange(ggarrange(figSPb + theme(legend.position = "none"),
                    
                    figSPc,
                    get_legend(figSPb),
                    labels = c("A", "B", ""),
                    nrow =3, heights = c(3,3,2)), 
          figSPa + theme(legend.position = "none"),    
          # bar plot spaning two columns                       # box plot and scatter plot
          labels = c("", "C"),
          ncol = 2, nrow = 1, widths = c(3,5), 
          common.legend = F)

figSP

ggsave(filename = file.path(here::here(), "03_Results", "fig_Assignement_SP.png"), 
       plot = figSP, width = 8, height = 9, units = "in", bg = "white")


# Stats
RES.total %>%   filter( Levels == "species") %>%# head()
  left_join(ESV.reads %>% select(QueryAccVer = ID, Nreads.tot)) %>% 
  group_by(Phylum, Taxon,method.threshold, method, Validity, DB) %>% summarise(Ntot = sum(Nreads.tot)) %>% 
  filter(Ntot > 0) %>% 
  mutate(Taxon = factor(Taxon, levels = SP.order)) %>%
  group_by(DB) %>% summarise(N = length(unique(Taxon)))



# NCBI Blast over GSL - Dataset ---------------------------------------------------------

#load(file.path("01_Raw_data/MS_NCBI.test.assignments.Rdata"))
#
#head(RES.postONncbi)
#
#write_csv(RES.postONncbi, file.path(here::here(), "01_Raw_data", "MS_NCBI.TestAssignments.csv"))

RES.NCBI.test <- read_csv(file.path(here::here(), "01_Raw_data", "MS_NCBI.TestAssignments.csv"))



RES.NCBI.test <-RES.NCBI.test %>% filter(str_ends(Taxon.ori, " sp", negate = T)) %>% 
  mutate(Taxo.group = ifelse(Phylum.ori %in% c("Sipuncula", "Nemertea", "Brachiopoda", "Annelida"), "Annelida<br>Brachipoda<br>Nemertea",
                             ifelse(Phylum.ori %in% c("Porifera", "Cnidaria"), "Cnidaria<br>Porifera", Phylum.ori)),
         Taxo.group = factor(Taxo.group, levels = c("Cnidaria<br>Porifera", "Annelida<br>Brachipoda<br>Nemertea", "Mollusca", "Arthropoda",  "Echinodermata","Chordata" )),
         Levels.group = ifelse(Levels %in% c("order", "phylum", "family", "class", "kingdom", "Unassigned"), "> genus or unassigned", Levels),
         Levels.group = factor(Levels.group, levels = c("species", "genus", "> genus or unassigned")),    
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

RES.NCBI.test  %>%  group_by(method, threshold, Levels.group, Similar) %>% summarise(N = n()) %>%
  mutate(freq = N / sum(N)) %>%
  filter(Similar == "Wrong identification", Levels.group == "genus") %>% 
  group_by(method) %>% 
  summarise(min = min(freq),
            max = max(freq))

RES.NCBI.test %>% group_by(method, threshold,  Levels.group, Similar) %>% summarise(N = n()) %>%
  mutate(freq = N / sum(N)) %>%
  #filter(Similar == "Wrong identification", Levels.group == "species", threshold == 97) %>% 
  arrange(method, freq) %>% View()

RES.NCBI.test %>% pull(QueryAccVer) %>% unique() %>% length()

RES.NCBI.test 

grapha <-RES.NCBI.test %>% filter(threshold == 97) %>% 
  group_by(method, threshold, Taxo.group, Levels.group) %>% summarise(N = n()) %>%
  mutate(freq = N / sum(N),
         Levels.group = factor(Levels.group, levels =  c("> genus or unassigned", "genus", "species"))) %>% #View()
  #filter(Levels %in% c("species", "genus")) %>% 
  ggplot(aes(y=freq, x = method, fill = Levels.group)) +
  #geom_point(position=position_dodge(width=0.3), size = 3) +
  scale_y_continuous(limits = c(0,1)) +
  geom_bar(stat= "identity",col = "darkgray") +
  scale_fill_manual(values = c("gray15", "gray60", "gray95"), limits = c("species","genus",  "> genus or unassigned"))+
  
  facet_grid(. ~Taxo.group , scale = "free") + 
  xlab("") +
  ylab(expression(paste("Proportion of assigments"))) +  theme_bw() + 
  theme_bw() + 
  
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom",
        strip.text = element_markdown(),
        legend.title = element_blank())

grapha


grapha.acc <- RES.NCBI.test %>% filter(threshold == 97) %>% 
  group_by(method, threshold, Taxo.group, Levels.group) %>% summarise(Nok = length(SeqName[Similar == "Right identification"]),
                                                                      N = n()) %>%
  mutate(freq = Nok / sum(N),
         Levels.group = factor(Levels.group, levels =  c("> genus or unassigned", "genus", "species"))) %>% #View()
  #filter(Levels %in% c("species", "genus")) %>% 
 # filter(Similar == "Right identification") %>% 
  ggplot(aes(y=freq, x = method, fill = Levels.group)) +
  #geom_point(position=position_dodge(width=0.3), size = 3) +
  scale_y_continuous(limits = c(0,1)) +
  geom_bar(stat= "identity", col = "darkgray")+
  #geom_bar(stat= "identity",col = "darkgray", position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values = c("gray15", "gray60", "gray95"), limits = c("species","genus",  "> genus or unassigned"))+
  facet_grid(. ~Taxo.group , scale = "free") + 
  xlab("") +
  ylab(expression(paste("Proportion of accurate assigments"))) +  theme_bw() + 
  theme_bw() + 
  
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom",
        strip.text = element_markdown(),
        legend.title = element_blank())

grapha.acc

# 
graphb <-RES.NCBI.test %>% filter(Levels %in% c("species", "genus"),
                                  threshold == 97) %>% 
  group_by(method, threshold, Taxo.group, Levels.group, Similar) %>% summarise(N = n()) %>%
  mutate(freq = N / sum(N)) %>% 
         #Levels.group = factor(Levels.group, levels =  c("> genus or unassigned", "genus", "species"))) %>% #View()
  filter(Similar != "Wrong identification") %>% 
  ggplot(aes(y=freq, x = method, fill = Levels.group, group = Levels.group)) +
  #geom_point(position=position_dodge(width=0.3), size = 3) +
  #scale_y_continuous(limits = c(0,0.15)) +
  geom_bar(stat= "identity", position = "dodge") +
  scale_shape_manual(values = c(21,22))+
  scale_fill_manual(values = c("gray15", "gray60", "gray95"), limits = c("species","genus",  "> genus or unassigned"))+
  
  facet_grid(. ~ Taxo.group , scale = "free") + 
  xlab("") +
  ylab(expression(paste("Accuracy"))) +  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom",
        strip.text = element_markdown(),
        legend.title = element_blank())

graphb 


grapha.over <- RES.NCBI.test %>% 
  group_by(method.graph, threshold, Levels.group) %>% summarise(N = n()) %>%
  mutate(freq = N / sum(N),
         threshold = factor(threshold),
         Levels.group = factor(Levels.group, levels =  c("> genus or unassigned", "genus", "species"))) %>% #View()
  #filter(Levels %in% c("species", "genus")) %>% 
  ggplot(aes(y=freq, x = threshold, fill = Levels.group)) +
  #geom_point(position=position_dodge(width=0.3), size = 3) +
  scale_y_continuous(limits = c(0,1)) +
  geom_bar(stat= "identity",col = "darkgray") +
  scale_fill_manual(values = c("gray15", "gray60", "gray95"), limits = c("species","genus",  "> genus or unassigned"))+
  
  facet_grid(. ~method.graph , scale = "free") + 
  xlab("") +
  ylab(expression(paste("Proportion of assigments"))) +  theme_bw() + 
  theme_bw() + 
  
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom",
        strip.text = element_markdown(),
        legend.title = element_blank())

grapha.over

grapha.acc.over <- RES.NCBI.test %>% 
  group_by(method.graph, threshold, Levels.group) %>% summarise(Nok = length(SeqName[Similar == "Right identification"]),
                                                            N = n()) %>%
  mutate(freq = Nok / sum(N),
         threshold = factor(threshold),
         Levels.group = factor(Levels.group, levels =  c("> genus or unassigned", "genus", "species"))) %>% #View()
  #filter(Levels %in% c("species", "genus")) %>% 
  ggplot(aes(y=freq, x = threshold, fill = Levels.group)) +
  #geom_point(position=position_dodge(width=0.3), size = 3) +
  scale_y_continuous(limits = c(0,1)) +
  geom_bar(stat= "identity",col = "darkgray") +
  scale_fill_manual(values = c("gray15", "gray60", "gray95"), limits = c("species","genus",  "> genus or unassigned"))+
  
  facet_grid(. ~method.graph , scale = "free") + 
  xlab("") +
  ylab(expression(paste("Proportion of accurate assigments"))) +  theme_bw() + 
  theme_bw() + 
  
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom",
        strip.text = element_markdown(),
        legend.title = element_blank())

grapha.acc.over


graphb.over <- RES.NCBI.test %>% filter(Levels %in% c("species", "genus")) %>% 
    group_by(method.graph, threshold, Levels.group, Similar) %>% summarise(N = n()) %>%
  mutate(freq = N / sum(N),
         threshold = factor(threshold)) %>% 
        # Levels.group = factor(Levels.group, levels =  c("> genus or unassigned", "genus", "species"))) %>% #View()
  filter(Similar != "Wrong identification") %>% 
  ggplot(aes(y=freq, x = threshold, fill = Levels.group, group = Levels.group)) +
  #geom_point(position=position_dodge(width=0.3), size = 3) +
  #scale_y_continuous(limits = c(,1)) +
  geom_bar(stat= "identity", position = "dodge") +
  scale_shape_manual(values = c(21,22))+
  scale_fill_manual(values = c("gray15", "gray60", "gray95"))+
  facet_grid(. ~ method.graph , scale = "free") + 
  xlab("") +
  ylab(expression(paste("Accuracy"))) +  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom",
        strip.text = element_markdown(),
        legend.title = element_blank())

graphb.over




graphc.over <- RES.NCBI.test %>% 
  group_by(method.graph, threshold, Levels.group, Similar) %>% summarise(N = n()) %>%
  mutate(freq = N / sum(N),
         threshold = factor(threshold)) %>% #View()
  filter(Similar == "Wrong identification", Levels.group == "genus") %>% 
  ggplot(aes(y=freq, x = threshold, fill = Levels.group, group = Levels.group)) +
  #geom_point(position=position_dodge(width=0.3), size = 3) +
  scale_y_continuous(limits = c(0,0.50)) +
  geom_bar(stat= "identity", position = "dodge") +
  scale_shape_manual(values = c(21,22))+
  scale_fill_manual(values = c("gray60"))+
  facet_grid(. ~ method.graph ) + 
  xlab("") +
  ylab(expression(paste("Error rate"))) +  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none",
        strip.text = element_markdown(),
        legend.title = element_blank())

graphc.over 



graphc <-RES.NCBI.test %>% filter(threshold == 97) %>% 
  group_by(method.graph, threshold, Taxo.group, Levels.group, Similar) %>% summarise(N = n()) %>%
  mutate(freq = N / sum(N)) %>% #View()
  filter(Similar == "Wrong identification",
         Levels.group == "genus") %>% 
  bind_rows(data.frame(method = "LCA", threshold = 97, Taxo.group = "Arthropoda", Levels.group = "genus")) %>% 
  ggplot(aes(y=freq, x = method.graph, fill = Levels.group, group = Levels.group)) +
  #geom_point(position=position_dodge(width=0.3), size = 3) +
  scale_y_continuous(limits = c(0,0.5)) +
  geom_bar(stat= "identity", position = "dodge") +
  scale_shape_manual(values = c(21,22))+
  scale_fill_manual(values = c("gray60"))+
  facet_grid(. ~ Taxo.group) + 
  xlab("") +
  ylab(expression(paste("Error rate"))) +  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none",
        strip.text = element_markdown(),
        legend.title = element_blank())

graphc 


ggarrange(graphc.over, graphc,
          labels = c("A", "B"),
          nrow = 1, ncol = 2, widths = c(2,5), align = "hv",
          common.legend = TRUE, legend = "none"
)


# Fig supp overall

graphd <- RES.postONncbi %>% filter(str_ends(Taxon.ori, " sp", negate = T),
                                    #Levels %in% c("species"),
                                    threshold == 97) %>% 
  mutate(Taxo.group = ifelse(Phylum.ori %in% c("Sipuncula", "Nemertea", "Brachiopoda", "Annelida"), "Annelida<br>Brachipoda<br>Nemertea",
                             ifelse(Phylum.ori %in% c("Porifera", "Cnidaria"), "Cnidaria<br>Porifera", Phylum.ori)),
         Taxo.group = factor(Taxo.group, levels = c("Cnidaria<br>Porifera", "Annelida<br>Brachipoda<br>Nemertea", "Mollusca", "Arthropoda", "Echinodermata", "Chordata" )),
         Levels.group = ifelse(Levels %in% c("order", "phylum", "family", "class", "kingdom", "Unassigned"), "> genus or unassigned", Levels),
         Similar = ifelse(Taxon.ori == Taxon & Levels == "species", "Right identification",
                          ifelse(Levels == "species", "Wrong identification",
                                 ifelse(Genus.ori == genus & Levels == "genus", "Right identification",
                                        ifelse(Levels == "genus", "Wrong identification",
                                               #ifelse(Order.ori == order, "Order",
                                               "PROBLEMS"
                                        ))))) %>%
  group_by(method, Taxo.group,  Levels.group, Similar) %>% summarise(N = n()) %>%
  group_by(method, Taxo.group) %>% 
  summarise(Ntot = sum(N),
            freqGoodSp = N[Levels.group == "species" & Similar=="Right identification"]/Ntot ) %>% 
  
  #mutate(freq = N / sum(N)) %>% #View()
  #filter(Similar != "Wrong identification",
  #       Levels.group == "species") %>% 
  ggplot(aes(y=freqGoodSp, x = method)) +
  #geom_point(position=position_dodge(width=0.3), size = 3) +
  scale_y_continuous(limits = c(0,1)) +
  geom_bar(stat= "identity", position = "dodge", fill = "gray15") +
  scale_shape_manual(values = c(21,22))+
  #scale_fill_manual(values = c("gray15", "gray60", "gray95"))+
  facet_grid(. ~ Taxo.group , scale = "free") + 
  xlab("") +
  ylab(expression(paste("Proportion of accurate assigments"))) +  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom",
        strip.text = element_markdown(),
        legend.title = element_blank())

graphd 

graphd.over <- RES.postONncbi %>% filter(str_ends(Taxon.ori, " sp", negate = T)
                                         #Levels %in% c("species"),
) %>% 
  mutate(Taxo.group = ifelse(Phylum.ori %in% c("Sipuncula", "Nemertea", "Brachiopoda", "Annelida"), "Annelida<br>Brachipoda<br>Nemertea",
                             ifelse(Phylum.ori %in% c("Porifera", "Cnidaria"), "Cnidaria<br>Porifera", Phylum.ori)),
         Taxo.group = factor(Taxo.group, levels = c("Cnidaria<br>Porifera", "Annelida<br>Brachipoda<br>Nemertea", "Mollusca", "Arthropoda", "Echinodermata", "Chordata" )),
         Levels.group = ifelse(Levels %in% c("order", "phylum", "family", "class", "kingdom", "Unassigned"), "> genus or unassigned", Levels),
         Similar = ifelse(Taxon.ori == Taxon & Levels == "species", "Right identification",
                          ifelse(Levels == "species", "Wrong identification",
                                 ifelse(Genus.ori == genus & Levels == "genus", "Right identification",
                                        ifelse(Levels == "genus", "Wrong identification",
                                               #ifelse(Order.ori == order, "Order",
                                               "PROBLEMS"
                                        )))), 
         method = ifelse(method == "LCA", "<br>LCA<br>", method)) %>%
  group_by(method, threshold,  Levels.group, Similar) %>% summarise(N = n()) %>%
  group_by(method, threshold) %>% 
  summarise(Ntot = sum(N),
            freqGoodSp = N[Levels.group == "species" & Similar =="Right identification"]/Ntot ) %>% 
  
  #mutate(freq = N / sum(N)) %>% #View()
  #filter(Similar != "Wrong identification",
  #       Levels.group == "species") %>% 
  ggplot(aes(y=freqGoodSp, x = factor(threshold))) +
  #geom_point(position=position_dodge(width=0.3), size = 3) +
  scale_y_continuous(limits = c(0,1)) +
  geom_bar(stat= "identity", position = "dodge", fill = "gray15") +
  scale_shape_manual(values = c(21,22))+
  #scale_fill_manual(values = c("gray15", "gray60", "gray95"))+
  facet_grid(. ~ method , scale = "free") + 
  xlab("") +
  ylab(expression(paste("Proportion of accurate assigments"))) +  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom",
        strip.text = element_markdown(),
        legend.title = element_blank())

graphd.over 

ggarrange(graphd.over, graphd,
          labels = c("A", "B"),
          nrow = 1, ncol = 2, widths = c(2,5), align = "hv",
          common.legend = TRUE, legend = "bottom"
)


graph.comparaison <- ggarrange(grapha.over, grapha,
          graphb.over, graphb,
          #graphc.over, graphc,
          labels =  LETTERS[1:6],
          nrow = 2, ncol = 2, widths = c(2,5), align = "hv",
          common.legend = TRUE, legend = "bottom"
)
graph.comparaison

graph.comparaison.v2 <- ggarrange(grapha.over, grapha,
                               graphb.over, graphb,
                               grapha.acc.over, grapha.acc,
                               #graphc.over, graphc,
                               labels =  LETTERS[1:6],
                               nrow = 3, ncol = 2, widths = c(2,5), align = "hv",
                               common.legend = TRUE, legend = "bottom"
)
graph.comparaison.v2

ggsave(filename = file.path(here::here(), "03_Results", "fig_NCBI_tests.png"), 
       plot = graph.comparaison.v2, width = 9, height = 10, units = "in", bg = "white")
