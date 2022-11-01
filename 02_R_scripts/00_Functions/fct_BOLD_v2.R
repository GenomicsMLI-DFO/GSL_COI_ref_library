# Library -----------------------------------------------------------------

if(!require(tidyverse)){install.packages("tidyverse")}
library(tidyverse)

if(!require(bold)){install.packages("bold")}
library(bold)
packageVersion("bold")
citation("bold")

library(Biostrings)
library(msa)



# Values ------------------------------------------------------------------

  # Info on geographical location


  #min.l.pb = 650
  
# Function ----------------------------------------------------------------

#  Species <- c("Artediellus uncinatus", "Salmo salar", "Acipenser oxyrinchus", "Alepocephalus bairdii", "Colga villosa", 
#               "Isodictya deichmannae", 
#               "Isodictya palmata",
#               "Lissodendoryx indistincta",
#               "Mycale lingua",
#               "Mycale loveni",
#               "Myxilla incrustans",
# "Aulactinia stella"
#  ) 

#x <- "Leptasterias polaris"

BOLD.sratch <- function(Species, taxo.level = "Species", min.l.pb){

  # To give a priority to some seq over others
  GSL_abbr = c("St-Laurent", "Saint-Laurent", "St Lawrence", "St. Lawrence")
  NWA_abbr = c("Fundy", "Flemish Cap", "Maritimes", "Scotian Shelf", "Georges Bank", "Grand Bank", "Baffin Bay", "Durban Harbour")
  
  CAN_province = c("Quebec", "Nova Scotia", "New Brunswick", "Newfoundland and Labrador", "Prince Edward Island")
  US_state = c("Maine", "New Hampshire", "Massachusetts", "Rhode Island", "Connecticut", "New Jersey", "Delaware", "Maryland", "Virginia")
  
  Lat.value = c(35,78)
  Long.value = c(-78,-40)
  
  NOT_abbr = c("Hudson Bay", "Hudson Strait","Beaufort", "Nord-du-Quebec", "Ungava")
  
  
  
  
  # Empty results
  summary.stats <- data.frame(species = character(),
                              n.BOLD.raw = numeric(),
                              n.flag.raw = numeric(),
                              n.BOLD.filt = numeric(),
                              n.BOLD.filt.NWA = numeric(),
                              #n.DNA.final = numeric(),
                              # Compute stats on BINS
                              n.BIN.raw = numeric(),
                              n.BIN.raw.NWA = numeric()
                              

                              #n.BIN.final = numeric(),
                              #n.BIN.final.NWA = numeric(),                              
                              # Compute stats on SEQ
                              #Seq.start = numeric(),
                              #Seq.stop = numeric(),
                              #Alig.length = numeric()
                              )
  
  #DNA.unique <- DNAStringSet()
  
  res.raw <- data.frame()
  res.final <- data.frame()
  
  for(x in Species) {
    
    cat("\nProcessing ", x, "\n", sep= "")
    # Bold research
    res <- bold_seqspec(taxon = x)

    # If BOLD do not return a result
    if(is.null(nrow(res))) {
      
      summary.stats <- rbind(summary.stats,
                             data.frame(species = x,
                                        n.BOLD.raw = 0,
                                        n.flag.raw = 0,
                                        n.BOLD.filt = 0,
                                        n.BOLD.filt.NWA = 0,
                                        #n.DNA.final = 0,
                                        
                                        n.BIN.raw = 0,
                                        n.BIN.raw.NWA = 0
                                        
                                        #n.BIN.final = 0,
                                        #n.BIN.final.NWA = 0,   
                                        

                                        #Seq.start = NA,
                                        #Seq.stop = NA,
                                        #Alig.length = NA
                                        
                             )
      )
      
    }
    
    # If BOLD return a result
    if(isTRUE(nrow(res)>=1)){
      
      
      if(taxo.level == "Genus"){
        res <- res %>% filter(genus_name == x)
      }
      
      
      res <- res %>% filter(markercode == "COI-5P"#,
                     #species_name == species | subspecies_name == species
      ) %>% 
        
        mutate(nucleotides_ok = str_remove_all(nucleotides, "-"), 
               #Compute fragment length
               N_nuc = nchar(nucleotides_ok),
               #Compute ambiguous sites
               N_nuc_A = str_count(nucleotides_ok, "A"),
               N_nuc_C = str_count(nucleotides_ok, "C"),
               N_nuc_T = str_count(nucleotides_ok, "T"),
               N_nuc_G = str_count(nucleotides_ok, "G"),
               N_nuc_ambi = N_nuc - (N_nuc_A + N_nuc_C + N_nuc_T + N_nuc_G),
               Perc_nuc_ambi = N_nuc_ambi / N_nuc,
               # Give a priority given the % of ambiguous site
               Nuc_code = ifelse(Perc_nuc_ambi == 0, 1, 
                                 ifelse(Perc_nuc_ambi < 0.01,2,3)),
               # Give a priority given the geographic location
               Geo_code =   ifelse(# Voir les mots-clés a éviter en premier
                                   str_detect(region, paste(NOT_abbr, collapse = "|")) | 
                                   str_detect(sector, paste(NOT_abbr, collapse = "|")) |
                                   str_detect(exactsite, paste(NOT_abbr, collapse = "|")), 2,
                            ifelse(# Voir les mots-clés chosis       
                                   str_detect(region, paste(c(NWA_abbr, GSL_abbr), collapse = "|")) | 
                                   str_detect(sector, paste(c(NWA_abbr, GSL_abbr), collapse = "|")) |
                                   str_detect(exactsite, paste(c(NWA_abbr, GSL_abbr), collapse = "|"))|
                                   str_detect(province_state, paste(c(CAN_province, US_state),collapse = "|")), 1,
                            2)),
               lat = as.numeric(as.character(lat)),
               lon = as.numeric(as.character(lon)),
               Coor_code = ifelse( between(lat, Lat.value[1], Lat.value[2]) &
                                   between(lon, Long.value[1], Long.value[2]), 1,2  
                                     ),
               NWA_code = ifelse(Coor_code == 1 | Geo_code == 1, 1,2),
               NWA_code = ifelse(is.na(NWA_code),2, NWA_code),
               Image_code = ifelse(image_ids == "",2,
                                   ifelse(!is.na(image_ids),1,2)),
               Check_BIN =  ifelse(str_detect(bin_uri, "BOLD"),1,0)
               ) %>% # Because something a have a few problems        
               mutate(image_ids = as.character(image_ids),
                                      copyright_years = as.character(copyright_years),
                                      site_code = as.character(site_code),
                                      collectiontime = as.character(collectiontime),
                                      exactsite = as.character(exactsite),
                                      coord_source = as.character(coord_source),
                                      sex = as.character(sex),
                                      extrainfo = as.character(extrainfo),
                                      fieldnum = as.character(fieldnum),
                                      trace_ids = as.character(trace_ids),
                                      image_ids = as.character(image_ids),
                                      media_descriptors = as.character(media_descriptors),
                                      directions = as.character(directions),
                                      Species.search = x)
      
      #res %>% select(lat, lon, country, province_state, region, sector, exactsite,Geo_code, Coor_code, NWA_code) %>% View()
      
      # Stats
      n.BOLD.raw.int <- res %>% nrow()  
      
      n.BIN.raw.int <- res %>% filter(str_detect(bin_uri, "BOLD")) %>%  pull(bin_uri) %>% unique() %>% length()
      n.BIN.NWA.raw.int <- res %>% filter(str_detect(bin_uri, "BOLD"),
                                          NWA_code %in% c(1)) %>%  
                                   pull(bin_uri) %>% unique() %>% length()
      
      
      n.flag.raw.int <- res %>% filter(str_detect(bin_uri, "BOLD", negate = TRUE)) %>% nrow()
      
      # Filter for the right marker
      res.filter <- res %>% 
        filter(N_nuc >= min.l.pb,
               Check_BIN == 1,
               Nuc_code %in% c(1:2)) %>%
        # Order species
        arrange(NWA_code, Nuc_code, Perc_nuc_ambi, Image_code, desc(N_nuc)) # desc to be from larger to shorter ...
      
      #
      
      # What to do if filtration remove everything
      if(nrow(res.filter) == 0){
        
        summary.stats <- rbind(summary.stats,
                               data.frame(species = x,
                                          n.BOLD.raw = n.BOLD.raw.int,
                                          n.flag.raw = n.flag.raw.int,
                                          n.BOLD.filt = 0,
                                          n.BOLD.filt.NWA = 0,
                                          #n.DNA.final = 0,
                                          
                                          n.BIN.raw = n.BIN.raw.int,
                                          n.BIN.raw.NWA = n.BIN.NWA.raw.int
                                          
                                          #n.BIN.final = 0,
                                          #n.BIN.final.NWA = 0, 
                                          
                                          #Seq.start = NA,
                                          #Seq.stop = NA,
                                          #Alig.length = NA
                               )
        )
      }
      
      if(nrow(res.filter) >= 1){ 
        
        # Add a column to keep an eye on the priority - from arrange function above
        res.filter$PriorityOrder <- 1:nrow(res.filter)
        
        # Compute info on N sequences following the filtration
        n.BOLD.filt.int <- res.filter %>% nrow()
        n.BOLD.filt.NWA.int <- res.filter %>% filter(NWA_code %in% c(1)) %>% nrow()
        
        
        res.raw   <- bind_rows(res.raw, res)
        
        res.final <- bind_rows(res.final, res.filter)
        
        
        
        # Compute stats
        summary.stats <- rbind(summary.stats,
                               data.frame(species = x,
                                          n.BOLD.raw = n.BOLD.raw.int,
                                          n.flag.raw = n.flag.raw.int,
                                          n.BOLD.filt = n.BOLD.filt.int,
                                          n.BOLD.filt.NWA = n.BOLD.filt.NWA.int,
                                          #n.DNA.final = n.DNA.final.int,
                                          
                                          n.BIN.raw = n.BIN.raw.int,
                                          n.BIN.raw.NWA = n.BIN.NWA.raw.int
                                          
                                          #n.BIN.final = n.BIN.final.int,
                                          #n.BIN.final.NWA = n.BIN.final.NWA.int, 
                                          
                                          #Seq.start = cut.F.int,
                                          #Seq.stop = cut.R.int,
                                          #Alig.length = Alig.length.int
                                          
                               )
        )
        
      }
      
    } }
  # End of the loop over species
  
  return(list(Summary = summary.stats, Meta.Raw = res.raw, Meta.filt = res.final))
}

  

BOLD.seq <- function(res.filt, min.l.pb){
  
  
  # Info on the locus 
  
  LCO1490 <- "GGTCAACAAATCATAAAGATATTGG"
  HCO2198 <- "TAAACTTCAGGGTGACCAAAAAATCA"
  
  n.mist = 7
  
  
    # Empty results
    summary.stats <- data.frame(species = character(),

                                n.DNA.final = numeric(),
                                # Compute stats on BINS
                               
                                n.BIN.final = numeric(),
                                n.BIN.final.NWA = numeric(),                              
                                # Compute stats on SEQ
                                Seq.start = numeric(),
                                Seq.stop = numeric(),
                                Alig.length = numeric()
    )
    
    DNA.unique <- DNAStringSet()
    
    res.final <- data.frame()
    
    Species <- res.filt %>% pull(Species.search) %>% unique()
    
    for(x in Species) {
      
      cat("\nProcessing ", x, "\n", sep= "")
  
      res.filter <- res.filt %>% filter(Species.search == x)
      
      
  # HERE IS A SECOND STEP
  # split by data
  
  DNA.int.NWA <-        DNAStringSet(res.filter %>% filter(NWA_code %in% c(1)) %>% pull(nucleotides_ok))
  if(length(DNA.int.NWA) > 0){
    names(DNA.int.NWA) <- res.filter %>% filter(NWA_code %in% c(1)) %>% pull(PriorityOrder)       
    
  }
  
  DNA.int.OTHER <-        DNAStringSet(res.filter %>% filter(NWA_code %in% c(2)) %>% pull(nucleotides_ok))
  if(length(DNA.int.OTHER) > 0){
    names(DNA.int.OTHER) <- res.filter %>% filter(NWA_code %in% c(2)) %>% pull(PriorityOrder)       
  }
  
  if(length(DNA.int.NWA) > 50){
    DNA.int.NWA <- DNA.int.NWA[1:50] 
  }
  
  if(length(DNA.int.OTHER) > 20){
    DNA.int.OTHER <- DNA.int.OTHER[1:20] 
  }        
  
  if(length(DNA.int.NWA) > 0){
    DNA.int <- DNA.int.NWA
  } else {
    DNA.int <- DNA.int.OTHER         
  }
  

  
  
  # The alignment if more than 1 sequence
  if(length(DNA.int) == 1){
    DNA.alig <- DNA.int
  }
  
  if(length(DNA.int) > 1){
    DNA.alig <- msa(DNA.int, method = "ClustalW")
    
  }
  
  # The search for primers
  res.F <- vmatchPattern(DNAString(LCO1490), DNAStringSet(DNA.alig), max.mismatch = n.mist)
  #print(res.F@ends %>% unlist() %>% table())
  
  if(length(res.F@ends %>% unlist() %>% table()) == 1){
    cut.F.int <- dimnames(res.F@ends %>% unlist() %>% table())[[1]]
    cut.F.int <- as.numeric(as.character(cut.F.int)) + 1 # to exclude the primer
  } else {
    cut.F.int <- 1 
  }
  
  res.R <- vmatchPattern(reverseComplement(DNAString(HCO2198)), DNAStringSet(DNA.alig), max.mismatch = n.mist)
  
  if(length(res.R@ends %>% unlist() %>% table()) == 1){
    cut.R.int <- dimnames(res.R@ends %>% unlist() %>% table())[[1]]
    cut.R.int <- as.numeric(as.character(cut.R.int)) - nchar(HCO2198)  # to exclude the primer
    
  } else {
    cut.R.int <- DNAStringSet(DNA.alig)@ranges@width %>% unique() 
  }
  
  # What if there is a potential error and the fragment will be too short
  if(cut.R.int - cut.F.int < min.l.pb){
    cut.F.int <- 1
    cut.R.int <- DNAStringSet(DNA.alig)@ranges@width %>% unique() 
  }
  
  DNA.cut <- subseq(DNAStringSet(DNA.alig), start = cut.F.int, end = cut.R.int)
  DNA.cut <- DNA.cut[base::order(names(DNA.cut) %>% as.numeric())]
  
  Alig.length.int <- DNA.cut@ranges@width %>% unique()
  
  # Check if some sequences are now too short
  
  new.length.int <- DNA.cut %>% as.data.frame %>% mutate(x = str_remove_all(x, "-"),
                                                         N = nchar(x)) %>% pull(N)
  DNA.cut.long <- DNA.cut[new.length.int >= min.l.pb]
  
  # Create a temporary shorter sequence stuff to check for duplicates
  DNA.dup.temp <- subseq(DNA.cut.long, start = 3, end = Alig.length.int - 3)
  DNA.index.1 <- table(DNA.dup.temp) %>% as.data.frame() %>% 
    mutate(SEQ_name = paste("SEQ", 1:nrow(.), sep = ""))
  
  DNA.index.2 <-  DNA.dup.temp %>% as.data.frame() %>% 
    mutate(PriorityOrder = names(DNA.dup.temp)) %>% 
    left_join(DNA.index.1,
              by = c("x" = "DNA.dup.temp")) %>% 
    mutate(SEQ_details = paste(SEQ_name, paste0("n", Freq), sep = "_"),
           PriorityOrder = as.numeric(as.character(PriorityOrder))) %>% 
    select(-x)
  
  
  DNA.unique.int <- DNA.cut.long[DNA.dup.temp %>% duplicated() == F]
  
  
  
  n.DNA.final.int <- length(DNA.unique.int)
  
  
  # Get info back 
  res.unique <- res.filter %>% filter(PriorityOrder %in% names(DNA.unique.int)) %>% 
    left_join(DNA.index.2, by = "PriorityOrder")
  
  
  
  if(identical(as.character(names(DNA.unique.int)), as.character(res.unique$PriorityOrder))){
    names(DNA.unique.int) <- paste(res.unique$processid, x, res.unique$bin_uri, res.unique$SEQ_details, sep = "_")
    
  } else {
    print("There is a problem with the names of sequences (wrong order")
  }
  
  n.BIN.final.int <- res.unique %>% filter(str_detect(bin_uri, "BOLD")) %>%  pull(bin_uri) %>% unique() %>% length()
  n.BIN.final.NWA.int <- res.unique %>% filter(str_detect(bin_uri, "BOLD"),
                                               Geo_code %in% c(1:2)) %>%  
                                         pull(bin_uri) %>% unique() %>% length()
  
  # Compute stats
  summary.stats <- rbind(summary.stats,
                         data.frame(species = x,
                                    #n.BOLD.raw = n.BOLD.raw.int,
                                    #n.flag.raw = n.flag.raw.int,
                                    #n.BOLD.filt = n.BOLD.filt.int,
                                    #n.BOLD.filt.NWA = n.BOLD.filt.NWA.int,
                                    n.DNA.final = n.DNA.final.int,
                                    
                                    #n.BIN.raw = n.BIN.raw.int,
                                    #n.BIN.raw.NWA = n.BIN.NWA.raw.int
                                    
                                    n.BIN.final = n.BIN.final.int,
                                    n.BIN.final.NWA = n.BIN.final.NWA.int, 
                                    
                                    Seq.start = cut.F.int,
                                    Seq.stop = cut.R.int,
                                    Alig.length = Alig.length.int
                                    
                         )
  )
  
  
  res.final <- bind_rows(res.final, res.unique)
  
  DNA.unique <- c(DNA.unique, DNA.unique.int)
  
    }
  
  return(list(Summary = summary.stats, Meta.seq = res.final, DNA = DNA.unique))
}
