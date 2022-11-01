# FUnctions to play around sequences data

# Check start and stop position - as long as they are A C T G

start.pos <- function(DNA.vec){
  SEQ.start <- data.frame(
    id = seq_along(DNA.vec),
    A.pos = DNA.vec %>% str_locate(c("A")) %>% as_tibble() %>% pull("start"),
    C.pos = DNA.vec %>% str_locate(c("C")) %>% as_tibble() %>% pull("start"),
    T.pos = DNA.vec %>% str_locate(c("T")) %>% as_tibble() %>% pull("start"),
    G.pos = DNA.vec %>% str_locate(c("G")) %>% as_tibble() %>% pull("start")
  )
  
  SEQ.start <- SEQ.start %>% group_by(id) %>% mutate(Start.pos = min(A.pos, C.pos, T.pos, G.pos))
  
  res <- SEQ.start %>% pull(Start.pos)
  
  return(res)
  
}

stop.pos <- function(DNA.vec){
  SEQ.stop <- data.frame(
    id = seq_along(DNA.vec),
    A.pos = DNA.vec %>%  stringi::stri_locate_last_coll(c("A")) %>% as_tibble() %>% pull("start"),
    C.pos = DNA.vec %>%  stringi::stri_locate_last_coll(c("C")) %>% as_tibble() %>% pull("start"),
    T.pos = DNA.vec %>%  stringi::stri_locate_last_coll(c("T")) %>% as_tibble() %>% pull("start"),
    G.pos = DNA.vec %>%  stringi::stri_locate_last_coll(c("G")) %>% as_tibble() %>% pull("start")
  )
  
  SEQ.stop <- SEQ.stop %>% group_by(id) %>% mutate(Stop.pos = max(A.pos, C.pos, T.pos, G.pos))
  
  res <- SEQ.stop %>% pull(Stop.pos)
  
  return(res)
  
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}