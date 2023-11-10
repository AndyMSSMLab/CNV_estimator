suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(fst))

args <- commandArgs(trailingOnly=TRUE)
input <- args[1]

convertToFST <- function(file){
  # converts the output of mosdepth into FST to reduce the file size
  x <- fread(file, header = F, sep = "\t") %>% as.data.frame %>%
    select(V4, V5) %>%
    mutate(V5 = round(V5*100) ) %>%
    rename(key = 1, RD = 2) 
  
  outfile <- sub(".gz$",".fst",basename(file))
  write_fst(x, outfile, compress = 100, uniform_encoding = TRUE)
}

addCoordToKey <- function(file, key_to_coord){
  x <- read_fst(file) %>% rename(this_key = key) %>% arrange(this_key)
  combined <- cbind(key_to_coord, x)
  test <- combined %>% filter(key!=this_key)
  if(nrow(test)!=0 && nrow(x)!=nrow(key_to_coord)){
    stop(paste0("Error: ", file, " and ", key_to_coord, " are not compatible"))
  }else {
    combined %>% select(-this_key) %>% 
      mutate(RD=RD/100)
  #    mutate(RD = as.numeric(formatC(RD/100, format = 'f', flag='0', digits = 2)))
  }
  
}

calGCcorrFac <- function(data,which.chr = c("Auto", "X", "Y")){

  ## Calculate GC correction factors for 100bp bins on Autosomes and chrX/Y
  if(which.chr == "Auto"){
    chrToSelect <- data %>% select(1) %>% distinct %>% 
       filter(grepl("^chr[0-9]+$",chr)) %>% pull
  }else if(which.chr == "X") {
    chrToSelect <- data %>% select(1) %>% distinct %>% 
      filter(grepl("^chrX$",chr)) %>% pull
  }else if(which.chr == "Y") {
    chrToSelect <- data %>% select(1) %>% distinct %>% 
      filter(grepl("^chrY$",chr)) %>% pull
  }
  x <- data %>% filter(chr %in% chrToSelect) %>% select(-key)
  
  meanRD_highGC <- x %>% filter(GC> 89) %>% summarize(median = median(RD), n = n()) %>%
    ungroup %>% as.data.frame
  
  y <- x %>%  
    group_by(GC) %>% summarize(median = median(RD), n = n()) %>%
    ungroup %>% as.data.frame %>% arrange(GC) %>%
    mutate(median = ifelse(GC>89, meanRD_highGC$median,median)) %>%
    mutate(n = ifelse(GC>89, meanRD_highGC$n, n)) 
  
  if(which.chr == "Auto"){ 
    y = y %>% filter(n>=100)
  }else {
    y = y %>% filter(n>=100 | GC>89)
  }

  missing_bins <- as.numeric(setdiff(as.character(1:89), as.character(y$GC)))
  
  missing_bins_rd <- lapply(missing_bins, function(bin){
    upper = y %>% filter(GC>bin) %>% arrange(GC) %>% slice(1)
    lower = y %>% filter(GC<bin) %>% arrange(desc(GC)) %>% slice(1)
    if(nrow(lower) == 0){ lower = data.frame(GC=-1,median=-1,n=-1)}
    if(nrow(upper) == 0){ upper = data.frame(GC=1000000,median=-1,n=-1)}
    
    x  %>%
      filter(GC >=lower$GC & GC<=upper$GC) %>%
      summarize(median = median(RD), n = n()) %>%
      mutate(GC = bin) %>% select(GC, median, n) %>% as.data.frame
  })
  
  missing_bins_rd <- do.call(rbind, missing_bins_rd)
  
  missing_bins_rd_highGC <- cbind(
    GC = as.numeric(setdiff(as.character(90:100), as.character(y$GC))),
    y %>% filter(GC>89) %>% select(-GC) %>% distinct)
  
  if(length(missing_bins)>0){
    missing_bins_rd <- missing_bins_rd %>% mutate(Type = "P")
  }
  
  z <- rbind(y %>% mutate(Type = "O"),
                 missing_bins_rd,
                 missing_bins_rd_highGC %>% mutate(Type = "P")) %>%
    arrange(GC) %>%
    mutate(corr_fact = round(median(x$RD)/median,4)) %>%
    mutate(global_mean = round(mean(x$RD),4)) %>%
    mutate(type = which.chr)
  z
}


medianPerChr <- function(x, chrXY_100bin_file) {
  #Calculate stats for Autosomes and non PAR regions on chrX
  x1 <- x  %>%
    group_by(chr) %>% 
    summarize(
      median = median(RD), 
      mean = round(mean(RD),3),
      n = n()
    ) %>%
    ungroup %>% as.data.frame
  
  ## removing PAR1 and PAR2
  chrXYbins <- fread(chrXY_100bin_file, header = F, sep = "\t") %>% 
    as.data.frame %>% rename(chr=1, start=2, end=3, GC =4) %>%
    select(chr, start, end)
  x2 <- x %>% 
    inner_join(chrXYbins, by = c("chr","start","end")) %>%
    group_by(chr) %>% 
    summarize(
      median = median(RD), 
      mean = round(mean(RD),3),
      n = n()
    ) %>%
    ungroup %>% as.data.frame %>%
    mutate(chr = ifelse(chr == "chrX", "chrX_noPAR", 
                        ifelse(chr == "chrY", "chrY_noPAR", chr)))
  
  x3 <- x  %>%
    filter(grepl("^chr[0-9]+$",chr)) %>% 
    summarize(
      chr = "Autosomes",
      median = median(RD), 
      mean = round(mean(RD),3),
      n = n()
    ) %>%
    ungroup %>% as.data.frame
  out <- rbind(x1, x2,x3)
  out
}


normalizeRD <- function(x, sample_name, hap=2, which.chr = c("Auto", "X", "Y"),
                        gc_file, regionsToMerge){
  ### Peform normalization on raw read depth to generate CN estimates
  if(which.chr == "Auto"){
    chrToSelect <- data %>% select(1) %>% distinct %>% 
      filter(grepl("^chr[0-9U]+",chr)) %>% pull
  }else if(which.chr == "X") {
    chrToSelect <- data %>% select(1) %>% distinct %>% 
      filter(grepl("^chrX$",chr)) %>% pull
  }else if(which.chr == "Y") {
    chrToSelect <- data %>% select(1) %>% distinct %>% 
      filter(grepl("^chrY$",chr)) %>% pull
  } else { stop(paste0("Error: option ",which.chr," not found !!!") )}
  
  x <- x %>% rename(chr_100 = chr, start_100 = start, end_100 = end) %>%
    mutate(GC = as.character(GC))
  
  gc <- fread(gc_file,
             header =T, sep ="\t") %>% as.data.frame %>% 
    filter(type == ifelse(which.chr == "Auto", "Auto", "X")) %>%
    select(GC, corr_fact, global_mean) %>%
    mutate(GC = as.character((GC))) 
  
  out <- fread(regionsToMerge, header = F, sep = "\t") %>% as.data.frame %>%
    rename(chr = 1, start =2, end = 3, geneNames = 4, Type = 5,
           chr_100 = 6, start_100 = 7, end_100 =8, GC = 9) %>%
    select(chr:end_100) %>%
    filter(chr %in% chrToSelect) %>%
    inner_join(x,by = c("chr_100", "start_100", "end_100")) %>%
    mutate(GC = as.character(GC)) %>%
    inner_join(gc, by = "GC") %>%
    mutate(CN = round(hap* RD * corr_fact/global_mean,4)) %>%
    mutate(start_100 = ifelse(start_100<start, start, start_100)) %>%
    mutate(end_100 = ifelse(end_100>end, end, end_100)) %>%
    group_by(geneNames, Type) %>% mutate(n = n(), end_100 = end_100-1) %>%
    mutate(frac = ifelse(end == end_100, (end_100- start_100)/100,(end_100- start_100 + 1)/100)) %>%
    mutate(CN_frac = CN* frac) %>%
    summarize(CN = sum(CN_frac)/(sum(frac))) %>% ungroup %>%
    as.data.frame %>%
    mutate(sampleID = sample_name) %>%
    select(geneNames, sampleID, CN, Type) %>%
    mutate(CN = round(CN,4))
  out
}


getRawRD <- function(x, regionsToMerge, type, sample_name){
  x <- x %>% rename(chr_100 = chr, start_100 = start, end_100 = end) %>%
    mutate(GC = as.character(GC))
  
  data_raw_ribo <- fread(regionsToMerge, header = F, sep = "\t") %>% as.data.frame %>%
    rename(chr = 1, start =2, end = 3, geneNames = 4, Type = 5,
           chr_100 = 6, start_100 = 7, end_100 =8, GC = 9) %>%
    filter(Type == type) %>%
    inner_join(x,by = c("chr_100", "start_100", "end_100")) %>%
    mutate(CN = RD) %>%
    mutate(start_100 = ifelse(start_100<start, start, start_100)) %>%
    mutate(end_100 = ifelse(end_100>end, end, end_100)) %>%
    group_by(geneNames, Type) %>% mutate(n = n(), end_100 = end_100-1) %>%
    mutate(frac = ifelse(end == end_100, (end_100- start_100)/100,(end_100- start_100 + 1)/100)) %>%
    mutate(CN_frac = CN* frac) %>%
    summarize(CN = sum(CN_frac)/(sum(frac))) %>% ungroup %>%
    as.data.frame %>%
    mutate(sampleID = sample_name) %>%
    select(geneNames, sampleID, CN, Type) %>%
    mutate(CN = round(CN,4)) %>% rename(meanRD = CN)
}


input.mosdepth.file <- args[1]
key_to_coord_file <- args[2]
chrXY_100bin_file <- args[3]
regionsToMerge <- args[4]
gender <- args[5]

### Check if input file is .gz format
### if yes, convert them to fst format
if(grepl(".gz", input.mosdepth.file)){
  convertToFST(input.mosdepth.file)
}

input.mosdepth.file <- sub(".gz$",".fst", basename(input.mosdepth.file))

#key_to_coord_file <- "../chr_100bp.gc.bed.gz"
key_to_coord <- fread(key_to_coord_file, header = F, sep = "\t") %>%
  as.data.frame %>%
  rename(chr = 1, start =2, end = 3, GC = 4, key = 5) %>%
  arrange(key)

data <- addCoordToKey(input.mosdepth.file, key_to_coord)


#### Calculating GC correction factors
auto_gc_corfac <- calGCcorrFac(data, which.chr = "Auto")
chrx_gc_corfac <- calGCcorrFac(data, which.chr = "X")
outfile <- sub(".bed.fst",".GCbins.bed.gz", basename(input.mosdepth.file))
fwrite(
  rbind(auto_gc_corfac,chrx_gc_corfac),
  outfile,
  row.names = F, sep = "\t", quote = F
)

#### Calculating median and per per chromosome
stats <- medianPerChr(data, chrXY_100bin_file)
outfile <- sub(".bed.fst",".meanRD.perChr.bed.gz", basename(input.mosdepth.file))
fwrite(
  stats,
  outfile,
  row.names = F, sep = "\t", quote = F
)


#### Performing normaliztion
gc_file <- sub(".bed.fst",".GCbins.bed.gz", basename(input.mosdepth.file))
sample_name <- sub(".regions.bed.fst","", basename(input.mosdepth.file))
auto <- normalizeRD(data, sample_name, hap=2, which.chr = "Auto", gc_file, regionsToMerge)
if(gender == "female"){
  chrx <- normalizeRD(data, sample_name, hap=2, which.chr = "X", gc_file, regionsToMerge)
}else{
  chrx <- normalizeRD(data, sample_name, hap=1, which.chr = "X", gc_file, regionsToMerge)
}

chry <- normalizeRD(data, sample_name, hap=1, which.chr = "Y", gc_file, regionsToMerge)
normCN <- rbind(auto, chrx, chry) 
outfile <- sub(".regions.bed.fst",".norm.txt.gz", basename(input.mosdepth.file))
fwrite(
  normCN,
  outfile,
  row.names = F, sep = "\t", quote = F
)
