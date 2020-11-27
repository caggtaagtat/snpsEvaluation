getAffectedSpliceSites <- function(surroundingSeq, altSurroundingSeq ){
  
  splicesiteInfo <- list()
  surroundingSeqSDs <- substr(surroundingSeq, 40, nchar(surroundingSeq)-39)
  altSurroundingSeqSDs <- substr(altSurroundingSeq, 40, nchar(altSurroundingSeq)-39)
  
  ## Searching overlapping SDs
  ## Create 11nt long subsequences
  durchzahl <-  calculateHZEIperNT(surroundingSeqSDs)
  durchzahl$Sequence <- "reference"
  durchzahl2 <-  calculateHZEIperNT(altSurroundingSeqSDs)
  durchzahl2$Sequence <- "Mutated"
  
  ## Calculate HBS
  durchzahl$hbs <- hbg$hbs[match(durchzahl$seq9,hbg$seq)]
  durchzahl2$hbs <- hbg$hbs[match(durchzahl2$seq9,hbg$seq)]
  durchzahl[is.na(durchzahl)] <- 0    
  durchzahl2[is.na(durchzahl2)] <- 0    
  
  splicesiteInfo$sd_max_ref <- max(durchzahl$hbs)
  splicesiteInfo$sd_max_alt <- max(durchzahl2$hbs)
  splicesiteInfo$sd_max_diff <- splicesiteInfo$sd_max_alt - splicesiteInfo$sd_max_ref
  
  rm(durchzahl, durchzahl2)
  
  ## Searching overlapping SDs for GC sites
  ## Create 11nt long subsequences
  durchzahl <- strsplit(altSurroundingSeqSDs, "")[[1]]
  d <- getOverlappingVectorsFromVector(durchzahl, 9, 8)
  d <- as.character(lapply(d, function(x){paste(x, collapse = "")}))
  durchzahl <- data.frame(seq9 = d, hbs=0)
  durchzahl2 <- durchzahl[nchar(durchzahl$seq9) == 9,]
  
  durchzahl <- strsplit(surroundingSeqSDs, "")[[1]]
  d <- getOverlappingVectorsFromVector(durchzahl, 9, 8)
  d <- as.character(lapply(d, function(x){paste(x, collapse = "")}))
  durchzahl <- data.frame(seq9 = d, hbs=0)
  durchzahl <- durchzahl[nchar(durchzahl$seq9) == 9,]
  
  ## Calculate MES
  durchzahl$hbs <- as.numeric(calculateMaxEntScanScore(durchzahl$seq9, 5))
  durchzahl2$hbs <- as.numeric(calculateMaxEntScanScore(durchzahl2$seq9, 5))
  
  ## Only select GC sites
  durchzahl$seq9  <- substr(durchzahl$seq9 , 4,5)
  durchzahl2$seq9 <- substr(durchzahl2$seq9, 4,5)
  durchzahl$hbs[durchzahl$seq9   != "GC"] <- 0
  durchzahl2$hbs[durchzahl2$seq9 != "GC"] <- 0

  splicesiteInfo$sd_max_refGC <- max(durchzahl$hbs)
  splicesiteInfo$sd_max_altGC  <- max(durchzahl2$hbs)
  splicesiteInfo$sd_max_diffGC  <- splicesiteInfo$sd_max_altGC - splicesiteInfo$sd_max_refGC
  
  
  ## Searching overlapping SAs
  ## Calculate SA strengths
  
  ##Reduce string length
  surroundingSeqSAs <- substr(surroundingSeq, 27, nchar(surroundingSeq)-26)
  altsurroundingSeqSAs <- substr(altSurroundingSeq, 27, nchar(altSurroundingSeq)-26)
  
  durchzahl <- strsplit(altsurroundingSeqSAs, "")[[1]]
  d <- getOverlappingVectorsFromVector(durchzahl, 23, 22)
  d <- as.character(lapply(d, function(x){paste(x, collapse = "")}))
  durchzahl <- data.frame(seq9 = d, hbs=0)
  durchzahl2 <- durchzahl[nchar(durchzahl$seq9) == 23,]
  
  durchzahl <- strsplit(surroundingSeqSAs, "")[[1]]
  d <- getOverlappingVectorsFromVector(durchzahl, 23, 22)
  d <- as.character(lapply(d, function(x){paste(x, collapse = "")}))
  durchzahl <- data.frame(seq9 = d, hbs=0)
  durchzahl <- durchzahl[nchar(durchzahl$seq9) == 23,]
  
  ## Calculate MaxEntScan scores
  durchzahl2$maxent3  <-  as.numeric(calculateMaxEntScanScore(durchzahl2$seq9, 3))
  durchzahl$maxent3 <-  as.numeric(calculateMaxEntScanScore(durchzahl$seq9, 3))
  
  splicesiteInfo$sa_max_ref <- max(durchzahl$maxent3 )
  splicesiteInfo$sa_max_alt <- max(durchzahl2$maxent3 )
  splicesiteInfo$sa_max_diff <-splicesiteInfo$sa_max_alt - splicesiteInfo$sa_max_ref
  
  return(splicesiteInfo)
}
