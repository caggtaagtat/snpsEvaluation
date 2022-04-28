# surroundingSeq <-quickContextOutput[2]
# altSurroundingSeq <-quickContextOutput[3]
# coordinates <-quickContextOutput[4]
# coordinatesAlt <-quickContextOutput[5]


getAffectedSpliceSites <- function(surroundingSeq, altSurroundingSeq,
                                   coordinates, coordinatesAlt, SDcor, SAcor){
  
  splicesiteInfo <- list()
  surroundingSeqSDs <- substr(surroundingSeq, 40, nchar(surroundingSeq)-39)
  altSurroundingSeqSDs <- substr(altSurroundingSeq, 40, nchar(altSurroundingSeq)-39)
  subcoordSD <- strsplit(coordinates," ")[[1]]
  subcoordSD <- subcoordSD[c(40:(nchar(surroundingSeq)-39))]
  subcoordSDAlt <- strsplit(coordinatesAlt," ")[[1]]
  subcoordSDAlt <- subcoordSDAlt[c(40:(nchar(surroundingSeq)-39))]
  
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
  
  ## Remove annotated SD
  durchzahl$coordinate <- as.numeric(subcoordSD[1:13])+2
  durchzahl2$coordinate <- as.numeric(subcoordSDAlt[1:13])+2
  durchzahl <- durchzahl[durchzahl$coordinate != SDcor,]
  durchzahl2 <- durchzahl2[durchzahl2$coordinate != SDcor,]
  
  ## Save coordinate for SRP support calculation
  durchzahl2 <- durchzahl2[order(durchzahl2$hbs, decreasing = T),]
  splicesiteInfo$sd_coord <- durchzahl2$coordinate[1]
  
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
  
  durchzahl$coordinate <- as.numeric(subcoordSD[1:15])+2
  
  ## Remove annotated SD
  durchzahl$coordinate <- as.numeric(subcoordSD[1:15])+2
  durchzahl2$coordinate <- as.numeric(subcoordSD[1:15])+2
  durchzahl <- durchzahl[durchzahl$coordinate != SDcor,]
  durchzahl2 <- durchzahl2[durchzahl2$coordinate != SDcor,]
  
  ## Only select GC sites
  durchzahl$seq9  <- substr(durchzahl$seq9 , 4,5)
  durchzahl2$seq9 <- substr(durchzahl2$seq9, 4,5)
  durchzahl$hbs[durchzahl$seq9   != "GC"] <- 0
  durchzahl2$hbs[durchzahl2$seq9 != "GC"] <- 0
  
  ## Save coordinate for SRP support calculation
  durchzahl2 <- durchzahl2[order(durchzahl2$hbs, decreasing = T),]
  splicesiteInfo$sd_coordGC <- durchzahl2$coordinate[1]
  
  splicesiteInfo$sd_max_refGC <- max(durchzahl$hbs)
  splicesiteInfo$sd_max_altGC  <- max(durchzahl2$hbs)
  splicesiteInfo$sd_max_diffGC  <- splicesiteInfo$sd_max_altGC - splicesiteInfo$sd_max_refGC
  
  
  ## Searching overlapping SAs
  ## Calculate SA strengths
  
  ##Reduce string length
  surroundingSeqSAs <- substr(surroundingSeq, 27, nchar(surroundingSeq)-26)
  altsurroundingSeqSAs <- substr(altSurroundingSeq, 27, nchar(altSurroundingSeq)-26)
  subcoordSA <- strsplit(coordinates," ")[[1]]
  subcoordSA <- subcoordSA[c(27:(nchar(surroundingSeq)-26))]
  subcoordSAalt <- strsplit(coordinatesAlt," ")[[1]]
  subcoordSAalt <- subcoordSAalt[c(27:(nchar(surroundingSeq)-26))]
  
  durchzahl <- strsplit(altsurroundingSeqSAs, "")[[1]]
  d <- getOverlappingVectorsFromVector(durchzahl, 23, 22)
  d <- as.character(lapply(d, function(x){paste(x, collapse = "")}))
  durchzahl <- data.frame(seq9 = d, hbs=0)
  durchzahl$coordinates <- as.numeric(subcoordSAalt)+20
  durchzahl2 <- durchzahl[nchar(durchzahl$seq9) == 23,]
  
  durchzahl <- strsplit(surroundingSeqSAs, "")[[1]]
  d <- getOverlappingVectorsFromVector(durchzahl, 23, 22)
  d <- as.character(lapply(d, function(x){paste(x, collapse = "")}))
  durchzahl <- data.frame(seq9 = d, hbs=0)
  durchzahl$coordinates <- as.numeric(subcoordSA)+20
  durchzahl <- durchzahl[nchar(durchzahl$seq9) == 23,]
  
  ## Calculate MaxEntScan scores
  durchzahl2$maxent3  <-  as.numeric(calculateMaxEntScanScore(durchzahl2$seq9, 3))
  durchzahl$maxent3   <-  as.numeric(calculateMaxEntScanScore(durchzahl$seq9, 3))
  
  ## Remove annotated SD
  durchzahl <- durchzahl[durchzahl$coordinate != SAcor,]
  durchzahl2 <- durchzahl2[durchzahl2$coordinate != SAcor,]
  
  ## Save coordinate for SRP support calculation
  durchzahl2 <- durchzahl2[order(durchzahl2$maxent3, decreasing = T),]
  splicesiteInfo$sa_coord <- durchzahl2$coordinate[1]
  
  splicesiteInfo$sa_max_ref <- max(durchzahl$maxent3 )
  splicesiteInfo$sa_max_alt <- max(durchzahl2$maxent3 )
  splicesiteInfo$sa_max_diff <- splicesiteInfo$sa_max_alt - splicesiteInfo$sa_max_ref
  
  return(splicesiteInfo)
}
