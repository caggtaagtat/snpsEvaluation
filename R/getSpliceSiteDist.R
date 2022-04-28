getSpliceSiteDist <- function(transCur, genomic_coordinate, strand){
  
  #Convert strand annotation
  if(strand== "+") strand <- "1"
  if(strand== "-") strand <- "-1"
  
  subset_transcr <- transCur
  subset_transcr <- subset_transcr[order(subset_transcr$Exon.rank.in.transcript, decreasing = F),]
  
  ## Find SDs
  SDs <- subset_transcr$Exon.region.start..bp.[-nrow(subset_transcr)]
  if(strand == 1)  SDs <- subset_transcr$Exon.region.end..bp.[-nrow(subset_transcr)]
  
  ## Find SAs
  SAs <- subset_transcr$Exon.region.end..bp.[-1]
  if(strand == 1)  SAs <- subset_transcr$Exon.region.start..bp.[-1]
  
  ## Get nearest splice sites
  splicesites <- data.frame(coordinate=c(SDs, SAs), type= c(rep("SD", length(SDs)), rep("SA", length(SDs))))
  splicesites$distance <-  genomic_coordinate - splicesites$coordinate 
  if(strand == -1) splicesites$distance <- splicesites$distance * -1
  splicesites$distance_abs <- abs(splicesites$distance)
  splicesites <- splicesites[order(splicesites$distance_abs, decreasing = F),]
  splicesites$number <- 1:nrow(splicesites)
  
  return(splicesites)
}
