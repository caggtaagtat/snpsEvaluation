convertTableForExcel <- function(snvTable){
  
  ## General info
  snvTable$Delta_HZEI_integral   <- round(snvTable$Delta_HZEI_integral,2)
  snvTable$SD_HBS_diff  <- round(snvTable$SD_HBS_diff,1)
  
  ## SRE support and changes
  snvTable$SA_SRE_support  <- round(snvTable$SA_SRE_support,2)
  snvTable$SA_SRE_support_alternative  <- round(snvTable$SA_SRE_support_alternative,2)
  snvTable$Delta_SA_SRE_support  <- round(snvTable$Delta_SA_SRE_support,2)
  snvTable$Delta_SA_SRE_support_percent  <- round(snvTable$Delta_SA_SRE_support_percent,2)
  snvTable$SD_SRE_support  <- round(snvTable$SD_SRE_support,2)
  snvTable$SD_SRE_support_alternative  <- round(snvTable$SD_SRE_support_alternative,2)
  snvTable$Delta_SD_SRE_support  <- round(snvTable$Delta_SD_SRE_support,2)
  snvTable$Delta_SD_SRE_support_percent  <- round(snvTable$Delta_SD_SRE_support_percent,2)
  
  ## Splice sites Info
  snvTable$MaxEnt_alt_overlapping_SAs  <-  round(snvTable$MaxEnt_alt_overlapping_SAs,2)
  snvTable$MaxEnt_diff_overlapping_SAs  <-  round(snvTable$MaxEnt_diff_overlapping_SAs,2)
  snvTable$MaxEnt_diff_overlapping_SA_vs_annotated  <- round(  snvTable$MaxEnt_diff_overlapping_SA_vs_annotated,2)
  snvTable$Maxent_alt_overlapping_SD_GC_sites  <-  round(snvTable$Maxent_alt_overlapping_SD_GC_sites,2)
  snvTable$Maxent_diff_overlapping_SD_GC_sites  <-  round(snvTable$Maxent_diff_overlapping_SD_GC_sites,2)

  snvTable$SA_MaxEnt_ref  <- round(snvTable$SA_MaxEnt_ref,2)
  snvTable$SA_MaxEnt_alt  <- round(snvTable$SA_MaxEnt_alt,2)
  
  snvTable$SRE_support_alt_overlapping_SD <- round(  snvTable$SRE_support_alt_overlapping_SD,2)
  snvTable$SRE_support_alt_overlapping_SA <- round(  snvTable$SRE_support_alt_overlapping_SA,2)
  snvTable$SRE_support_alt_overlapping_GCsite <- round(snvTable$SRE_support_alt_overlapping_GCsite,2)


  ## Replace "." in every column with an ",", since excle transforms 
  ## everything it can into a date for whatever reason
  ## General info
  #snvTable$Delta_HZEI_integral   <- gsub("\\.", ",", snvTable$Delta_HZEI_integral)
  # snvTable$HBS_alt_overlapping_SDs  <- snvTable$HBS_alt_overlapping_SDs)
  # snvTable$HBS_diff_overlapping_SDs  <- snvTable$HBS_diff_overlapping_SDs)
  # snvTable$HBS_diff_overlapping_SDs_vs_annotated <- snvTable$HBS_diff_overlapping_SDs_vs_annotated )
  # snvTable$SD_HBS_ref  <-  snvTable$SD_HBS_ref)
  # snvTable$SD_HBS_alt  <-  snvTable$SD_HBS_alt)
  
  snvTable$gene_strand[snvTable$gene_strand=="1"] <- "+"
  snvTable$gene_strand[snvTable$gene_strand=="-1"] <- "-"
  
  snvTable$chromosome <- paste0("chr",snvTable$chromosome)
  
  return(snvTable)
}
