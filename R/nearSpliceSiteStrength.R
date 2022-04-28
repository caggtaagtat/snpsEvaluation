nearSpliceSiteStrength <- function(transCur, gen_cord, strand, context){
  
  #Convert strand annotation
  if(strand== "+") strand <- "1"
  if(strand== "-") strand <- "-1"
  
  second_up_ss_strength   <- 9999
  second_down_ss_strength <- 9999
  
  ## If located within exons
  if(context == "exonic"){

    ## Find current Exon
    withinIndex <- gen_cord >= transCur$Exon.region.start..bp. & 
    gen_cord <= transCur$Exon.region.end..bp.
    cur_rank <-  transCur$Exon.rank.in.transcript[withinIndex]
      
    ## Access saved ss strengths
    upIndex <- transCur$Exon.rank.in.transcript == cur_rank-1 
    downIndex <- transCur$Exon.rank.in.transcript == cur_rank+1
    if(cur_rank != 1 ) second_up_ss_strength   <- paste0("SD: ",transCur$donor_seq_HBS[upIndex])
    if(cur_rank !=  max(transCur$Exon.rank.in.transcript)) second_down_ss_strength  <- paste0("SA: ",transCur$acc_seq_MAXENT[downIndex])

  }
  
  ## If located within intron
  if(context == "intronic"){
    
    ## If located on plus strand
    if(strand==1){
      
      ## Get upstream exon rank 
      indexEnd    <-  gen_cord > transCur$Exon.region.end..bp. 
      cur_rank_up <-  max(transCur$Exon.rank.in.transcript[indexEnd])
      
      ## Get downstream exon rank 
      indexStart    <-  gen_cord < transCur$Exon.region.start..bp.
      cur_rank_down <-  min(transCur$Exon.rank.in.transcript[indexStart])
      
      indexUp   <- transCur$Exon.rank.in.transcript == cur_rank_up
      indexDown <- transCur$Exon.rank.in.transcript == cur_rank_down
      second_up_ss_strength   <- paste0("SA: ",transCur$acc_seq_MAXENT[indexUp])
      second_down_ss_strength  <- paste0("SD: ",transCur$donor_seq_HBS[indexDown])
      
    ## If on minus strand
    }else{
      
      ## Get upstream exon rank 
      indexEnd    <-  gen_cord < transCur$Exon.region.end..bp.
      indexStart  <-  gen_cord < transCur$Exon.region.start..bp.
      cur_rank_up <-  max(transCur$Exon.rank.in.transcript[indexEnd & indexStart])
      
      ## Get downstream exon rank 
      indexEnd    <-  gen_cord > transCur$Exon.region.end..bp.
      indexStart  <-  gen_cord > transCur$Exon.region.start..bp.
      cur_rank_down <-  min(transCur$Exon.rank.in.transcript[indexEnd & indexStart])
      
      indexUp   <- transCur$Exon.rank.in.transcript == cur_rank_up
      indexDown <- transCur$Exon.rank.in.transcript == cur_rank_down 
      second_up_ss_strength   <- paste0("SA: ",transCur$acc_seq_MAXENT[indexUp])
      second_down_ss_strength  <- paste0("SD: ",transCur$donor_seq_HBS[indexDown])
    }
  }
  
  return(c(second_up_ss_strength, second_down_ss_strength))
}
