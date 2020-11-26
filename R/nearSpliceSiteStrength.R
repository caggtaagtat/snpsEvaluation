nearSpliceSiteStrength <- function(transCur, gen_cord, strand, location){
  
  #Convert strand annotation
  if(strand== "+") strand <- "1"
  if(strand== "-") strand <- "-1"
  
  ## If located within exons
  if(location == "exonic"){

    ## Find current Exon
    withinIndex <- gen_cord >= transCur$Exon.region.start..bp. & 
    gen_cord <= transCur$Exon.region.end..bp.
    cur_rank <-  transCur$Exon.rank.in.transcript[withinIndex]
      
    ## Access saved ss strengths
    if(cur_rank != 1 ) second_up_ss_strength   <- transCur$donor_seq_HBS[transCur$Exon.rank.in.transcript == cur_rank-1 ]
    if(cur_rank !=  max(transCur$Exon.rank.in.transcript)) second_down_ss_strength  <- transCur$acc_seq_MAXENT[transCur$Exon.rank.in.transcript == cur_rank+1 ]

  }
  
  ## If located within intron
  if(location == "intronic"){
    
    ## If located on plus strand
    if(strand==1){
      
      ## Get upstream exon rank 
      indexEnd    <-  gen_cord > transCur$Exon.region.end..bp. 
      indexStart  <-  gen_cord > transCur$Exon.region.start..bp.
      cur_rank_up <-  max(transCur$Exon.rank.in.transcript[indexEnd & indexStart])
      
      ## Get downstream exon rank 
      indexEnd      <-  gen_cord < transCur$Exon.region.end..bp.
      indexStart    <-  gen_cord < transCur$Exon.region.start..bp.
      cur_rank_down <-  min(transCur$Exon.rank.in.transcript[indexEnd & indexStart])
      
      indexUp   <- transCur$Exon.rank.in.transcript == cur_rank_up
      indexDown <- transCur$Exon.rank.in.transcript == cur_rank_down
      second_up_ss_strength   <- transCur$donor_seq_HBS[indexUp]
      second_down_ss_strength  <- transCur$acc_seq_MAXENT[indexDown]
      
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
      second_up_ss_strength  <- transCur$acc_seq_MAXENT[indexUp]
      second_down_ss_strength  <- transCur$donor_seq_HBS[indexDown]  
    }
  }
  
  return(c(second_up_ss_strength, second_down_ss_strength))
}
