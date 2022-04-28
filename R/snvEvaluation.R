## Calculate various parameters per SNV
snvEvaluation <- function(startSNVdata=all_data, referenceDnaStringSet=referenceDnaStringSet,
                                      Transcript_cords_chr=Transcript_cords_chr){
  
    ## Summarize per SNP
    all_data$ID <- seq_len(nrow(all_data))
    snp <- data.frame(ID=unique(all_data$ID))
    snp$gene <- all_data$gene[match(snp$ID, all_data$ID)]
    snp$chromosome <- all_data$chromosome[match(snp$ID, all_data$ID)]
    snp$genome_position <- all_data$genome_position[match(snp$ID, all_data$ID)]
    snp$c.DNA <- all_data$c.DNA[match(snp$ID, all_data$ID)]
    snp$strand <-  all_data$strand[match(snp$ID, all_data$ID)]
    snp$ENSEMBL_transcriptID <- all_data$ENSEMBL_transcriptID[match(snp$ID, all_data$ID)]
    
    inputSNVs <- snp
    
    rm(all_data, new, snp)
    
    ## Insert collected information
    inputSNVs$SNP_ID <- paste(inputSNVs$chromosome, inputSNVs$genome_position)
    inputSNVs$genome_position <- as.numeric(inputSNVs$genome_position)
    
    #Convert strand annotation
    inputSNVs$strand[inputSNVs$strand == "+"]  <- "1"
    inputSNVs$strand[inputSNVs$strand == "-"]  <- "-1"
    
    ## Create results dataframe
    snvTable <- createOutputDF(nrow(inputSNVs))
    
    case <- 1
    for(case in 1:nrow(inputSNVs)){
        
        ## Get strand from  gene name  
        strand <-  inputSNVs$strand[case]
        snvTable$gene_strand[case] <- strand
        snvTable$gene[case] <- inputSNVs$gene[case]
        
        ## Save annotation 
        anno <- strsplit(as.character(inputSNVs$c.DNA[case]),"/")[[1]][[1]]
        snvTable$c.DNA[case] <- anno
        
        ## Select relevant transcript table
        transcript <- inputSNVs$ENSEMBL_transcriptID[case]
        snvTable$ENSEMBL_transcriptID[case] <- transcript
        transCur <- Transcript_cords_chr[Transcript_cords_chr$Transcript.stable.ID == transcript,]
        
        ## Get genomic coordiante and the reference and alternative NT
        output <- genomicPosAndNTs(anno, referenceDnaStringSet, transcript, transCur)
        
        ## Convert and save coordiantes and NT
        gen_cord <- strsplit(output, " ")[[1]][[1]]
        gen_cord <- as.numeric(gen_cord)
        snvTable$genome_position[case] <- gen_cord
        refNuc <- strsplit(output, " ")[[1]][[2]]
        altNuc <- strsplit(output, " ")[[1]][[3]]
        deletion_length <- as.numeric(strsplit(output, " ")[[1]][[4]] )
        insertion_length <- as.numeric(strsplit(output, " ")[[1]][[5]]) 
        
        snvTable$variation_type[case] <- "SNV"
        snvTable$alt[case] <- altNuc
        
        ## Define type of variation
        if(str_detect(anno, "del")) snvTable$variation_type[case] <- "DEL"
        if(str_detect(anno, "dup")) snvTable$variation_type[case] <- "DUP"
        if(str_detect(anno, "ins")) snvTable$variation_type[case] <- "INS"
        varType <- snvTable$variation_type[case]
        
        ## Save chromosome
        chromosome <- inputSNVs$chromosome[case]
        snvTable$chromosome[case] <- chromosome
        
        ## Get gerenal change in HZEI integral and genomic sequence
        quickContextOutput <- quickContext(chromosome, gen_cord, strand, 
                                           referenceDnaStringSet, varType, 
                                           deletion_length)
        
        snvTable$Delta_HZEI_integral[case] <- as.numeric(quickContextOutput[1])
        
        ## Calculate distance of variation to all transcript splice sites
        splicesites <- getSpliceSiteDist(transCur, gen_cord, strand)
        
        ## Enter Info about nearest SD
        SDinfo <- getInfoSDanno(splicesites, varType, altNuc, gen_cord,
                                chromosome, strand, referenceDnaStringSet,
                                deletion_length, insertion_length)
        
        ## Save gathered SD information
        snvTable$SNP_distance_to_SD_nt[case] <- SDinfo$`SNP distance to SD NR`
        snvTable$SD_HBS_diff[case]           <- SDinfo$`SD HBS delta`
        snvTable$SD_HBS_ref[case]            <- SDinfo$`SD HBS ref`
        snvTable$SD_HBS_alt[case]            <- SDinfo$`SD HBS alt`
        if(SDinfo$`SD HBS alt` == -999 &
           SDinfo$`SD HBS ref` != 999) snvTable$annoSD_GT_dinucleotide_affected[case]  <- "yes"
        SDcor <- SDinfo$SDcor
        snvTable$SNV_localized_relative_to_SD[case] <- SDinfo$SDloc
        rm(SDinfo)
        
        ## Enter distance to next SA
        SAinfo <- getInfoSAanno(splicesites, varType, altNuc, gen_cord,
                                chromosome, strand, referenceDnaStringSet,
                                deletion_length, insertion_length)
        
        ## Save gathered SA information
        snvTable$SNP_distance_to_SA_nt[case] <- SAinfo$`SNP distance to SA NR`
        snvTable$SA_MaxEnt_ref[case]            <- SAinfo$`SA MES ref`
        snvTable$SA_MaxEnt_alt[case]            <- SAinfo$`SA MES alt`
        if(SAinfo$`SA MES delta` == 9999) snvTable$annoSA_AG_dinucleotide_affected[case]  <- "yes"
        snvTable$Delta_MaxEnt_anno_SA[case]  <- SAinfo$`SA MES alt`-SAinfo$`SA MES ref`
        SAcor <- SAinfo$SAcor
        snvTable$SNV_localized_relative_to_SA[case] <- SAinfo$SAloc
        rm(SAinfo)
        
        ## Access changes in HBS or MES of overlapping splice sites
        overlappingInfo <- getAffectedSpliceSites(quickContextOutput[2],
                                                  quickContextOutput[3],
                                                  quickContextOutput[4],
                                                  quickContextOutput[5],
                                                  SDcor, SAcor)
        
        ## Enter information
        snvTable$HBS_alt_overlapping_SDs[case]  <- overlappingInfo$sd_max_alt
        snvTable$HBS_diff_overlapping_SDs[case] <- overlappingInfo$sd_max_diff
        snvTable$MaxEnt_alt_overlapping_SAs[case] <- overlappingInfo$sa_max_alt
        snvTable$MaxEnt_diff_overlapping_SAs[case]  <- overlappingInfo$sa_max_diff
        snvTable$Maxent_alt_overlapping_SD_GC_sites[case] <- overlappingInfo$sd_max_altGC
        snvTable$Maxent_diff_overlapping_SD_GC_sites[case]  <- overlappingInfo$sd_max_diffGC
        
        ## Get SRE support of potential GT-,GC- AG-sites
        GTcor <- overlappingInfo$sd_coord
        GCcor <- overlappingInfo$sd_coordGC
        AGcor <- overlappingInfo$sa_coord
        
        ## Asses SRE support for annotated nearest GT-, GC- and SA-site
        GTcorSREsupport <- assesSREsupportSD(GTcor, varType, altNuc, gen_cord, chromosome, 
                                             strand, referenceDnaStringSet, posToSD)
        GCcorSREsupport <- assesSREsupportSD(GCcor, varType, altNuc, gen_cord, chromosome, 
                                             strand, referenceDnaStringSet, posToSD)
        AGcorSREsupport <- assesSREsupportSA(AGcor, varType, altNuc, gen_cord, chromosome, 
                                             strand, referenceDnaStringSet, posToSD)
        
        GTcorSRE <- GTcorSREsupport[2]
        GCcorSRE <- GCcorSREsupport[2]
        AGcorSRE <- AGcorSREsupport[2]
        
        ## Enter the information
        snvTable$SRE_support_alt_overlapping_SD[case] <- GTcorSRE
        snvTable$SRE_support_alt_overlapping_GCsite[case] <- GCcorSRE
        snvTable$SRE_support_alt_overlapping_SA[case] <- AGcorSRE
        
        ## Compare with the annotated SD and SA
        altHBS <- snvTable$HBS_alt_overlapping_SDs[case]
        snvTable$HBS_diff_overlapping_SDs_vs_annotated[case] <- altHBS - snvTable$SD_HBS_ref[case]
        
        altMES <- snvTable$MaxEnt_alt_overlapping_SAs[case]
        snvTable$MaxEnt_diff_overlapping_SA_vs_annotated[case] <- altMES - snvTable$SA_MaxEnt_ref[case]
        
        
        ## Define relative positin of variation to next annotated splice sites
        posToSD <- snvTable$SNV_localized_relative_to_SD[case]
        posToSA <- snvTable$SNV_localized_relative_to_SA[case]
        
        ## Asses SRE support for annotated nearest SD
        sdSREsupport <- assesSREsupportSD(SDcor, varType, altNuc, gen_cord, chromosome, 
                                          strand, referenceDnaStringSet, posToSD)
        
        ## save gathered SRE info about annotated SD                 
        snvTable$SD_SRE_support[case] <- sdSREsupport[1]                                  
        snvTable$SD_SRE_support_alternative[case] <- sdSREsupport[2]  
        snvTable$Delta_SD_SRE_support[case] <- sdSREsupport[2] - sdSREsupport[1]
        snvTable$Delta_SD_SRE_support_percent[case] <- ((sdSREsupport[2] - sdSREsupport[1])/sdSREsupport[1])*100
        snvTable$SD_sequence_surrounding[case] <- sdSREsupport[3]
        snvTable$SD_sequence_surrounding_with_SNP[case] <- sdSREsupport[4]
        
        ## Asses SRE support for annotated nearest SA
        saSREsupport <- assesSREsupportSA(SAcor, varType, altNuc, gen_cord, chromosome, 
                                          strand, referenceDnaStringSet, posToSA)
        
        ## save gathered SRE info about annotated SA                
        snvTable$SA_SRE_support[case] <- saSREsupport[1]                                  
        snvTable$SA_SRE_support_alternative[case] <- saSREsupport[2]  
        snvTable$Delta_SA_SRE_support[case] <- saSREsupport[2] - saSREsupport[1]
        snvTable$Delta_SA_SRE_support_percent[case] <- ((saSREsupport[2] - saSREsupport[1])/saSREsupport[1])*100
        snvTable$SA_sequence_surrounding[case] <- saSREsupport[3]
        snvTable$SA_sequence_surrounding_with_SNP[case] <- saSREsupport[4]
        
        ## Define transcript context
        snvTable$transcript_context[case] <- "intronic"
        
        ## Test where the variation lies in the gene context
        for(inn in 1:nrow(transCur)){
            
            if(gen_cord %in% c(transCur$Exon.region.start..bp.[inn]:
                               transCur$Exon.region.end..bp.[inn])) snvTable$transcript_context[case] <- "exonic"
            
        }
        
        ## Retrieve Length and sequence of nearest exon
        transCur$exon_length <- abs(transCur$Exon.region.end..bp.-transCur$Exon.region.start..bp.)
        indexStart <- transCur$Exon.region.start..bp. == splicesites$coordinate[1] 
        indexEnd <- transCur$Exon.region.end..bp. == splicesites$coordinate[1]
        snvTable$length_of_next_exon[case] <- transCur$exon_length[indexStart | indexEnd]+1
        exonlength <- transCur$exon_length[indexStart | indexEnd]
        
        ## Retrive sequence
        upborder <-  transCur$Exon.region.start..bp.[indexStart | indexEnd]
        if(strand == -1) upborder <-  transCur$Exon.region.end..bp.[indexStart | indexEnd]
        
        snvTable$Exon_sequence[case] <- getReferenceSequence(chromosome = chromosome, 
                                                             indexCoordinate = upborder,
                                                             upRange = 0, downRange = exonlength, 
                                                             strand = strand, referenceDnaStringSet)
        
        context <- snvTable$transcript_context[case]
        transCur$acc_seq_MAXENT <- as.numeric(as.character(transCur$acc_seq_MAXENT))
        
        ssInfo <- nearSpliceSiteStrength(transCur, gen_cord, 
                                         strand, context)
        
        snvTable$second_up_ss_strength[case]    <- ssInfo[1]
        snvTable$second_down_ss_strength[case]  <- ssInfo[2]
        
        ## Retrieve SA surrounding sequences
        SAcontext <- getSAcontext(SAcor, varType, altNuc, gen_cord, chromosome, 
                                  strand, referenceDnaStringSet, posToSA)
        
        snvTable$SA_sequence_surrounding[case] <- SAcontext[1]
        snvTable$SA_sequence_surrounding_with_SNP[case] <- SAcontext[2]
        
        ## Retrieve SD surrounding sequences 
        SDcontext <- getSDcontext(SDcor, varType, altNuc, gen_cord, chromosome, 
                                  strand, referenceDnaStringSet, posToSD)
        
        snvTable$SD_sequence_surrounding[case] <- SDcontext[1]
        snvTable$SD_sequence_surrounding_with_SNP[case] <- SDcontext[2]
        
        ## Now mark potentially interesting variations which
        ## substantially alter splicing elements
        SDdown <- snvTable$annoSD_GT_dinucleotide_affected[case]
        SAdown <- snvTable$annoSA_AG_dinucleotide_affected[case]
        
        ## Values checking for new/strengthened alternative GT SD sites
        HBSaltOver <- snvTable$HBS_alt_overlapping_SDs[case]
        HBSaltOverDiff <- snvTable$HBS_diff_overlapping_SDs[case]
        HBSaltOverDiffAno <- snvTable$HBS_diff_overlapping_SDs_vs_annotated[case]
        sdHBS <- snvTable$SD_HBS_ref[case]
        sdHBSalt <- snvTable$SD_HBS_alt[case]
        saMES <- snvTable$SA_MaxEnt_ref[case]
        saMESalt <- snvTable$SA_MaxEnt_alt[case]
        
        ## Values checking for new/strengthened alternative AG SA sites
        MESaltOver <- snvTable$MaxEnt_alt_overlapping_SAs[case]
        MESaltOverDiff  <- snvTable$MaxEnt_diff_overlapping_SAs[case]
        MESaltOverDiffAno  <- snvTable$MaxEnt_diff_overlapping_SA_vs_annotated[case]
        
        ## SD SRE changes drastically
        sdSRE <- snvTable$SD_SRE_support[case]
        sdSREalt <- snvTable$SD_SRE_support_alternative[case]
        sdSREdiffperc <- snvTable$Delta_SD_SRE_support_percent[case]
        saSRE <- snvTable$SA_SRE_support[case]
        saSREalt <- snvTable$SA_SRE_support_alternative[case]
        saSREdiffperc <- snvTable$Delta_SA_SRE_support_percent[case]
        
        ## Distance to splice sites
        dist_to_SD  <- snvTable$SNP_distance_to_SD_nt[case]
        dist_to_SA  <- snvTable$SNP_distance_to_SD_nt[case]
        
        impactReport <- markInterestingVariations(SDdown, SDcor, SAdown, SAcor, sdHBS, saMES, 
                                                  HBSaltOver, HBSaltOverDiff, HBSaltOverDiffAno,
                                                  MESaltOver, MESaltOverDiff, MESaltOverDiffAno,
                                                  posToSD, posToSA, sdHBSalt, saMESalt, sdSREdiffperc,
                                                  sdSRE, sdSREalt, saSRE, saSREalt, saSREdiffperc,
                                                  dist_to_SD, dist_to_SA, GTcorSRE, GCcorSRE,
                                                  AGcorSRE )
        
        ## Save impac message
        snvTable$Potentially_interesting[case] <- impactReport
        snvTable$Next_SA_coordinate[case] <- SAcor
        snvTable$Next_SD_coordinate[case] <- SDcor
        
    }
    
    
    ## Transform numeric values to character and replace "." for excel
    snvTable <- convertTableForExcel(snvTable)
    snvTableDocu <- createOutputDocu()
  
  ## Return dataframe at the end
  return(snvTable)
}
