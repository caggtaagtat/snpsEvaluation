createOutputDF <- function(rowsDF){
  
  snvTable <- data.frame( gene="",
                          c.DNA="",
                          category="",
                          Delta_HZEI_integral=0,
                          Delta_MaxEnt_anno_SA=0,
                          annoSA_AG_dinucleotide_affected="no",
                          SD_HBS_diff=0,
                          annoSD_GT_dinucleotide_affected="no",
                          SNP_distance_to_SA_nt=0,
                          SNP_distance_to_SD_nt=0,
                          transcript_context="",
                          SA_SRE_support=0,
                          SA_SRE_support_alternative=0,
                          Delta_SA_SRE_support=0,
                          Delta_SA_SRE_support_percent=0,
                          SD_SRE_support=0,
                          SD_SRE_support_alternative=0,
                          Delta_SD_SRE_support=0,
                          Delta_SD_SRE_support_percent=0,
                          HBS_alt_overlapping_SDs=0,
                          HBS_diff_overlapping_SDs=0,
                          HBS_diff_overlapping_SDs_vs_annotated=0,
                          SRE_support_alt_overlapping_SD=0,
                          MaxEnt_alt_overlapping_SAs=0,
                          MaxEnt_diff_overlapping_SAs=0,
                          MaxEnt_diff_overlapping_SA_vs_annotated=0,
                          SRE_support_alt_overlapping_SA=0,
                          Maxent_alt_overlapping_SD_GC_sites=0,
                          Maxent_diff_overlapping_SD_GC_sites=0,
                          SRE_support_alt_overlapping_GCsite=0,
                          SD_HBS_ref=0,
                          SD_HBS_alt=0,
                          SA_MaxEnt_ref=0,
                          SA_MaxEnt_alt=0,
                          SNV_localized_relative_to_SD="",
                          SNV_localized_relative_to_SA="",
                          length_of_next_exon=0,
                          variation_type="",
                          second_up_ss_strength=0,
                          second_down_ss_strength=0,
                          gene_strand=0,
                          chromosome="",
                          ENSEMBL_transcriptID="",
                          genome_position=0,
                          alt="",
                          Exon_sequence="",
                          SD_sequence_surrounding="",
                          SD_sequence_surrounding_with_SNP="",
                          SA_sequence_surrounding="",
                          SA_sequence_surrounding_with_SNP="",
                          category_details="",
                          Next_SA_coordinate=0,
                          Next_SD_coordinate=0,
                          ID=seq_len(rowsDF))
  
  return(snvTable)
}
