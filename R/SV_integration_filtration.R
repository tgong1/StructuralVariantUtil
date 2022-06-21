#' Convert VCF format to bed format
#'
#' This function read vcf file and convert it to bed format
#'
#' @param vcf_file names of vcf file
#' @return data frame
#' @export
vcf_to_bed <- function(vcf_file){
  vcf <- VariantAnnotation::readVcf(vcf_file)

  fixed_df <- vcf@fixed
  gr <- vcf@rowRanges
  info <- vcf@info

  for(field in c("SVTYPE","SVLEN","END", "STRANDS","CT","INV5","INV3","MATEID")){
    tmp <- eval(parse(text=paste0("info$", field)))
    if (length(tmp) == 0){
      assign(paste0("INFO_", field), NA)
    }else{
      idx <- !(sapply(tmp, length))
      tmp[idx] <- NA
      if(field == "STRANDS"){
        assign(paste0("INFO_", field), sapply(tmp, paste, collapse = ","))
      }else{
        assign(paste0("INFO_", field), unlist(tmp))
      }
    }
  }

  bed <- data.frame(CHROM = gr@seqnames,
                    POS = gr@ranges@start,
                    ID_caller = gr@ranges@NAMES,
                    REF = as.character(fixed_df$REF),
                    ALT = data.frame(fixed_df$ALT)$value,
                    QUAL = fixed_df$QUAL,
                    FILTER = fixed_df$FILTER,
                    INFO_END,
                    INFO_SVTYPE,
                    INFO_SVLEN,
                    INFO_STRANDS,
                    INFO_CT = INFO_CT,
                    INFO_INV5 = INFO_INV5,
                    INFO_INV3 = INFO_INV3,
                    INFO_MATEID_caller = INFO_MATEID,
                    stringsAsFactors = FALSE
  )

  return(bed)
}

#' Classify SV types based on VCF
#'
#' This function read bed format
#'
#' @param vcf_file names of vcf file
#' @param caller_name name of caller to be used in ID
#' @return data frame
#' @export
simple_SVTYPE_classification <- function(vcf_file, caller_name){
  bed <- vcf_to_bed(vcf_file)
  bedpe <- bed_to_bedpe(bed)
  if(length(bedpe$ID_caller)!=0){
    bedpe <- bedpe[is.na(match(bedpe$ID_caller, bedpe$INFO_MATEID_caller)) | ### either don't have mate (i.e. not BND)
                     (c(1:nrow(bedpe)) < match(bedpe$ID_caller, bedpe$INFO_MATEID_caller)),] ###  or present first as BND
  }

  SV_index <- c(1:nrow(bedpe))
  event_index <- SV_index

  SVTYPE <- bedpe$INFO_SVTYPE
  SVTYPE[bedpe$INFO_SVTYPE == "INS"] <- "INS"
  SVTYPE[(bedpe$strand1 == "+") & (bedpe$strand2 == "-") &
           (bedpe$chrom1 == bedpe$chrom2) & (bedpe$INFO_SVTYPE != "INS")] <- "DEL"
  SVTYPE[(bedpe$strand1 == "-") & (bedpe$strand2 == "+") &
           (bedpe$chrom1 == bedpe$chrom2) & (bedpe$INFO_SVTYPE != "INS")] <- "DUP"

  SVTYPE[(bedpe$strand1 == bedpe$strand2) &
           (bedpe$chrom1 == bedpe$chrom2) &
           (bedpe$INFO_SVTYPE != "INS")] <- "INV"
  SVTYPE[(bedpe$chrom1 != bedpe$chrom2)] <- "TRA"
  SVTYPE[is.na(bedpe$strand1)] <- "INS"


  standard_bedpe <- data.frame(chrom1 = as.character(bedpe$chrom1),
                               pos1 = as.integer(bedpe$pos1),
                               chrom2 = as.character(bedpe$chrom2),
                               pos2 = as.integer(bedpe$pos2),
                               SVTYPE = SVTYPE,
                               strand1 = bedpe$strand1,
                               strand2 = bedpe$strand2,
                               ID = paste0(caller_name,"_",SV_index,"_1_",event_index),
                               ID_mate = paste0(caller_name,"_",SV_index,"_2_",event_index),
                               #ALT = bedpe$ALT,
                               stringsAsFactors = FALSE)

  bedpe_SVTYPE_classified <- data.frame(standard_bedpe, bedpe[,!(colnames(bedpe) %in% c("chrom1","chrom2","pos1","pos2","strand1","strand2"))])
  colnames(bedpe_SVTYPE_classified) <- c("chrom1", "pos1","chrom2","pos2","SVTYPE",
                                         "strand1","strand2",
                                         "ID","ID_mate",
                                         colnames(bedpe)[!(colnames(bedpe) %in% c("chrom1","chrom2","pos1","pos2","strand1","strand2"))])
  SVTYPE_stat <- SVTYPE_stat_generate(bedpe_SVTYPE_classified)
  rownames(SVTYPE_stat) <- vcf_file
  return(list(bedpe_SVTYPE_classified, SVTYPE_stat))
}

#' Integrate SV call sets and write output
#'
#' This function read bed format
#'
#' @param sampleID sample ID
#' @param SVCaller_name name of callers
#' @param vcf_list list of VCF files from different callers
#' @param bkpt_T_callers threshold of breakpoint difference
#' @param SVTYPE_ignore whether ignore SV type for integration
#' @param bedtools_dir directory of bedtools
#' @return data frame
#' @export
SV_integration <- function(sampleID, SVCaller_name, vcf_list, bkpt_T_callers, SVTYPE_ignore, bedtools_dir){
  BND_diff <- 2000
  directory <- "./"
  sub_directory <- paste0("./", paste0(SVCaller_name,collapse = ""))
  dir.create(sub_directory)
  SVTYPE_ignore_text <- ifelse(SVTYPE_ignore, "SVTYPE_ignore", "SVTYPE_same")

  for(i in c(1:length(SVCaller_name))){
    results <- simple_SVTYPE_classification(vcf_file = vcf_list[i], caller_name = SVCaller_name[i])
    tmp <- results[[1]]
    tmp <- tmp[tmp$chrom1 %in% paste0("chr", c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y")),]
    assign(paste0(SVCaller_name[i],"_standard_bedpe"), tmp)
  }

  SVCaller_bed_union <- SVCaller_union_intersect_generate(sampleID, SVCaller_name, SVCaller_bed_name, BND_diff, bkpt_T_callers, SVTYPE_ignore, bedtools_dir)
  index <- which(colnames(SVCaller_bed_union) %in% SVCaller_name)
  SVCaller_bed_union <- SVCaller_bed_union[,c(1:9, index)]
  write.table(SVCaller_bed_union,
              file = paste0(directory,"/",sub_directory,"/",
                            sampleID, "_", paste0(SVCaller_name,collapse = "_"),
                            "_combine_all_",bkpt_T_callers,"bp","_",SVTYPE_ignore_text,".bed"),
              row.names = FALSE,col.names = TRUE, quote = FALSE, sep = "\t")
  #head(SVCaller_bed_union)
  #SVCaller_bed_intersection <- SVCaller_bed_union[(!is.na(SVCaller_bed_union$ID_overlap_Caller)),]
  #write.table(SVCaller_bed_intersection,
  #            file = paste0(directory,"/",sub_directory,"/",
  #                          sampleID, "_", paste0(SVCaller_name,collapse = "_"),
  #                         "_combine_intersection_",bkpt_T_callers,"bp","_",SVTYPE_ignore_text,".bed"),
  #            row.names = FALSE,col.names = TRUE, quote = FALSE, sep = "\t")
  return(SVCaller_bed_union)
}
