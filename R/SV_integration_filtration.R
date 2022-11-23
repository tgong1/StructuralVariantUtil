#' Convert VCF format to bed format
#'
#' This function read vcf file and convert it to bed format
#'
#' @param vcf_file names of vcf file
#' @return data frame
#' @export
vcf_to_dataframe <- function(vcf_file){
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
#' @param SV_data SV VCF path or SV in data.frame
#' @param caller_name name of caller to be used in ID
#' @return data frame
#' @export
simple_SVTYPE_classification <- function(SV_data, caller_name="caller1"){
  if(is.data.frame(SV_data)){
    df <- SV_data
    cat(paste0("Read SV data called from ", caller_name, " in R data.frame format.\n"))
  }else{
    vcf_file <- SV_data
    df <- vcf_to_dataframe(vcf_file)
    cat(paste0("Read SV data called from ", caller_name, " in VCF format.\n"))
  }

  bedpe <- bed_to_bedpe(df)
  if(length(bedpe$ID_caller)!=0){
    bedpe <- bedpe[is.na(match(bedpe$ID_caller, bedpe$INFO_MATEID_caller)) | ### either don't have mate (i.e. not BND)
                     (c(1:nrow(bedpe)) < match(bedpe$ID_caller, bedpe$INFO_MATEID_caller)),] ###  or present first as BND
  }

  SV_index <- c(1:nrow(bedpe))
  event_index <- SV_index

  SVTYPE <- bedpe$INFO_SVTYPE
  SVTYPE[bedpe$INFO_SVTYPE == "INS"] <- "INS"
  bedpe$end2 <- as.numeric(bedpe$end2)
  SVTYPE[(bedpe$strand1 == "+") & (bedpe$strand2 == "-") &
           (bedpe$chrom1 == bedpe$chrom2) & (bedpe$INFO_SVTYPE != "INS") & (bedpe$end1 < bedpe$end2)] <- "DEL"
  SVTYPE[(bedpe$strand1 == "-") & (bedpe$strand2 == "+") &
           (bedpe$chrom1 == bedpe$chrom2) & (bedpe$INFO_SVTYPE != "INS") & (bedpe$end1 > bedpe$end2)] <- "DEL"

  SVTYPE[(bedpe$strand1 == "-") & (bedpe$strand2 == "+") &
           (bedpe$chrom1 == bedpe$chrom2) & (bedpe$INFO_SVTYPE != "INS")& (bedpe$end1 < bedpe$end2)] <- "DUP"
  SVTYPE[(bedpe$strand1 == "+") & (bedpe$strand2 == "-") &
           (bedpe$chrom1 == bedpe$chrom2) & (bedpe$INFO_SVTYPE != "INS")& (bedpe$end1 > bedpe$end2)] <- "DUP"

  SVTYPE[(bedpe$strand1 == bedpe$strand2) &
           (bedpe$chrom1 == bedpe$chrom2) &
           (bedpe$INFO_SVTYPE != "INS")] <- "INV"
  SVTYPE[(bedpe$chrom1 != bedpe$chrom2)] <- "TRA"
  SVTYPE[is.na(bedpe$strand1)] <- "INS"


  standard_bedpe <- data.frame(chrom1 = as.character(bedpe$chrom1),
                               #pos1 = as.integer(bedpe$pos1),
                               start1 = as.integer(bedpe$start1),
                               end1 = as.integer(bedpe$end1),
                               chrom2 = as.character(bedpe$chrom2),
                               #pos2 = as.integer(bedpe$pos2),
                               start2 = as.integer(bedpe$start2),
                               end2 = as.integer(bedpe$end2),
                               SVTYPE = SVTYPE,
                               strand1 = bedpe$strand1,
                               strand2 = bedpe$strand2,
                               ID = paste0(caller_name,"_",SV_index,"_1_",event_index),
                               ID_mate = paste0(caller_name,"_",SV_index,"_2_",event_index),
                               #ALT = bedpe$ALT,
                               stringsAsFactors = FALSE)

  # if(sum(abs(bedpe_SVTYPE_classified$pos1 - bedpe_SVTYPE_classified$pos2)<2 &
  #        bedpe_SVTYPE_classified$chrom1 == bedpe_SVTYPE_classified$chrom2)!=0){
  #   bedpe_SVTYPE_classified[abs(bedpe_SVTYPE_classified$pos1 - bedpe_SVTYPE_classified$pos2)<2 &
  #                             bedpe_SVTYPE_classified$chrom1 == bedpe_SVTYPE_classified$chrom2,]$SVTYPE <- "INS"
  # }
  #
  # bedpe_SVTYPE_classified <- data.frame(standard_bedpe,
  #                                       bedpe[,!(colnames(bedpe) %in% c("chrom1","chrom2","pos1","pos2","strand1","strand2"))])
  # colnames(bedpe_SVTYPE_classified) <- c("chrom1", "pos1","chrom2","pos2","SVTYPE",
  #                                        "strand1","strand2",
  #                                        "ID","ID_mate",
  #                                        colnames(bedpe)[!(colnames(bedpe) %in% c("chrom1","chrom2","pos1","pos2","strand1","strand2"))])

  if(sum(abs(standard_bedpe$end1 - standard_bedpe$end2)<2 &
         standard_bedpe$chrom1 == standard_bedpe$chrom2)!=0){
    standard_bedpe[abs(standard_bedpe$end1 - standard_bedpe$end2)<2 &
                     standard_bedpe$chrom1 == standard_bedpe$chrom2,]$SVTYPE <- "INS"
  }

  bedpe_SVTYPE_classified <- data.frame(standard_bedpe,
                                        bedpe[,!(colnames(bedpe) %in% c("chrom1","chrom2",
                                                                        "start1","start2",
                                                                        "end1","end2",
                                                                        "strand1","strand2"))])
  colnames(bedpe_SVTYPE_classified) <- c("chrom1", "start1", "end1",
                                         "chrom2","start2", "end2",
                                         "SVTYPE",
                                         "strand1","strand2",
                                         "ID","ID_mate",
                                         colnames(bedpe)[!(colnames(bedpe) %in%  c("chrom1","chrom2",
                                                                                   "start1","start2",
                                                                                   "end1","end2",
                                                                                   "strand1","strand2"))])
  return(bedpe_SVTYPE_classified)
}

#' Integrate SV call sets and write output
#'
#' This function read bed format
#'

#' @param vcf_files list of VCF files from different callers
#' @param SVCaller_names names of callers
#' @param sampleID unique identifier of sample, default as "sample_1".
#' @param bkpt_T_callers threshold of breakpoint difference. Default as 100 bp
#' @param PASS_filter filtering based on FILTER field of two calls in VCF: "both" to require both calls with "PASS", one to require one of the two calls with "PASS", "none" to ignore this filtering. Default is "both".
#' @param svtype_ignore whether ignore SV type for integration. Default as FALSE.
#' @param bedtools_dir directory of bedtools
#' @return data frame
#' @export
SV_integration <- function(vcf_files, SVCaller_names, sampleID = "sample_1", bkpt_T_callers =100, PASS_filter="both", svtype_ignore =FALSE, bedtools_dir=NULL){
  if(is.null(bedtools_dir)){bedtools_dir <- Check_bedtools(x = "bedtools")}else{cat(paste0("Provided path for bedtools ... \n", bedtools_dir,"\n"))}
  if(is.null(bedtools_dir) | bedtools_dir == ""){cat(paste0("ERROR: Please provide the bedtools path.\n"))}else{
    df_SV <- do.call("rbind", lapply(vcf_files, vcf_to_dataframe))
    for(i in c(1:length(vcf_files))){
      vcf_files_tmp <- c("all",vcf_files[i])
      SVCaller_names_tmp <- c("all", SVCaller_names[i])
      caller1_bedpe_tmp <- SV_integration_pairwise(sampleID, df_SV, vcf_files_tmp, bkpt_T_callers, SVCaller_names_tmp, PASS_filter, svtype_ignore,bedtools_dir)
      if(i == 1){
        caller1_bedpe <- caller1_bedpe_tmp
      }else{
        caller1_bedpe <- cbind(caller1_bedpe, caller1_bedpe_tmp[,24])
        colnames(caller1_bedpe)[ncol(caller1_bedpe)] <- paste0(SVCaller_names[i],"_ID")
      }
    }
    union_ID <- caller1_bedpe[, (ncol(caller1_bedpe)-length(SVCaller_names)+1):ncol(caller1_bedpe)]
    integrated_bedpe <- caller1_bedpe[!duplicated(union_ID),]
    union_ID <- integrated_bedpe[, (ncol(integrated_bedpe)-length(SVCaller_names)+1):ncol(integrated_bedpe)]
    for(i in c(1: length(SVCaller_names))){
      integrated_bedpe <- integrated_bedpe[!(duplicated(union_ID[,i]) & (!is.na(union_ID[,i]))),]
    }
    return(integrated_bedpe)
  }
}
#' Integrate SV call sets and write output
#'
#' This function read bed format
#'
#' @param sampleID unique identifier of sample
#' @param df_SV name of callers
#' @param vcf_files list of VCF files from different callers
#' @param bkpt_T_callers threshold of breakpoint difference
#' @param SVCaller_names names of callers
#' @param PASS_filter filtering based on FILTER field of two calls in VCF: "both" to require both calls with "PASS", one to require one of the two calls with "PASS", "none" to ignore this filtering. Default is "both".
#' @param svtype_ignore whether ignore SV type for integration
#' @param bedtools_dir directory of bedtools
#' @return data frame
#' @export
SV_integration_pairwise <- function(sampleID, df_SV, vcf_files, bkpt_T_callers, SVCaller_names, PASS_filter, svtype_ignore,bedtools_dir){
  DIR <- "./"
  dir.create(paste0(DIR, sampleID))
  pairtopair <- SVCaller_union_intersect_generate(sampleID, df_SV, vcf_files, bkpt_T_callers, SVCaller_names, bedtools_dir, DIR)

  caller1_bedpe <- StructuralVariantUtil::simple_SVTYPE_classification(df_SV, caller= SVCaller_names[1])
  main_chroms <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
                   "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                   "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
  caller1_bedpe <- caller1_bedpe[(caller1_bedpe$chrom1 %in% main_chroms) &
                                   (caller1_bedpe$chrom2 %in% main_chroms),]

  if(length(pairtopair)==0){
    caller1_bedpe$caller2_ID <- NA
  }else{
    pairtopair_filtered <- pairtopair_filter(pairtopair,PASS_filter, svtype_ignore)
    caller1_bedpe_filtered <- caller1_bedpe[caller1_bedpe$ID %in% pairtopair_filtered$caller1_ID, ]

    caller1_ID <- unique(pairtopair_filtered$caller1_ID)
    caller2_ID <- c()
    for(i in c(1:length(unique(pairtopair_filtered$caller1_ID)))){
      caller2_ID <- c(caller2_ID,
                      paste(pairtopair_filtered[pairtopair_filtered$caller1_ID == caller1_ID[i],]$caller2_ID, collapse = ","))
    }
    caller1_bedpe$caller2_ID <- NA
    caller1_bedpe$caller2_ID[match(caller1_ID, caller1_bedpe$ID)] <- caller2_ID
  }
  colnames(caller1_bedpe)[which(colnames(caller1_bedpe)=="caller2_ID")] <- paste0(SVCaller_names[2],"_ID")
  return(caller1_bedpe)
}


