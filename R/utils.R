#' Convert strands in VCF into BED
#'
#' This function read SV in bed format and standardise the strands
#'
#' @param bed SV in bed format
#' @return SV in bed format
#' @export
strands_standardisation <- function(bed){
  if(length(bed$ID_caller) == 0){bed$ID_caller = c(1:nrow(bed))}
  strand1 <- rep(NA, nrow(bed)); strand2 <- rep(NA, nrow(bed))
  if(sum(!is.na(bed$INFO_STRANDS))!=0){
    tmp <- data.frame(stringr::str_split_fixed(bed$INFO_STRANDS,",",2), bed$ID_caller)
    tmp$X1 <- gsub(":.*","", tmp$X1)
    tmp$X2 <- gsub(":.*","", tmp$X2)
    idx <- rep(1:nrow(tmp), rowSums(tmp[,c(1,2)]!=""))
    tmp2 <- tmp[idx,]

    tmp2[duplicated(tmp2$bed.ID_caller),]$X1 <- tmp2[duplicated(tmp2$bed.ID_caller),]$X2
    tmp3 <- tmp2[,-2]
    tmp3 <- cbind(tmp3, stringr::str_split_fixed(tmp3$X1,"",2))
    colnames(tmp3) <- c("strands","ID_caller","strand1","strand2")
  }else if(sum(!is.na(bed$INFO_CT)) != 0){
    tmp <- data.frame(INFO_CT = bed$INFO_CT, ID_caller = bed$ID_caller)
    strand1[tmp$INFO_CT=="3to5"] <- "+"; strand2[tmp$INFO_CT=="3to5"] <- "-"
    strand1[tmp$INFO_CT=="5to3"] <- "-"; strand2[tmp$INFO_CT=="5to3"] <- "+"
    strand1[tmp$INFO_CT=="3to3"] <- "+"; strand2[tmp$INFO_CT=="3to3"] <- "+"
    strand1[tmp$INFO_CT=="5to5"] <- "-"; strand2[tmp$INFO_CT=="5to5"] <- "-"
    strand1[tmp$INFO_CT=="NtoN"] <- NA; strand2[tmp$INFO_CT=="NtoN"] <- NA
    tmp3 <- cbind(tmp, strand1, strand2)
  }else if(sum(!is.na(bed$INFO_INV5)) !=0){
    tmp <- data.frame(INFO_INV5 = bed$INFO_INV5, INFO_INV3 = bed$INFO_INV3, ID_caller = bed$ID_caller)
    strand1[bed$INFO_INV5] <- "-";strand2[bed$INFO_INV5] <- "-"
    strand1[bed$INFO_INV3] <- "+";strand2[bed$INFO_INV3] <- "+"
    tmp3 <- cbind(tmp, strand1, strand2)
  }else{
    tmp3 <- cbind(bed, strand1, strand2)
  }

  bed2 <- bed[match(tmp3$ID_caller, bed$ID_caller),]
  bed2$strand1 <- tmp3$strand1
  bed2$strand2 <- tmp3$strand2

  if((sum(bed2$INFO_SVTYPE %in% c("DEL","DUP")) != 0) & (sum(is.na(bed2$strand1))!=0)){
    if(length(bed2[is.na(bed2$strand1) & bed2$INFO_SVTYPE == "DEL",]$strand1) != 0){
      bed2[is.na(bed2$strand1) & bed2$INFO_SVTYPE == "DEL",]$strand1 <- "+"
      bed2[is.na(bed2$strand2) & bed2$INFO_SVTYPE == "DEL",]$strand2 <- "-"
    }
    if(length(bed2[is.na(bed2$strand1) & bed2$INFO_SVTYPE == "DUP",]$strand1) != 0){
      bed2[is.na(bed2$strand1) & bed2$INFO_SVTYPE == "DUP",]$strand1 <- "-"
      bed2[is.na(bed2$strand2) & bed2$INFO_SVTYPE == "DUP",]$strand2 <- "+"
    }
  }
return(bed2)
}

#' Convert bed format to bedpe format
#'
#' This function read SV data in bed format and convert it into bedpe format
#'
#' @param bed SV data in bed format
#' @return SV data in bedpe format
#' @export
bed_to_bedpe <- function(bed){
  bed <- strands_standardisation(bed)

  ALT <- bed$ALT
  tmp <- gsub("\\:.*",'', ALT[grepl('\\[',ALT) | grepl('\\]',ALT)])
  tmp1 <- gsub(".*\\[",'', tmp)
  tmp2 <- gsub(".*\\]",'', tmp1)
  chrom1 = as.character(bed$CHROM)
  if(sum(grepl("chr", chrom1))!=nrow(bed)){
    chrom1[!grepl("chr", chrom1)] <- paste0("chr", chrom1[!grepl("chr", chrom1)])
  }

  chrom2 <- chrom1
  chrom2[grepl('\\[',ALT) | grepl('\\]',ALT)] <- tmp2

  if(sum(grepl("CHR", chrom2))!=0){
    chrom2[grepl("CHR", chrom2)] <- paste0("chr", substring(chrom2[grepl("CHR", chrom2)],4))
  }

  if(sum(grepl("chr", chrom2))!=nrow(bed)){
    chrom2[!grepl("chr", chrom2)] <- paste0("chr", chrom2[!grepl("chr", chrom2)])
  }

  tmp <- gsub(".*\\:",'', ALT[grepl('\\[',ALT) | grepl('\\]',ALT)])
  tmp1 <- gsub("\\[.*",'', tmp)
  tmp2 <- gsub("\\].*",'', tmp1)
  pos2 <- bed$INFO_END
  pos2[grepl('\\[',ALT) | grepl('\\]',ALT)] <- tmp2
  pos2[is.na(pos2)] <- bed$POS[is.na(pos2)] ###if no POS2 then use #CHROM in VCF

  strand1 <- bed$strand1; strand2 <- bed$strand2
  ###[p[t
  strand1[grepl('\\[', bed$ALT) & substr(bed$ALT,1,1) == "["] <- "-"
  strand2[grepl('\\[', bed$ALT) & substr(bed$ALT,1,1) == "["] <- "-"

  ###t[p[
  strand1[grepl('\\[', bed$ALT) & substr(bed$ALT,1,1) != "["] <- "+"
  strand2[grepl('\\[', bed$ALT) & substr(bed$ALT,1,1) != "["] <- "-"

  ###t]p]
  strand1[grepl(']', bed$ALT) & substr(bed$ALT,1,1) != "]"] <- "+"
  strand2[grepl(']', bed$ALT) & substr(bed$ALT,1,1) != "]"] <- "+"

  ###]p]t
  strand1[grepl(']', bed$ALT) & substr(bed$ALT,1,1) == "]"] <- "-"
  strand2[grepl(']', bed$ALT) & substr(bed$ALT,1,1) == "]"] <- "+"

  bedpe <- data.frame(chrom1=chrom1,
                      start1=bed$POS-1,
                      end1=bed$POS,
                      #pos1=bed$POS,
                      chrom2=paste0(tolower(substring(chrom2,1,3)), substring(chrom2,4,nchar(chrom2))),
                      start2 = as.numeric(pos2)-1,
                      end2=pos2,
                      #pos2=pos2,
                      strand1 = strand1,
                      strand2 = strand2,
                      ALT=ALT,
                      stringsAsFactors = FALSE)

  bedpe <- data.frame(bedpe, bed[,!(colnames(bed) %in% c("CHROM","POS","INFO_END","ALT","strand1","strand2"))])
  colnames(bedpe) <- c("chrom1", "start1","end1","chrom2","start2","end2","strand1","strand2","ALT",
                    # c("chrom1", "pos1","chrom2","pos2","strand1","strand2","ALT",
                       colnames(bed)[!(colnames(bed) %in% c("CHROM","POS","INFO_END","ALT","strand1","strand2"))])
  return(bedpe)
}

#' Prepare data for SV circos plot
#'
#' This function read SV data in bedpe format and convert to list of input data for circlize package
#'
#' @param bedpe SV data in bedpe format, for example output from simple_SVTYPE_classification
#' @return list of input data for circlize package
#' @export
prepare_SV_for_circos <- function(bedpe){
  bed1 <- data.frame(chr = bedpe$chrom1,
                     start = bedpe$start1,
                     end = bedpe$end1,
                     SVTYPE = bedpe$SVTYPE)
  bed2 <- data.frame(chr = bedpe$chrom2,
                     start = bedpe$start2,
                     end = bedpe$end2,
                     SVTYPE = bedpe$SVTYPE)
  bed3 <- bed1[bed1$SVTYPE!="TRA",]
  bed4 <- bed2[bed2$SVTYPE!="TRA",]

  bed5 <- bed1[bed1$SVTYPE =="TRA",]
  bed6 <- bed2[bed2$SVTYPE =="TRA",]
  return(list(bed3, bed4, bed5, bed6))
}

#' Prepare breakpoint to use for bedtools pairtopair
#'
#' This function read SV data and convert it to bedpe, which can be used in bedtools pairtopair
#'
#' @param SV_data VCF file or SV data in bed format
#' @param bkpt_T_caller the threshold of breakpoint difference
#' @param caller caller name or unique identifier of VCF
#' @return SV data in bedpe format
#' @export
standard_bedtool_prepare_bkpt <- function(SV_data, bkpt_T_callers,caller){
  caller_bedpe <- simple_SVTYPE_classification(SV_data, caller)
  bedpe <- data.frame(chrom1 = caller_bedpe$chrom1,
                      start1 = caller_bedpe$start1-(bkpt_T_callers/2),
                      end1 = caller_bedpe$end1+(bkpt_T_callers/2),
                      chrom2 = caller_bedpe$chrom2,
                      start2 = caller_bedpe$start2-(bkpt_T_callers/2),
                      end2 = caller_bedpe$end2+(bkpt_T_callers/2),
                      ID = paste0(caller_bedpe$ID,";",caller_bedpe$SVTYPE,";",caller_bedpe$FILTER))
  if (sum(bedpe$start1<0)!=0){bedpe[bedpe$start1<0,]$start1 <- 0}
  if (sum(bedpe$start2<0)!=0){bedpe[bedpe$start2<0,]$start2 <- 0}
  return(bedpe)
}


#' Filter SV types and breakpoint positions
#'
#' This function read bed format
#'
#' @param pairtopair file name
#' @param PASS_filter filtering of FILTER field in VCF of two calls
#' @param svtype_ignore whether ignore SV type for integration
#' @return data frame
#' @export
pairtopair_filter <- function(pairtopair,PASS_filter, svtype_ignore){
  if(PASS_filter=="both" & (!svtype_ignore)){
    pairtopair_filtered <- pairtopair[pairtopair$caller1_FILTER=="PASS" & pairtopair$caller2_FILTER=="PASS" &
                                        pairtopair$caller1_SVTYPE == pairtopair$caller2_SVTYPE ,]
  }else if(PASS_filter=="one"& (!svtype_ignore)){
    pairtopair_filtered <- pairtopair[(pairtopair$caller1_FILTER=="PASS" | pairtopair$caller2_FILTER=="PASS") &
                                        pairtopair$caller1_SVTYPE == pairtopair$caller2_SVTYPE,]
  }else{pairtopair_filtered <- pairtopair[pairtopair$caller1_SVTYPE == pairtopair$caller2_SVTYPE,]}
  return(pairtopair_filtered)
}

#' Integrate SV call sets
#'
#' This function read bed format
#'
#' @param sampleID sample ID
#' @param df_SV name of callers
#' @param vcf_files names of call sets
#' @param bkpt_T_callers threshold of breakpoint difference
#' @param SVCaller_names whether ignore SV type for integration
#' @param bedtools_dir directory of bedtools
#' @param DIR directory to write temporary files in
#' @return data frame
#' @export
SVCaller_union_intersect_generate <- function(sampleID, df_SV, vcf_files, bkpt_T_callers, SVCaller_names, bedtools_dir, DIR){
  bedpe <- standard_bedtool_prepare_bkpt(df_SV, bkpt_T_callers, caller = SVCaller_names[1])
  write.table(bedpe, paste0(DIR, sampleID, "/",SVCaller_names[1],"_sv.bedpe"), quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
  bedpe <- standard_bedtool_prepare_bkpt(vcf_files[2], bkpt_T_callers, caller = SVCaller_names[2])
  write.table(bedpe, paste0(DIR, sampleID, "/",SVCaller_names[2],"_sv.bedpe"), quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)

  cat(paste(sampleID, SVCaller_names[1], "SV data in bedpe write.\n"))
  cat(paste(sampleID, SVCaller_names[2], "SV data in bedpe write.\n"))

  system(paste(bedtools_dir,"pairtopair -a", paste0(DIR, sampleID, "/",SVCaller_names[1],"_sv.bedpe"),
               "-b", paste0(DIR, sampleID, "/",SVCaller_names[2],"_sv.bedpe"), ">",
               paste0(DIR, sampleID, "/",SVCaller_names[1],"_vs_",SVCaller_names[2],"_sv.bedpe")))

  if(file.info(paste0(DIR, sampleID, "/",SVCaller_names[1],"_vs_",SVCaller_names[2],"_sv.bedpe"))$size==0){
    pairtopair <- c()
  }else{
    pairtopair <- read.table(paste0(DIR, sampleID, "/",SVCaller_names[1],"_vs_",SVCaller_names[2],"_sv.bedpe"), header =FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
    pairtopair <- pairtopair[!duplicated(pairtopair),] ###Some duplicated rows, not sure why, so here remove duplicates
    pairtopair <- cbind(pairtopair[,c(1:6)],
                        stringr::str_split_fixed(pairtopair[,7],";",3)[,1],
                        stringr::str_split_fixed(pairtopair[,7],";",3)[,2],
                        stringr::str_split_fixed(pairtopair[,7],";",3)[,3],
                        pairtopair[,c(8:13)],
                        stringr::str_split_fixed(pairtopair[,14],";",3)[,1],
                        stringr::str_split_fixed(pairtopair[,14],";",3)[,2],
                        stringr::str_split_fixed(pairtopair[,14],";",3)[,3])

    colnames(pairtopair) <- c("caller1_chrom1","caller1_start1","caller1_end1",
                              "caller1_chrom2","caller1_start2","caller1_end2",
                              "caller1_ID","caller1_SVTYPE","caller1_FILTER",
                              "caller2_chrom1","caller2_start1","caller2_end1",
                              "caller2_chrom2","caller2_start2","caller2_end2",
                              "caller2_ID","caller2_SVTYPE","caller2_FILTER")

    main_chroms <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
                     "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                     "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
    pairtopair <- pairtopair[(pairtopair$caller1_chrom1 %in% main_chroms) &
                               (pairtopair$caller1_chrom2 %in% main_chroms) &
                               (pairtopair$caller2_chrom1 %in% main_chroms) &
                               (pairtopair$caller2_chrom2 %in% main_chroms),]
  }

  return(pairtopair)
}



#' Check the directory of bedtools
#'
#' This function check the binary of bedtools
#'
#' @param x "bedtools"
#' @return bedtools path
#' @export
Check_bedtools <- function(x = "bedtools"){
  # check if binary is in path
  cat(paste0('Checking path for ', x, '... ') );
  if (Sys.which(x) == "") {
    cat(paste0('FAIL\n') )
    #path <- Sys.getenv('PATH')
  }
  else {
    cat(paste0('PASS\n    ', Sys.which(x), '\n') )
  }
  return(Sys.which(x))
}


#' Define frequency of gene fusions
#'
#' Define frequency of gene fusions
#'
#' @param bed_geneAnnotated data frame
#' @return data frame of summarised gene fusions
#' @export
Summary_gene_fusions <- function(SV_geneAnnotated){
  if(!("sampleID" %in% colnames(SV_geneAnnotated))){
    sampleID <- "sample_1"
    df_All_gene_fusions <- cbind(sampleID, SV_geneAnnotated)
  }else{
    df_All_gene_fusions <- SV_geneAnnotated
  }

  df_All_gene_fusions$gene_fusions <- paste0(df_All_gene_fusions$pos1_overlap_gene,";",df_All_gene_fusions$pos2_overlap_gene)
  df_All_gene_fusions <- unique(df_All_gene_fusions) ###in the gene bed, there are same gene names but different gene id, so keep unique here
  df_All_gene_fusions2 <- df_All_gene_fusions[!(is.na(df_All_gene_fusions$pos1_overlap_gene)) & !(is.na(df_All_gene_fusions$pos2_overlap_gene)),]
  df_All_gene_fusions2 <- df_All_gene_fusions2[df_All_gene_fusions2$pos1_overlap_gene!= df_All_gene_fusions2$pos2_overlap_gene,]
  length(unique(paste0(df_All_gene_fusions2$sampleID, df_All_gene_fusions2$ID)))

  df_summary_gene_fusions <- data.frame(table(df_All_gene_fusions2$gene_fusions))

  tmp <- stringr::str_split_fixed(as.character(df_summary_gene_fusions$Var1), ";", 2)
  df_summary_gene_fusions <- cbind(data.frame(tmp), df_summary_gene_fusions)
  colnames(df_summary_gene_fusions) <- c("pos1_overlap_gene","pos2_overlap_gene","gene_fusions", "breakpoint_count")
  sample_count <- c()
  for(i in c(1: nrow(df_summary_gene_fusions))){
    tmp <- df_All_gene_fusions2[df_All_gene_fusions2$gene_fusions == df_summary_gene_fusions[i,]$gene_fusions,]
    sample_count <- c(sample_count, length(unique(tmp$sampleID)))
  }
  df_summary_gene_fusions$sample_count <- sample_count
  return(df_summary_gene_fusions)
}

