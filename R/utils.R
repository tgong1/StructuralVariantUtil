#' Convert strands into standard form
#'
#' This function read bed format
#'
#' @param bed dataframe
#' @return data frame
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
  tmp <- data.frame(bed$INFO_CT, ID_caller = bed$ID_caller)
  strand1[tmp$INFO_CT=="3to5"] <- "+"; strand2[tmp$INFO_CT=="3to5"] <- "-"
  strand1[tmp$INFO_CT=="5to3"] <- "-"; strand2[tmp$INFO_CT=="5to3"] <- "+"
  strand1[tmp$INFO_CT=="3to3"] <- "+"; strand2[tmp$INFO_CT=="3to3"] <- "+"
  strand1[tmp$INFO_CT=="5to5"] <- "-"; strand2[tmp$INFO_CT=="5to5"] <- "-"
  strand1[tmp$INFO_CT=="NtoN"] <- NA; strand2[tmp$INFO_CT=="NtoN"] <- NA
  tmp3 <- cbind(tmp, strand1, strand2)
}else if(sum(!is.na(bed$INFO_INV5)) !=0){
  tmp <- data.frame(bed$INFO_INV5, bed$INFO_INV3, ID_caller = bed$ID_caller)
  strand1[bed$INFO_INV5] <- "-";strand2[bed$INFO_INV5] <- "-"
  strand1[bed$INFO_INV3] <- "+";strand2[bed$INFO_INV3] <- "+"
  tmp3 <- cbind(tmp, strand1, strand2)
}else{
  tmp3 <- cbind(bed, strand1, strand2)
}

bed2 <- bed[match(tmp3$ID_caller, bed$ID_caller),]
bed2$strand1 <- tmp3$strand1
bed2$strand2 <- tmp3$strand2

if((sum(bed2$INFO_SVTYPE %in% c("DEL","DUP")) != 0) & sum(is.na(bed2$strand1))!=0){
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
#' This function read bed format
#'
#' @param bed data frame
#' @return data frame
#' @export
bed_to_bedpe <- function(bed){
  bed <- strands_standardisation(bed)

  ALT <- bed$ALT
  tmp <- gsub("\\:.*",'', ALT[grepl('\\[',ALT) | grepl('\\]',ALT)])
  tmp1 <- gsub(".*\\[",'', tmp)
  tmp2 <- gsub(".*\\]",'', tmp1)
  chrom1 = as.character(bed$CHROM)
  chrom2 <- chrom1
  chrom2[grepl('\\[',ALT) | grepl('\\]',ALT)] <- tmp2

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

  bedpe <- data.frame(chrom1=bed$CHROM,
                      pos1=bed$POS,
                      chrom2=paste0(tolower(substring(chrom2,1,3)), substring(chrom2,4,nchar(chrom2))),
                      pos2=pos2,
                      strand1 = strand1,
                      strand2 = strand2,
                      ALT=ALT,
                      stringsAsFactors = FALSE)

  bedpe <- data.frame(bedpe, bed[,!(colnames(bed) %in% c("CHROM","POS","INFO_END","ALT","strand1","strand2"))])
  colnames(bedpe) <- c("chrom1", "pos1","chrom2","pos2","strand1","strand2","ALT",
                       colnames(bed)[!(colnames(bed) %in% c("CHROM","POS","INFO_END","ALT","strand1","strand2"))])
  return(bedpe)
}

#' Prepare breakpoint to use for bedtools
#'
#' This function read bed format
#'
#' @param standard_bedpe data frame
#' @param BND_diff the maximun value of breakpoint difference
#' @return data frame
#' @export
Standard_bedtool_prepare_bkpt <- function(standard_bedpe, BND_diff){
  diff <- BND_diff/2
  standard_bed_tmp <- data.frame(
    chrom = c(standard_bedpe$chrom1, standard_bedpe$chrom2),
    start = c(standard_bedpe$pos1-diff, standard_bedpe$pos2-diff),
    end = c(standard_bedpe$pos1+diff, standard_bedpe$pos2+diff),
    SVTYPE = c(standard_bedpe$SVTYPE, standard_bedpe$SVTYPE),
    ID = c(as.character(standard_bedpe$ID), as.character(standard_bedpe$ID_mate)),
    ID_mate = c(as.character(standard_bedpe$ID_mate), as.character(standard_bedpe$ID)), stringsAsFactors = FALSE)
  if (sum(standard_bed_tmp$start<0)!=0){standard_bed_tmp[standard_bed_tmp$start<0,]$start <- 0}
  return(standard_bed_tmp)
}

#' Check SV types and breakpoint positions
#'
#' This function read bed format
#'
#' @param intersect_file file name
#' @param SVTYPE_ignore whether ignore SV type for integration
#' @return data frame
#' @export
TypePosfilter <- function(intersect_file, SVTYPE_ignore){
  intersect <- read.table(intersect_file,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
  colnames(intersect) <- c("Caller1_CHROM", "Caller1_POS","Caller1_END","Caller1_SVTYPE","Caller1_ID","Caller1_ID_mate","Caller2",
                           "Caller2_CHROM", "Caller2_POS","Caller2_END","Caller2_SVTYPE","Caller2_ID","Caller2_ID_mate","overlap")
  if(SVTYPE_ignore){
    intersect_Typefilter <- intersect
  }else{
    intersect_Typefilter <- intersect[(intersect$Caller1_SVTYPE == intersect$Caller2_SVTYPE)|
                                        grepl("BND",intersect$Caller1_SVTYPE) | grepl("BND",intersect$Caller2_SVTYPE),]

  }
  if(nrow(intersect_Typefilter)!=0){
    ### Check Caller_1 ID_mate also listed in Caller_1 ID
    tmp_index <- lapply(intersect_Typefilter$Caller1_ID_mate,function(x) which(x==intersect_Typefilter$Caller1_ID))
    ### Further check Caller_2 ID_mate also listed in Caller_2 ID, which matched with Caller_1 ID
    BND_ID_match <- vector(length=nrow(intersect_Typefilter))
    for (i in 1: nrow(intersect_Typefilter)){
      BND_ID_match[i] <- intersect_Typefilter$Caller2_ID_mate[i] %in% intersect_Typefilter$Caller2_ID[tmp_index[[i]]]
    }
    intersect_TypePosBNDfilter <- intersect_Typefilter[BND_ID_match,]
  }else{
    intersect_TypePosBNDfilter <- intersect_Typefilter
  }
  return(intersect_TypePosBNDfilter)
}

#' Integrate SV call sets
#'
#' This function read bed format
#'
#' @param sampleID sample ID
#' @param SVCaller_name name of callers
#' @param SVCaller_bed_name names of call sets
#' @param BND_diff maximun breakpoint difference
#' @param bkpt_T_callers threshold of breakpoint difference
#' @param SVTYPE_ignore whether ignore SV type for integration
#' @param bedtools_dir directory of bedtools
#' @return data frame
#' @export
SVCaller_union_intersect_generate <- function(sampleID, SVCaller_name,SVCaller_bed_name,BND_diff,bkpt_T_callers,SVTYPE_ignore,bedtools_dir){
  ### Each bed, convert to bed_tmp and written to bed_tmp file
  for (i in 1:length(SVCaller_name)){
    #assign(paste0(SVCaller_name[i],"_standard_bedpe"), eval(parse(text=SVCaller_bed_name[i])))
    assign(paste0(SVCaller_name[i],"_bed_tmp"), Standard_bedtool_prepare_bkpt(eval(parse(text=paste0(SVCaller_name[i],"_standard_bedpe"))),BND_diff))
    write.table(eval(parse(text=paste0(SVCaller_name[i],"_bed_tmp"))), paste0(sampleID, "_", SVCaller_name[i],"_tmp.bed"), quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
  }

  ### Union set
  # generate bed and bed with new name, bed with new name bed_tmp and written to bed_tmp file
  SVCaller_name_all <- paste0(paste0(SVCaller_name,collapse = ""),"ALL")
  SVCaller_bed_all <- do.call("rbind", lapply(paste0(SVCaller_name,"_standard_bedpe"),function(s) eval(parse(text=s))))
  SVCaller_bed_all_newID <- SVCaller_bed_all
  #SVCaller_bed_all_newID <- SVCaller_bed_newID_generate2(SVCaller_bed_all,SVCaller_name_all)
  SVCaller_bed_all_newID_tmp <- Standard_bedtool_prepare_bkpt(SVCaller_bed_all_newID,BND_diff)
  write.table(SVCaller_bed_all_newID_tmp, paste0(sampleID, "_", SVCaller_name_all,"_tmp.bed"), quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)

  ### Intersect set
  # bedtools intersect union set bed_tmp with all SV caller bed_tmp
  intersect_file <- paste0(sampleID, "_", "all_",paste0(SVCaller_name,collapse = "_"),"_intersect.bed")
  overlap_f <- (BND_diff - bkpt_T_callers)/BND_diff
  system(paste(bedtools_dir,"intersect -a", paste0(sampleID, "_", SVCaller_name_all,"_tmp.bed"),
               "-b", paste(paste0(sampleID, "_", SVCaller_name,"_tmp.bed"), collapse = " "),
               "-names",paste(SVCaller_name,collapse = " "), "-f",overlap_f, "-wo >", intersect_file))

  intersect_filter <- TypePosfilter(intersect_file, SVTYPE_ignore)
  ### remove ID in tmp bed, only the ID in original bed
  intersect_filter <- intersect_filter[intersect_filter$Caller1_ID %in% SVCaller_bed_all_newID$ID,]

  ############################################
  ### work for more than two sv caller overlapping ###
  data <- c()
  for(i in 1: length(unique(intersect_filter$Caller1_ID))){
    tmp <- paste(intersect_filter[intersect_filter$Caller1_ID == unique(intersect_filter$Caller1_ID)[i],]$Caller2_ID)
    #data <- rbind(data,
    #              c(unique(intersect_filter$Caller1_ID)[i],
    #                unlist(lapply(SVCaller_name,function(s)paste(tmp[grepl(s,tmp)],collapse = ",")))))
    tmp2 <- unlist(lapply(SVCaller_name, function(s)paste(tmp[grepl(s,tmp)],collapse = ",")))
    tmp2[tmp2 == ""] <- NA
    data <- rbind(data,
                  c(unique(intersect_filter$Caller1_ID)[i], tmp2))
  }

  #### New simpler method WORK!!!!
  data <- data.frame(data)
  colnames(data) <- c("All_callers", SVCaller_name)
  tmp <- c()
  for(i in 1:length(SVCaller_name)){
    tmp <- rbind(tmp, eval(parse(text=paste0(SVCaller_name[i],"_standard_bedpe")))[as.character(eval(parse(text=paste0(SVCaller_name[i],"_standard_bedpe$ID"))))
                                                                                   %in% data$All_callers,])
  }
  union_tmp <- cbind(tmp, data[,-1])
  union <- union_tmp[!duplicated(data[,-1]),]
  for(i in (ncol(tmp)+1) :ncol(union)){
    union <- union[!(!is.na(union[,i]) & duplicated(union[,i])),]
  }

  SVCaller_bed_union <- union
  SVCaller_bed_combine_all <- SVCaller_bed_union
  return (SVCaller_bed_combine_all)
}

