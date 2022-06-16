#' SVtype count
#'
#' This function read bed format
#'
#' @param bedpe bed format
#' @return data frame
#' @export
SVTYPE_stat_generate <- function(bedpe){
  All_SVTYPE <- unique(bedpe$SVTYPE)
  for(i in c(1: length(All_SVTYPE))){
    assign(paste0("N_", All_SVTYPE[i]), sum(bedpe$SVTYPE == All_SVTYPE[i]))
  }
  STAT_bed <- data.frame(do.call("cbind", lapply(paste0("N_", All_SVTYPE),function(s) eval(parse(text=s)))))
  colnames(STAT_bed) <- All_SVTYPE
  return(STAT_bed)
}

#' Summary SV count
#'
#' This function read bed format
#'
#' @param All_sampleID sample ID
#' @param vcf_list vcf file names
#' @return data frame
#' @export
Summary_SV_type <- function(All_sampleID, vcf_list){
  summary_results <- c()
  for(i in c(1:length(All_sampleID))){
    sampleID <- All_sampleID[i]
    vcf_file <- vcf_list[i]
    SVTYPE_count <- simple_SVTYPE_classification(vcf_file, caller_name = "StructuralVariantUtil")[[2]]

    all_colnames <- unique(c(colnames(summary_results), colnames(SVTYPE_count)))

    SVTYPE_count_tmp <- data.frame(matrix(0, nrow = nrow(SVTYPE_count), ncol = length(all_colnames)))
    colnames(SVTYPE_count_tmp) <- all_colnames

    index <- match(colnames(SVTYPE_count), colnames(SVTYPE_count_tmp))
    index_SVTYPE <- index[!(is.na(index))]
    SVTYPE_count_tmp[,index_SVTYPE] <- SVTYPE_count

   # SVTYPE_count_tmp[colnames(SVTYPE_count_tmp) %in% colnames(SVTYPE_count)] <- SVTYPE_count

    if(i > 1){
      summary_results_tmp <- data.frame(matrix(0, nrow = nrow(summary_results), ncol = length(all_colnames)))
      colnames(summary_results_tmp) <- all_colnames

      index <- match(colnames(summary_results), colnames(summary_results_tmp))
      index_SVTYPE <- index[!(is.na(index))]
      summary_results_tmp[,index_SVTYPE] <- summary_results


      #summary_results_tmp[colnames(summary_results_tmp) %in% colnames(summary_results)] <- summary_results
      summary_results <- summary_results_tmp
    }
    summary_results <- rbind(summary_results, SVTYPE_count_tmp)
  }
  summary_results <- data.frame(sampleID =  All_sampleID, summary_results)
  return(summary_results)
}

#' Spectrum of SV count
#'
#' This function read bed format
#'
#' @param All_sampleID sample ID
#' @param vcf_list vcf file names
#' @param threshold_total threshold of minimum total count of SVs per sample
#' @param threshold_relative_freq threshold of minimum relative frequency of one SV type
#' @return data frame of hyper SV
#' @export
Spectrum_SV_type <- function(All_sampleID, vcf_list, threshold_total, threshold_relative_freq){
  input_SV_count <- Summary_SV_type(All_sampleID, vcf_list)
  theme1 <-  ggplot2::theme(axis.text=ggplot2::element_text(size=12,face="bold"),
                            axis.title=ggplot2::element_text(size=14,face="bold"),
                            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 12),
                            plot.title = ggplot2::element_text(size=14),
                            legend.text = ggplot2::element_text(size=12,face="bold"),
                            #legend.title = element_text(size=12,face="bold"),
                            legend.title = ggplot2::element_blank(),
                            legend.position="top")

  theme2 <-  ggplot2::theme(axis.text = ggplot2::element_text(size=12,face="bold"),
                            axis.title=ggplot2::element_text(size=14,face="bold"),
                            plot.title = ggplot2::element_text(size=14),
                            legend.text = ggplot2::element_text(size=12,face="bold"),
                            #legend.title = element_text(size=12,face="bold"),
                            legend.title = ggplot2::element_blank(),
                            legend.position="top")
  df <- data.frame(SVTYPE = rep(colnames(input_SV_count)[2:ncol(input_SV_count)], each = nrow(input_SV_count)),
                   sampleID = rep(input_SV_count$sampleID, (ncol(input_SV_count)-1)),
                   Count = as.vector(unlist(input_SV_count[,2:ncol(input_SV_count)])))

  p1 <- ggplot2::ggplot(data=df, ggplot2::aes(x=sampleID, y=Count, fill=SVTYPE)) +
    ggplot2::geom_bar(stat="identity")+
    ggplot2::scale_fill_manual(values = RColorBrewer::brewer.pal(12, "Set3"), drop=F)+
    theme1
  pdf(file="./Spectrum_SV_across_sample.pdf", width=19, height=9)
  print(p1)
  dev.off()

  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = SVTYPE, y=Count)) +
    ggplot2::geom_boxplot()+
    theme2

  N_total <- rowSums(input_SV_count[,2:ncol(input_SV_count)])
  if(missing(threshold_total)){
    threshold_total <- mean(N_total)
  }
  if(missing(threshold_relative_freq)){
    threshold_relative_freq <- 0.5
  }

  df2 <- cbind(df,
               relative_freq = df$Count/(rep(N_total,(ncol(input_SV_count)-1))),
               N_total = rep(N_total,(ncol(input_SV_count)-1)))
  p2 <- ggplot2::ggplot(df2, ggplot2::aes(x = SVTYPE, y=relative_freq)) +
    ggplot2::geom_boxplot()+
    ggplot2::ylab("Relative frequency")+
    ggplot2::scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1))+
    theme2

  figure1 <- ggpubr::ggarrange(p1, p2,
                               labels = c("A", "B"), heights = c(1,1),
                               ncol = 2, nrow = 1)

  pdf(file="./Spectrum_SV_across_SVTYPE.pdf",width=9,height=6)
  print(figure1)
  dev.off()

  hyper_SV <- df2[df2$relative_freq > threshold_relative_freq & df2$N_total > threshold_total,]
  hyper_SV$HYPER_SVTYPE <- paste0("hyper-", hyper_SV$SVTYPE)
  return(list(input_SV_count, hyper_SV))
}

#' Integrate SV and CNV
#'
#' This function integrate deletion and duplication with copy number segments
#'
#' @param sampleID name of sample
#' @param SV_bed data frame of SV
#' @param CNV_bed data frame of CNV segment
#' @param bedtools_dir bedtools for use
#' @param overlap_f the fraction of minimum overlap required of CNV segment as a fraction of SV
#' @return data frame of SV set
#' @export
SV_CNV_integration <- function(sampleID, SV_bed, CNV_bed, bedtools_dir, overlap_f){
  directory <- "./"
  sub_directory <- paste0("./tmp/")
  dir.create(sub_directory)

  assign(paste0(sampleID, "_CNV_tmp.bed"), data.frame(chrom = CNV_bed$chrom,
                                                      start = CNV_bed$start,
                                                      end = CNV_bed$end,
                                                      cn = CNV_bed$cn))
  write.table(eval(parse(text=paste0(sampleID, "_CNV_tmp.bed"))), paste0(sub_directory,sampleID,"_CNV_tmp.bed"), quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
  assign(paste0(sampleID,"_SV_tmp.bed"), data.frame(chrom = SV_bed[SV_bed$SVTYPE %in% c("DEL","DUP"),]$chrom1,
                                                    start = SV_bed[SV_bed$SVTYPE %in% c("DEL","DUP"),]$pos1,
                                                    end = SV_bed[SV_bed$SVTYPE %in% c("DEL","DUP"),]$pos2,
                                                    ID = SV_bed[SV_bed$SVTYPE %in% c("DEL","DUP"),]$ID,
                                                    SVTYPE = SV_bed[SV_bed$SVTYPE %in% c("DEL","DUP"),]$SVTYPE))

  write.table(eval(parse(text=paste0(sampleID,"_SV_tmp.bed"))), paste0(sub_directory,sampleID,"_SV_tmp.bed"), quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)

  intersect_file <- paste0(sampleID, "_", "SV_CNV","_intersect.bed")
  system(paste(bedtools_dir,"intersect -a", paste0(sub_directory,sampleID,"_SV_tmp.bed"),
               "-b", paste0(sub_directory,sampleID,"_CNV_tmp.bed"),
               "-f", overlap_f,"-r", "-wo >", intersect_file))

  if(file.info(intersect_file)$size != 0){
    intersect <- read.table(intersect_file,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
    colnames(intersect) <- c("SV_chrom", "SV_start","SV_end","SV_ID","SV_type",
                             "CNV_chrom", "CNV_start","CNV_end","CNV_cn","overlap")
    df_SV_CNV <- cbind(sampleID = sampleID, intersect)
  }
  return(df_SV_CNV)
}


#' SV breakpoints gene annotation
#'
#' This function annotate SV breakpoints based on gene regions
#'
#' @param input_df_name name of input data frame
#' @param gene_bed data frame of gene regions
#' @param bedtools_dir bedtools for use
#' @return data frame of SV set
#' @export
SV_breakpoint_gene_annotation <- function(input_df_name, gene_bed, bedtools_dir){
  directory <- "./"
  sub_directory <- paste0("./tmp/")
  dir.create(sub_directory)

  gene_file <- paste0(sub_directory,"gene","_tmp.bed")
  write.table(gene_bed, gene_file, quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)
  bedpe <- eval(parse(text = input_df_name))
  if(nrow(bedpe) !=0){
    bedpe$chrom1 <- as.character(bedpe$chrom1)
    bedpe$chrom2 <- as.character(bedpe$chrom2)
    bedpe$ID <- as.character(bedpe$ID)
    bedpe$ID_mate <- as.character(bedpe$ID_mate)

    SV_bed <- data.frame(chrom = c(bedpe$chrom1, bedpe$chrom2),
                         start = c(bedpe$pos1-1, bedpe$pos2-1),
                         end = c(bedpe$pos1, bedpe$pos2),
                         SVTYPE = c(bedpe$SVTYPE, bedpe$SVTYPE),
                         ID = c(bedpe$ID, bedpe$ID_mate))

    write.table(SV_bed, paste0(sub_directory,"SV_tmp.bed"), quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)

    intersect_file <- paste0(directory,"SV_gene","_intersect.bed")
    system(paste(bedtools_dir,"intersect -a", paste0(sub_directory,"SV_tmp.bed"),
                 "-b", gene_file,
                 "-wao >", intersect_file))
    intersect <- read.table(intersect_file,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
    colnames(intersect) <- c("SV_chrom", "SV_start","SV_end","SVTYPE", "SV_ID",
                             colnames(gene_bed),"overlap")
  }else{
    intersect <- c()
  }
  return(intersect)
}

#' SV breakpoints gene annotation on input bed format
#'
#' This function annotate SV breakpoints based on gene regions
#'
#' @param input_df_name name of input data frame
#' @param gene_bed data frame of gene regions
#' @param bedtools_dir bedtools for use
#' @return data frame of SV set
#' @export
SV_bedpe_gene_annotation <- function(input_df_name, gene_bed, bedtools_dir){
  bedpe <- eval(parse(text = input_df_name))
  bedpe_bkpt_geneAnnotated <- SV_breakpoint_gene_annotation(input_df_name, gene_bed,  bedtools_dir)

  bedpe_bkpt_geneAnnotated <- bedpe_bkpt_geneAnnotated[bedpe_bkpt_geneAnnotated$type == "gene" & bedpe_bkpt_geneAnnotated$gene_biotype == "protein_coding",]
  if(nrow(bedpe) !=0){
    bedpe_geneAnnotated <- c()
    for(j in c(1:nrow(bedpe))){
      tmp1 <- bedpe_bkpt_geneAnnotated[which(match(bedpe_bkpt_geneAnnotated$SV_ID, bedpe$ID) == j),]
      #if(nrow(tmp1) == 0){pos1_overlap_gene <- NA}else{pos1_overlap_gene = paste(tmp1$gene_name, collapse= ",")}

      if(nrow(tmp1) == 0){pos1_overlap_gene <- NA}else{pos1_overlap_gene = tmp1$gene_name}
      tmp3 <- cbind(bedpe[j,], pos1_overlap_gene)

      tmp2 <- bedpe_bkpt_geneAnnotated[which(match(bedpe_bkpt_geneAnnotated$SV_ID, bedpe$ID_mate) == j),]
      #if(nrow(tmp2) == 0){pos2_overlap_gene <- NA}else{pos2_overlap_gene = paste(tmp2$gene_name, collapse= ",")}

      if(nrow(tmp2) == 0){pos2_overlap_gene <- NA}else{pos2_overlap_gene = tmp2$gene_name}
      if(length(pos2_overlap_gene)>1){
        tmp4 <- cbind(sapply(tmp3,rep.int,times=length(pos2_overlap_gene)), pos2_overlap_gene = rep(pos2_overlap_gene, each = nrow(tmp3)))
      }else{
        tmp4 <- cbind(tmp3, pos2_overlap_gene)
      }
      bedpe_geneAnnotated <- rbind(bedpe_geneAnnotated, tmp4)
    }
  }else{
    bedpe_geneAnnotated <- c()
  }
  return(bedpe_geneAnnotated)
}

#' SV breakpoints summary for all samples
#'
#' This function put all sample SV breakpoints together
#'
#' @param All_sampleID sample ID for all samples
#' @param All_input_df_name names of all bed in df
#' @return data frame of SV set
#' @export
Summary_SV_breakpoint <- function(All_sampleID, All_input_df_name){
  df_breakpoints <- c()
  for(i in c(1 : length(All_sampleID))){
    sampleID <- All_sampleID[i]
    bedpe <- eval(parse(text = All_input_df_name[i]))
    if(nrow(bedpe) != 0){
      df_breakpoints <- rbind(df_breakpoints,
                              data.frame(sampleID = sampleID,
                                        chrom = c(as.character(bedpe$chrom1), as.character(bedpe$chrom2)),
                                         pos = c(bedpe$pos1, bedpe$pos2),
                                          ID = c(bedpe$ID, bedpe$ID_mate)))
    }

  }
  return(df_breakpoints)
}

#' Genomic bins
#'
#' Put SV breakpoint into bins
#'
#' @param df_breakpoints data frame of SV breakpoints
#' @return data frame of SV set
#' @export
Spectrum_SV_bin <- function(df_breakpoints){
  df_bin_all <- c()
  for(chrom in c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                 "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")){
    df2 <- df_breakpoints[df_breakpoints$chrom %in% chrom,]
    # breaks = seq(0,250e6,1e5)
    breaks = seq(0, 250e6, 1e6)
    bin_labels <- cut(df2$pos, breaks=breaks, labels=c(1:(length(breaks)-1)))
    bin <- cut(df2$pos, breaks=breaks)
    df_bin <- cbind(bin_labels,bin, breaks = breaks[bin_labels],df2)
    #df_bin <- cbind(df2, bin_labels, bin, breaks = breaks[bin_labels])
    bin_table1 <- data.frame(table(df_bin$bin_labels))

    df_sample <- df_bin[!(duplicated(paste0(df_bin$sampleID,"_", df_bin$bin_labels))),]
    bin_table2 <- data.frame(table(df_sample$bin_labels))

    df_bin <- cbind(df_bin,
                    count_breakpoints = bin_table1[df_bin$bin_labels,]$Freq,
                    count_sample = bin_table2[df_bin$bin_labels,]$Freq)
    df_bin_all <- rbind(df_bin_all,df_bin)
  }

  df_bin_all <- data.frame(chrom_bin_labels = paste0(df_bin_all$chrom,"_", df_bin_all$bin_labels), df_bin_all)
  return(df_bin_all)
}

#' Genomic bins hotspots
#'
#' Define hotspots
#'
#' @param df_breakpoints data frame of SV breakpoints
#' @param  threshold_count_breakpoint threshold of number of SD
#' @param  threshold_count_sample threshold of number of samples
#' @return data frame of SV set
#' @export
Spectrum_SV_bin_define_hotspot <- function(df_breakpoints, threshold_count_breakpoint, threshold_count_sample){
  df_bin_all <- Spectrum_SV_bin(df_breakpoints)
  #head(df_bin_all)
  ##### add HOTSPOT information
  df2 <- df_bin_all
  df3 <- df2
  df3 <- df3[!duplicated(df3$chrom_bin_labels),]
  SV_hotspots <- rbind(df3[df3$count_breakpoints>mean(df3$count_breakpoints)+threshold_count_breakpoint*sd(df3$count_breakpoints),],
                       df3[df3$count_sample>mean(df3$count_sample)+threshold_count_sample*sd(df3$count_sample),])

  hotspots <- unique(SV_hotspots$chrom_bin_labels)
  df_bin_all_hotspot <- data.frame(df_bin_all,
                                   is_hotspot_breakpoint = df_bin_all$chrom_bin_labels %in% df3[df3$count_breakpoints>mean(df3$count_breakpoints)+threshold_count_breakpoint*sd(df3$count_breakpoints),]$chrom_bin_labels,
                                   is_hotspot_sample = df_bin_all$chrom_bin_labels %in% df3[df3$count_sample>mean(df3$count_sample)+threshold_count_sample*sd(df3$count_sample),]$chrom_bin_labels,
                                   is_hotspot = df_bin_all$chrom_bin_labels %in% hotspots)
  return(df_bin_all_hotspot)
}

#' Plot genomic bins hotspots
#'
#' Plot hotspots
#'
#' @param df_bin_all_hotspots data frame of SV breakpoints bins with hotspots defined
#' @return data frame of SV set
#' @export
plot_ideograms <- function(df_bin_all_hotspots){
  df2 <- df_bin_all_hotspots
  df2$chrom <- factor(df2$chrom, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                                            "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"))
  is_hotspot <- c(FALSE, TRUE)

  col_all <- RColorBrewer::brewer.pal(8, "Dark2")[c(8,1)]
  cex_all <- c(1,2)

  chrom_all <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                 "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")

  png(file="./SV_breakpoints_ideogram_byBins.png",width=2000,height=550)
  kp <- karyoploteR::plotKaryotype("hg38", chromosomes = chrom_all, plot.type = 4, cex = 1.5, srt = 30)
  karyoploteR::kpAddBaseNumbers(kp, tick.dist = 50000000)
  karyoploteR::kpDataBackground(kp, data.panel = 1)

  for(i in c(1: length(chrom_all))){
    for(j in c(1,2)){
      df_tmp <- df2[df2$chrom==chrom_all[i] & df2$is_hotspot_breakpoint == is_hotspot[j],]
      x <- df_tmp$breaks
      y <- df_tmp$count_breakpoints/120
      karyoploteR::kpPoints(kp, chr=chrom_all[i], x=x, y=y, data.panel = 1, col = col_all[j],  cex=cex_all[j])
    }
  }
  karyoploteR::kpAxis(kp, ymin=0, ymax=120,numticks = 3, data.panel = 1, cex = 1.5)
  dev.off()

  ######## ideogram with sample counts for each chromosome by bins
  df3 <- df2
  df3 <- df3[!duplicated(df3$chrom_bin_labels),]
  chrom_all <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                 "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")

  png(file="./SV_samples_ideogram_byBins.png",width=2000,height=500)
  kp <- karyoploteR::plotKaryotype("hg38", chromosomes = chrom_all, plot.type = 4, cex = 1.5, srt = 30)
  #kpAddBaseNumbers(kp)
  karyoploteR::kpAddBaseNumbers(kp, tick.dist = 50000000)
  karyoploteR::kpDataBackground(kp, data.panel = 2)

  for(i in c(1: length(chrom_all))){
    for(j in c(1,2)){
      df_tmp <- df3[df3$chrom==chrom_all[i] & df3$is_hotspot_sample == is_hotspot[j],]
      #x <- df_tmp$pos
      x <- df_tmp$breaks
      y <- df_tmp$count_sample/50
      karyoploteR::kpPoints(kp, chr=chrom_all[i], x=x, y=y, data.panel = 2, col = col_all[j], cex=cex_all[j])
    }
  }
  karyoploteR::kpAxis(kp, ymin=0, ymax=50,numticks = 3, data.panel = 2, cex = 1.5)
  dev.off()
}

#' Genomic bins hotspots
#'
#' Define hotspots
#'
#' @param All_sampleID sample ID for all samples
#' @param All_input_df_name names of all bed in df
#' @param  threshold_count_breakpoint threshold of number of SD
#' @param  threshold_count_sample threshold of number of samples
#' @return data frame of genomic bins with hotspots defined
#' @export
Spectrum_SV_bin_generate <- function(All_sampleID, All_input_df_name, threshold_count_breakpoint,threshold_count_sample){
  df_breakpoints <- Summary_SV_breakpoint(All_sampleID, All_input_df_name) ### put all bed df together
  df_bin_all <- Spectrum_SV_bin(df_breakpoints) ### generate genomic bins
  df_bin_all_hotspots <- Spectrum_SV_bin_define_hotspot(df_breakpoints, threshold_count_breakpoint,threshold_count_sample) ##### define genomic bins and add HOTSPOT information
  write.table(df_bin_all_hotspots,"./df_bin_all_hotspots.txt", quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)
  plot_ideograms(df_bin_all_hotspots) ### SV hotspot visualisation
  return(df_bin_all_hotspots)
}

#' Define frequency of gene fusions
#'
#' Define frequency of gene fusions
#'
#' @param bed_GeneAnnotated data frame
#' @return data frame of summarised gene fusions
#' @export
Summary_gene_fusions <- function(bed_GeneAnnotated){
  if(!("sampleID" %in% colnames(bed_GeneAnnotated))){
    sampleID <- "sample_1"
    df_All_gene_fusions <- cbind(sampleID, bed_GeneAnnotated)
  }else{
    df_All_gene_fusions <- bed_GeneAnnotated
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
