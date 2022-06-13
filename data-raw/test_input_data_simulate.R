### Test data for ShinySoSV prediction
#Data was generated randomly with normal distribution:
#VAF with mean of 0.5 and SD of 0.1;
#normal coverage with mean of 30 and SD of 10 and
#tumour coverage with mean of 60 and SD of 10.
#The BND threshold was set as 100 bp.
set.seed(1)
ShinySoSV_newdata <- data.frame(sampleID = paste0("sample_",c(1:100)),
                                VAF = round(rnorm(100, mean=0.5, sd=0.1),digits = 2),
                                N_coverage = round(rnorm(100, mean=30, sd=10),digits = 2),
                                T_coverage = round(rnorm(100, mean=60, sd=10),digits = 2),
                                BND_threshold = 100)
usethis::use_data(ShinySoSV_newdata,overwrite = TRUE)

### Test data for simple Sv type classification
vcf_file <- system.file("extdata",
                        "GRIDSS_SVEngine_TumorSV2.60x_NormalSV1.60x_0.5_somatic_PASS_annotated.vcf",
                        package = "ShinySoSV2")
bed <- vcf_to_bed(vcf_file)
bed_test <- data.frame(CHROM = bed$CHROM,
                       POS = bed$POS,
                       ALT = bed$ALT)
usethis::use_data(bed_test,overwrite = TRUE)


### Test data for SV integration
vcf_file <- system.file("extdata",
                        "lumpy_SVEngine_TumorSV2.60x_NormalSV1.60x_0.5_TumorminSU4.vcf",
                        package = "ShinySoSV2")
CallerA_bed <- simple_SVTYPE_classification(bed = vcf_to_bed(vcf_file), "CallerA")
usethis::use_data(CallerA_bed,overwrite = TRUE)

vcf_file <- system.file("extdata",
                        "manta_SVEngine_TumorSV2.60x_NormalSV1.60x_0.5.T.PASS.recode.vcf",
                        package = "ShinySoSV2")
CallerB_bed <- simple_SVTYPE_classification(bed = vcf_to_bed(vcf_file), "CallerB")
usethis::use_data(CallerB_bed,overwrite = TRUE)

vcf_file <- system.file("extdata",
                        "GRIDSS_SVEngine_TumorSV2.60x_NormalSV1.60x_0.5_somatic_PASS_annotated.vcf",
                        package = "ShinySoSV2")
CallerC_bed <- simple_SVTYPE_classification(bed = vcf_to_bed(vcf_file), "CallerC")
usethis::use_data(CallerC_bed,overwrite = TRUE)

###Test data for Sv type composition
generateRandomPos <- function(n,chr,chr.sizes,width,strand){
  random_chr <- sample(x=chr,size=n,prob=chr.sizes,replace=T)
  random_pos <- sapply(random_chr,function(chrTmp){sample(chr.sizes[chr==chrTmp],1)})
  res <- GenomicRanges::GRanges(random_chr,IRanges::IRanges(random_pos, random_pos+width), strand = strand)
  return(res)
}

All_sampleID <- paste0("sample_",c(1:100))
set.seed(1)
input_SV_count <- data.frame(sampleID = All_sampleID,
                             DEL = sample.int(300, 100, replace = TRUE),
                             DUP = sample.int(300, 100, replace = TRUE),
                             INS = sample.int(100, 100, replace = TRUE),
                             INV = sample.int(300, 100, replace = TRUE),
                             TRA = sample.int(300, 100, replace = TRUE))
for(i in c(1: nrow(input_SV_count))){
  sample_tmp_bed <- c()
  for(SVTYPE in colnames(input_SV_count)[2:ncol(input_SV_count)]){
    if(SVTYPE == "TRA"){
      n <- 2*input_SV_count[i, colnames(input_SV_count) == SVTYPE]
    }else{
      n <- input_SV_count[i, colnames(input_SV_count) == SVTYPE]
    }

    SV_length <-sample(c(50,100,500,1000,2000,5000,10000,15000,100000,1000000), n, replace = TRUE, prob = c(0.1,0.1,0.2,0.2,0.3,0.5,0.5,0.4,0.2,0.1))
    res <- generateRandomPos(n = n, chr = c(seq(1:22),"X","Y"), chr.sizes = seq(0.2,24), width = SV_length, strand = "+")
    sample_tmp_bed <- rbind(sample_tmp_bed, data.frame(chrom = paste0("chr",res@seqnames), res@ranges, strand = res@strand, SVTYPE = SVTYPE))
  }
  write.table(sample_tmp_bed, paste0(input_SV_count$sampleID[i],"_tmp.bed"), quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
  system(paste0("/opt/homebrew/bin/bedtools shuffle ",
                "-excl hg38_gaps_centromeres_Telomeres.bed -i ",
                input_SV_count$sampleID[i],"_tmp.bed ","-g hg38.genome ",
                "> ",
                "random_",input_SV_count$sampleID[i],".bed"))

  random_sample.bed <- read.table(paste0("random_",input_SV_count$sampleID[i],".bed"))
  colnames(random_sample.bed) <- c("chrom1","pos1","pos2","SVLEN","strand1","SVTYPE")

  n_TRA <- sum(random_sample.bed$SVTYPE == "TRA")/2
  n_total_tmp <- nrow(random_sample.bed)

  df <- data.frame(chrom1 = random_sample.bed$chrom1[1:(n_total_tmp-n_TRA)],
                   pos1 = random_sample.bed$pos1[1:(n_total_tmp-n_TRA)],
                   chrom2 = c(random_sample.bed$chrom1[1:(n_total_tmp-2*n_TRA)], random_sample.bed$chrom1[(n_total_tmp-n_TRA+1):(n_total_tmp)]),
                   pos2 = c(random_sample.bed$pos2[1:(n_total_tmp-2*n_TRA)], random_sample.bed$pos2[(n_total_tmp-n_TRA+1):(n_total_tmp)]),
                   ID = paste0(c(1:(n_total_tmp-n_TRA)),"_1"),
                   ID_mate = paste0(c(1:(n_total_tmp-n_TRA)),"_2"),
                   SVTYPE = random_sample.bed$SVTYPE[1:(n_total_tmp-n_TRA)],
                   SVLEN = random_sample.bed$SVLEN[1:(n_total_tmp-n_TRA)],
                   strand = random_sample.bed$strand1[1:(n_total_tmp-n_TRA)])
  assign(paste0(All_sampleID[i], "_df"), df[!(df$SVTYPE == "TRA" & (df$chrom1 == df$chrom2)),])
}
save(list = paste0(All_sampleID, "_df"), file = "./input_SV_bed.RData")


set.seed(1)
input_SV_count <- data.frame(sampleID = paste0("sample_",c(1:100)),
                             DEL = sample.int(300, 100, replace = TRUE),
                             DUP = sample.int(300, 100, replace = TRUE),
                             INS = sample.int(100, 100, replace = TRUE),
                             INV = sample.int(300, 100, replace = TRUE),
                             TRA = sample.int(300, 100, replace = TRUE))
usethis::use_data(input_SV_count, overwrite = TRUE)

###Test data for CNV integration, currently use sample UP2003 in HRPCa project
SCNV <- read.table(system.file("extdata",
                               "UP2003-T.final.call.threshold.cns",
                               package = "ShinySoSV2"), header = TRUE)
CNV_bed <- SCNV[,c(1,2,3,6)]
usethis::use_data(CNV_bed, overwrite = TRUE)

bedpe <- read.table(system.file("extdata",
                                "UP2003_Manta_GRIDSS_intersect_both_high_confidence.bedpe",
                                package = "ShinySoSV2"), header = TRUE)
SV_bed <- bedpe[,c(1:10)]
usethis::use_data(SV_bed, overwrite = TRUE)


#data(CallerA_bed)
#load(system.file("extdata", "gene_bed.Rdata", package = "ShinySoSV2"))
#bedtools_dir = "/opt/homebrew/bin/bedtools"
#CallerA_bed_GeneAnnotated <- SV_bedpe_gene_annotation(input_df_name = "CallerA_bed", gene_bed, bedtools_dir)
#usethis::use_data(CallerA_bed_GeneAnnotated, overwrite = TRUE)

df_All_gene_fusions <- read.csv("/Users/tingtinggong/Desktop/work_at_home/HYPER-DUP/SomaticSV_MantaGRIDSS/All_gene_fusions.csv")
df_All_gene_fusions <- df_All_gene_fusions[!(df_All_gene_fusions$sampleID %in% c("16599-1159140", "10651-1042378","UP2063","UP2037","UP2041")),]
df_All_gene_fusions$sampleID<- gsub('-.*','',as.character(df_All_gene_fusions$sampleID))
df_All_gene_fusions$SVTYPE[df_All_gene_fusions$SVTYPE=="TRA_INV"] <- "TRA"
ProstateCancer_SV_bed_GeneAnnotated <- df_All_gene_fusions[,c(1:11,18,19)]
usethis::use_data(ProstateCancer_SV_bed_GeneAnnotated, overwrite = TRUE)


