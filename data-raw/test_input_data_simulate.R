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
