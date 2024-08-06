# load packages
library(rehh)
library(tidyverse)

# read in data for each 

ST <- data2haplohh(hap_file = "~/Downloads/ShapeitOutputRinput/Daff_ST.vcf.gz",
                         polarize_vcf = FALSE)

SR1 <- data2haplohh(hap_file = "~/Downloads/ShapeitOutputRinput/Daff_SR1.vcf.gz",
                   polarize_vcf = FALSE)

SR2 <- data2haplohh(hap_file = "~/Downloads/ShapeitOutputRinput/Daff_SR2.vcf.gz",
                   polarize_vcf = FALSE)


# filter on MAF - here 0.05

ST_f <- subset(ST, min_maf = 0.05)
SR1_f <- subset(SR1, min_maf = 0.05)
SR2_f <- subset(SR2, min_maf = 0.05)

# perform scans

ST_scan <- scan_hh(ST_f, polarized = FALSE)
SR1_scan <- scan_hh(SR1_f, polarized = FALSE)
SR2_scan <- scan_hh(SR2_f, polarized = FALSE)

# perform iHS 

ST_ihs <- ihh2ihs(ST_scan, freqbin = 1)
SR1_ihs <- ihh2ihs(SR1_scan, freqbin = 1)
SR2_ihs <- ihh2ihs(SR2_scan, freqbin = 1)

# plot

library(ggplot2)

ggplot(ST_ihs$ihs, aes(POSITION, LOGPVALUE)) + geom_point()
ggplot(ST_ihs$ihs, aes(POSITION, IHS)) + geom_point()

ggplot(SR1_ihs$ihs, aes(POSITION, LOGPVALUE)) + geom_point()
ggplot(SR1_ihs$ihs, aes(POSITION, IHS)) + geom_point()

ggplot(SR2_ihs$ihs, aes(POSITION, LOGPVALUE)) + geom_point()
ggplot(SR2_ihs$ihs, aes(POSITION, IHS)) + geom_point()

# perform xp-ehh

ST_SR1 <- ies2xpehh(ST_scan, SR1_scan,
                    popname1 = "ST", popname2 = "SR1",
                    include_freq = T)
ST_SR2 <- ies2xpehh(ST_scan, SR2_scan,
                    popname1 = "ST", popname2 = "SR2",
                    include_freq = T)
SR1_SR2 <- ies2xpehh(SR1_scan, SR2_scan,
                    popname1 = "SR1", popname2 = "SR2",
                    include_freq = T)

# plot
ggplot(ST_SR1, aes(POSITION, XPEHH_ST_SR1)) + geom_point()
ggplot(ST_SR2, aes(POSITION, XPEHH_ST_SR2)) + geom_point()
ggplot(SR1_SR2, aes(POSITION, XPEHH_SR1_SR2)) + geom_point()


# find the highest hit

hit_ST_SR1 <- ST_SR1 %>% 
  arrange(desc(LOGPVALUE)) %>% top_n(1)
hit_ST_SR2 <- ST_SR2 %>% 
  arrange(desc(LOGPVALUE)) %>% top_n(1)
hit_SR1_SR2 <- SR1_SR2 %>% 
  arrange(desc(LOGPVALUE)) %>% top_n(1)

# get SNP position

x_ST_SR1 <- hit_ST_SR1$POSITION
x_ST_SR2 <- hit_ST_SR2$POSITION
x_SR1_SR2 <- hit_SR1_SR2$POSITION


marker_id_ST_SR1_1 <- which(ST_f@positions == x_ST_SR1)
marker_id_ST_SR1_2 <- which(SR1_f@positions == x_ST_SR1)

marker_id_ST_SR2_1 <- which(ST_f@positions == x_ST_SR2)
marker_id_ST_SR2_2 <- which(SR2_f@positions == x_ST_SR2)

marker_id_SR1_SR2_1 <- which(SR1_f@positions == x_SR1_SR2)
marker_id_SR1_SR2_2 <- which(SR2_f@positions == x_SR1_SR2)



ST_furcation_ST_SR1 <- calc_furcation(ST_f, 
                               mrk = marker_id_ST_SR1_1)
ST_haplen_ST_SR1 <- calc_haplen(ST_furcation_ST_SR1)
plot(ST_haplen_ST_SR1)

SR1_furcation_ST_SR1 <- calc_furcation(SR1_f, 
                                      mrk = marker_id_ST_SR1_2)
SR1_haplen_ST_SR1 <- calc_haplen(SR1_furcation_ST_SR1)
plot(SR1_haplen_ST_SR1)





SR1_furcation_SR1_SR2 <- calc_furcation(SR1_f, 
                                      mrk = marker_id_SR1_SR2_1)
SR1_haplen_SR1_SR2 <- calc_haplen(SR1_furcation_SR1_SR2)
plot(SR1_haplen_SR1_SR2)

SR2_furcation_SR1_SR2 <- calc_furcation(SR2_f, 
                                       mrk = marker_id_SR1_SR2_2)
SR2_haplen_SR1_SR2 <- calc_haplen(SR2_furcation_SR1_SR2)
plot(SR2_haplen_SR1_SR2)




# write out xpEHH

ST_SR1 <- tibble::as_tibble(ST_SR1)
colnames(ST_SR1) <- tolower(colnames(ST_SR1))

write_tsv(ST_SR1, "~/Downloads/ShapeitOutputRinput/ST_SR1_xpEHH.tsv")



ST_SR2 <- tibble::as_tibble(ST_SR2)
colnames(ST_SR2) <- tolower(colnames(ST_SR2))

write_tsv(ST_SR2, "~/Downloads/ShapeitOutputRinput/ST_SR2_xpEHH.tsv")



SR1_SR2 <- tibble::as_tibble(SR1_SR2)
colnames(SR1_SR2) <- tolower(colnames(SR1_SR2))

write_tsv(SR1_SR2, "~/Downloads/ShapeitOutputRinput/SR1_SR2_xpEHH.tsv")
