## Plot BayPass results (SNPs and inversions)
## Kaichi Huang 2019 Aug

library(tidyverse)
library(gridExtra)
library(grid)
select=dplyr::select

# Read in datasets
my.data <- read.table("all_soilGEA_summary_betai_reg.out", h=T) # SNP GEA
my.pod <- read.table("pod_soilGEA_summary_betai_reg.out", h=T) # POD GEA
my.snp <- read.table("GSD_RAD.HA412.good_snps2.baypass.snp", h=F) # SNP table
names(my.snp) <- c("CHR", "POS")
my.snp$POS <- my.snp$POS/(10^6)
my.cov <- read.table("soil.baypass.cov", h=F, sep="\t") # Covariables
my.inversion <- read.table("inv_soilGEA_summary_betai_reg.out", h=T) # Inversion GEA
my.region <- read.table("../inv.bp.txt", h=T) # Inversion position
my.region$start <- my.region$start/(10^6)
my.region$end <- my.region$end/(10^6)

bp <- my.snp$POS[my.snp$CHR == 1]
mid <- 0+max(bp)/2
for (i in 2:length(levels(my.snp$CHR))) {
  max.now <- max(bp)
  mid <- c(mid, max.now+max(my.snp$POS[my.snp$CHR == i])/2)
  bp <- c(bp, my.snp$POS[my.snp$CHR == i] + max.now)
  my.region$bp1[my.region$CHROM == i] = my.region$start[my.region$CHROM == i] + max.now
  my.region$bp2[my.region$CHROM == i] = my.region$end[my.region$CHROM == i] + max.now
}
my.snp$bp <- bp

pdf("BayPass_plot.pdf", width=8, height=3)
for (i in unique(my.data$COVARIABLE)) {
  my.plot <- my.data %>% 
    filter(COVARIABLE == i) %>% inner_join(my.snp) %>% select(CHR, bp, BF.dB.)
  my.plot2 <- my.inversion %>% 
    filter(COVARIABLE == i) %>% inner_join(my.region) %>% select(CHROM, bp1, bp2, BF.dB.)
  my.q <- my.pod %>% 
    filter(COVARIABLE == i) %>% summarize(q99=quantile(BF.dB.,probs=0.99)) %>% pull() # Significance threshold
  p <- ggplot() + theme_classic() +
    geom_point(data=my.plot, mapping=aes(x=bp, y=BF.dB., col=CHR), pch=20, size=1, na.rm=T, show.legend=FALSE) +
    geom_segment(data=my.plot2, mapping=aes(x=bp1,xend=bp2,y=BF.dB.,yend=BF.dB.), col="red", size=2) +
    geom_hline(aes(yintercept=my.q), col="red", alpha=0.7, linetype="dotted") +
    xlab("Chromosome") + ylab(expression(BF[is]*" "*(dB))) +
    scale_colour_manual(values=rep(c("#7FC97F","grey70"),9)) +
    scale_x_continuous(breaks=mid, labels=levels(my.plot$CHR), expand=c(0.02,0.02)) +
    scale_y_continuous(expand=c(0,0), limits=c(-12,35)) +
    theme(axis.line.x=element_blank()) +
    ggtitle(my.cov$V1[i])
  print(p)
}
dev.off()
