## Plot genotype frequency on map
## Kaichi Huang 2018 Nov

library(tidyverse)
library(reshape2)
library(sp)
library(raster)
library(scatterpie)
library(gridExtra)
select=dplyr::select

# Prepare genotype dataset
sample_info <- read_tsv("../sample_info.tsv")
cluster_genotypes <- read_tsv("cluster_genotypes.tsv")
mds_info <- read_tsv("mds_info.tsv")
genotype_info <- sample_info %>% inner_join(cluster_genotypes) %>% inner_join(mds_info) %>% 
  mutate(location=paste0(chromosome,":",start,"-",end)) %>% 
  select(name, pop, eco, genotype, lat, long, location)

# Map
bg_map <- raster("../data_code_map/bg_latlong.tif")
bbox <- extent(c(range(genotype_info$long)+c(-1,1)*0.05,range(genotype_info$lat)+c(-1,1)*0.05))
bg_box <- crop(bg_map,bbox)
bg_df <- rasterToPoints(bg_box)
bg_df <- data.frame(bg_df)
colnames(bg_df) <- c("long","lat","lc")
bg_df$lc <- factor(bg_df$lc)
# Color only two land cover types
lc_col <- rep("gray60",length(unique(bg_df$lc)))
lc_col[c(which(levels(bg_df$lc)==31),which(levels(bg_df$lc)==52))] <- c("burlywood1","yellowgreen") # "Barren Land (Rock/Sand/Clay)" and "Shrub/Scrub" 
pie_col <- c("red","purple","blue")
col_all <- c(lc_col,pie_col)

# Plotting
pdf("Map_plot.pdf",height=8,width=9)
for (inv in unique(genotype_info$location)) {
  pop_group <- genotype_info %>% filter(location == inv) %>% select(-location) %>%
    group_by(pop,lat,long,genotype) %>% summarize(count=n()) %>% spread(genotype, count, fill=0) %>%
    mutate(cluster_0_perc=`0`/(`0`+`1`+`2`), cluster_1_perc=`1`/(`0`+`1`+`2`), cluster_2_perc=`2`/(`0`+`1`+`2`))
  p <- ggplot(data = bg_df, aes(long, lat)) +
    geom_raster(aes(fill=lc),alpha=.9) + theme_bw() +
    geom_scatterpie(aes(x=long, y=lat, group=pop, r=0.01),
                    data=pop_group, cols=c("cluster_0_perc","cluster_1_perc","cluster_2_perc"),
                    color=NA, alpha=.8) +
    scale_fill_manual(values=col_all, name="Land Cover Classification\nand\nCluster") +
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + coord_fixed(ratio=1) +
    xlab("Longitude") + ylab("Latitude") +
    ggtitle(inv)
  print(p)
}
dev.off()