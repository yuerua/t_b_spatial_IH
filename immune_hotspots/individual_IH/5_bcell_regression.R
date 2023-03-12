##############################################################################=
## Compare immune cell percentages in individual IH assigned with B cell%-high/low
##############################################################################=
# library(spdep)
library(EBImage)
library(ggpubr)
library(raster)
library(plyr)
library(dplyr)
library(reshape2)
# library(splancs)
# library(mclust)
save_dir <- "/Users/hzhang/Documents/project/sum/final/plot/clump_l_cl/UD_only"
dir.create(save_dir)

cus_theme <-theme(plot.title=element_text(size=16),
                    axis.text=element_text(size=14),
                    axis.title=element_text(size=16),
                    legend.title = element_text(size = 14),
                    legend.text = element_text(size = 14))

#outside TLS ####
df_cell_cd8_foxp3_bcell <- read.csv("/Users/hzhang/Documents/project/sum/final/clump_l_cl/sum/sum_all_no_join_clump.csv")
df_cell_cd8_foxp3_bcell$X <- NULL
df_cell_cd8_foxp3_bcell$sum_cells <- rowSums(df_cell_cd8_foxp3_bcell[,c(which(colnames(df_cell_cd8_foxp3_bcell)=='cd4'):which(colnames(df_cell_cd8_foxp3_bcell)=='foxp3'),
                                                                        which(colnames(df_cell_cd8_foxp3_bcell)=='cd20'):which(colnames(df_cell_cd8_foxp3_bcell)=='cd79bCoexp'))])

df_cell_cd8_foxp3_bcell$bcells <- rowSums(df_cell_cd8_foxp3_bcell[,which(colnames(df_cell_cd8_foxp3_bcell)=='cd20'):which(colnames(df_cell_cd8_foxp3_bcell)=='cd79bCoexp')])

df_cell_cd8_foxp3_bcell_cl <- subset(df_cell_cd8_foxp3_bcell, Hotspots == "Cancer-immune")
df_cell_cd8_foxp3_bcell_l <- subset(df_cell_cd8_foxp3_bcell, Hotspots == "Immune")

df <- df_cell_cd8_foxp3_bcell_l %>% group_by(slide) %>% summarise(n = n())
mean(df$n)
sd(df$n)

df_cell_cd8_foxp3_bcell_cl$Hotspots <- 'Intratumoral IH'
df_cell_cd8_foxp3_bcell_l$Hotspots <- 'Peritumoral IH'


select_col <- c('slide', 'Hotspots', 'Region', 'clump_id', 'cd8', 'cd4', 'foxp3', 'cd20', 'cd20cxcr5', 'cd79bCoexp', 'sum_cells', 'bcells')

df_cl_l <- rbind(df_cell_cd8_foxp3_bcell_l[, select_col], df_cell_cd8_foxp3_bcell_cl[, select_col])


df_cl_l <- subset(df_cl_l, Region == "UD")


# df_cl_l[,c('cd8', 'cd4', 'foxp3', 'cd20', 'cd20cxcr5', 'cd79bCoexp')] <- df_cl_l[,c('cd8', 'cd4', 'foxp3', 'cd20', 'cd20cxcr5', 'cd79bCoexp')] * df_cl_l$sum_cells
df_cl_l <- df_cl_l %>% group_by(slide, Hotspots, clump_id) %>% summarise(cd8 = sum(cd8, na.rm = TRUE), cd4 = sum(cd4, na.rm = TRUE),
                                                                                 foxp3 = sum(foxp3, na.rm=TRUE), cd20 = sum(cd20, na.rm=TRUE), 
                                                                                 cd20cxcr5 = sum(cd20cxcr5, na.rm = TRUE), cd79bCoexp = sum(cd79bCoexp, na.rm=TRUE),
                                        
                                                                                 sum_cells = sum(sum_cells, na.rm = TRUE),
                                                                                  area = n())


df_cl_l$bcells <- rowSums(df_cl_l[,c('cd20', 'cd20cxcr5', 'cd79bCoexp')]) / df_cl_l$sum_cells

df_cl_l[,c('cd8', 'cd4', 'foxp3', 'cd20', 'cd20cxcr5', 'cd79bCoexp')] <- df_cl_l[,c('cd8', 'cd4', 'foxp3', 'cd20', 'cd20cxcr5', 'cd79bCoexp')] / df_cl_l$sum_cells


df_cl_l$cd8_foxp3 <- df_cl_l$cd8 / df_cl_l$foxp3
df_cl_l$cd4_foxp3 <- df_cl_l$cd4 / df_cl_l$foxp3

#CD8/FOXP3 between intra and peritumoral#====
df_cl_l$cd8_foxp3 <- as.numeric(df_cl_l$cd8_foxp3)

df <- subset(df_cl_l, Hotspots == "Intratumoral IH")
df <- subset(df, is.finite(cd8_foxp3))
sd(df$cd8_foxp3, na.rm = T)
mean(df$cd8_foxp3, na.rm = T)

df <- subset(df_cl_l, Hotspots == "Peritumoral IH")
df <- subset(df, is.finite(cd8_foxp3))
sd(df$cd8_foxp3, na.rm = T)
mean(df$cd8_foxp3, na.rm = T)

p <- ggboxplot(subset(df_cl_l, !is.na(cd8_foxp3)), y = 'cd8_foxp3', x = 'Hotspots', add = 'jitter', color = 'Hotspots', palette = c("firebrick4", "tan3")) + 
  stat_compare_means(label = "p.format", label.x = 1.4, paired = F, label.y = 50) + 
  xlab('') + ylab('') +
  cus_theme

pdf("/Users/hzhang/Documents/project/sum/final/plot/clump_l_cl/UD_only/cd8_foxp3_intra_peri.pdf", width = 6, height = 6)
print(p)
dev.off()

#cutoff: c/cl median
#cutoff: all hs median
df_cl_l_1 <- subset(df_cl_l, Hotspots == 'Intratumoral IH')
df_cl_l_2 <- subset(df_cl_l, Hotspots == 'Peritumoral IH')

# 
ggscatter(df_cl_l_2, x = 'cd8', y = 'foxp3', add = "reg.line",
          conf.int = TRUE,
          add.params = list(color = "gray40",
                            fill = "lightgray"),
          xlab = 'CD8+ T cells / CD4+FOXP3+ T cells', ylab= "% of B cells")+
  stat_cor(method = "pearson") +
  theme(plot.title=element_text(size=14),
        axis.text=element_text(size=7),
        axis.title=element_text(size=16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))
# 
# 
# 

#assign high/low for all IH together ####
df_cl_l_1$cd8_foxp3_class <- ifelse(df_cl_l_1$cd8_foxp3 >= median(df_cl_l_1$cd8_foxp3, na.rm = T), "High", "Low")
df_cl_l_2$cd8_foxp3_class <- ifelse(df_cl_l_2$cd8_foxp3 >= median(df_cl_l_2$cd8_foxp3, na.rm = T), "High", "Low")

df_cl_l_1$cd4_foxp3_class <- ifelse(df_cl_l_1$cd4_foxp3 >= median(df_cl_l_1$cd4_foxp3, na.rm = T), "High", "Low")
df_cl_l_2$cd4_foxp3_class <- ifelse(df_cl_l_2$cd4_foxp3 >= median(df_cl_l_2$cd4_foxp3, na.rm = T), "High", "Low")

df_cl_l_1$foxp3_class <- ifelse(df_cl_l_1$foxp3 >= median(df_cl_l_1$foxp3, na.rm = T), "High", "Low")
df_cl_l_2$foxp3_class <- ifelse(df_cl_l_2$foxp3 >= median(df_cl_l_2$foxp3, na.rm = T), "High", "Low")

df_cl_l_1$cd8_class <- ifelse(df_cl_l_1$cd8 >= median(df_cl_l_1$cd8, na.rm = T), "High", "Low")
df_cl_l_2$cd8_class <- ifelse(df_cl_l_2$cd8 >= median(df_cl_l_2$cd8, na.rm = T), "High", "Low")

df_cl_l_1$bcell_class <- ifelse(df_cl_l_1$bcells >= median(df_cl_l_1$bcells, na.rm = T), "High", "Low")
df_cl_l_2$bcell_class <- ifelse(df_cl_l_2$bcells >= median(df_cl_l_2$bcells, na.rm = T), "High", "Low")

df_cl_l_3 <- rbind(df_cl_l_1, df_cl_l_2)

df_cl_l_3$cd8_foxp3_class <- factor(df_cl_l_3$cd8_foxp3_class, levels = c("High", "Low"))
df_cl_l_3$cd4_foxp3_class <- factor(df_cl_l_3$cd4_foxp3_class, levels = c("High", "Low"))
df_cl_l_3$foxp3_class <- factor(df_cl_l_3$foxp3_class, levels = c("High", "Low"))
df_cl_l_3$cd8_class <- factor(df_cl_l_3$cd8_class, levels = c("High", "Low"))
df_cl_l_3$bcell_class <- factor(df_cl_l_3$bcell_class, levels = c("High", "Low"))

#write.csv(df_cl_l_3, "/Users/hzhang/Documents/report/lusc_b/plots_v4/fig4_cl_l_by_clump/df_cl_l_3_UD_only.csv")

#=========================================#
#assign high/low at the slide level ####
save_dir <- "/Users/hzhang/Documents/project/sum/final/plot/clump_l_cl/UD_only_median_per_slide"
dir.create(save_dir)

df_cl_l_1 <- subset(df_cl_l, Hotspots == 'Intratumoral IH')
df_cl_l_2 <- subset(df_cl_l, Hotspots == 'Peritumoral IH')

df_cl_l_1 <- df_cl_l_1 %>% group_by(slide) %>% mutate(cd8_m = median(cd8, na.rm = T), 
                                                      cd4_m = median(cd4, na.rm = T), 
                                                      foxp3_m = median(foxp3, na.rm = T), 
                                                      cd20_m = median(cd20, na.rm = T), 
                                                      cd20cxcr5_m = median(cd20cxcr5, na.rm = T), 
                                                      cd79bCoexp_m = median(cd79bCoexp, na.rm = T), 
                                                      cd8_foxp3_m = median(cd8_foxp3, na.rm = T), 
                                                      cd4_foxp3_m = median(cd4_foxp3, na.rm = T), 
                                                      bcells_m = median(bcells, na.rm = T))

df_cl_l_2 <- df_cl_l_2 %>% group_by(slide) %>% mutate(cd8_m = median(cd8, na.rm = T), 
                                                      cd4_m = median(cd4, na.rm = T), 
                                                      foxp3_m = median(foxp3, na.rm = T), 
                                                      cd20_m = median(cd20, na.rm = T), 
                                                      cd20cxcr5_m = median(cd20cxcr5, na.rm = T), 
                                                      cd79bCoexp_m = median(cd79bCoexp, na.rm = T), 
                                                      cd8_foxp3_m = median(cd8_foxp3, na.rm = T), 
                                                      cd4_foxp3_m = median(cd4_foxp3, na.rm = T), 
                                                      bcells_m = median(bcells, na.rm = T))

cell_class <- c('cd8_foxp3', 'cd4_foxp3', 'foxp3', 'cd8', 'bcells')

assign_high_low <- function(df, cell_class){
  for (c in cell_class){
    df[paste0(c, "_class")] <- ifelse(df[[c]] >= df[[paste0(c, "_m")]], "High", "Low")
    df[paste0(c, "_class")] <-  factor(df[[paste0(c, "_class")]], levels = c("High", "Low"))
  }
  return(df)
}

df_cl_l_1 <- assign_high_low(df_cl_l_1, cell_class)
df_cl_l_2 <- assign_high_low(df_cl_l_2, cell_class)


df_cl_l_3 <- rbind(df_cl_l_1, df_cl_l_2)

write.csv(df_cl_l_3, file.path(save_dir, "df_cl_l_3.csv"))
#=========================================#

#assign high/low to each slide, then apply the same label to each IH within the slide####
save_dir <- "/Users/hzhang/Documents/project/sum/final/plot/clump_l_cl/UD_only_median_slide_level"
dir.create(save_dir)

df_cl_l_1 <- subset(df_cl_l, Hotspots == 'Intratumoral IH')
df_cl_l_2 <- subset(df_cl_l, Hotspots == 'Peritumoral IH')

df_cl_l_1 <- df_cl_l_1 %>% group_by(slide) %>% mutate(cd8_m = median(cd8, na.rm = T), 
                                                      cd4_m = median(cd4, na.rm = T), 
                                                      foxp3_m = median(foxp3, na.rm = T), 
                                                      cd20_m = median(cd20, na.rm = T), 
                                                      cd20cxcr5_m = median(cd20cxcr5, na.rm = T), 
                                                      cd79bCoexp_m = median(cd79bCoexp, na.rm = T), 
                                                      cd8_foxp3_m = median(cd8_foxp3, na.rm = T), 
                                                      cd4_foxp3_m = median(cd4_foxp3, na.rm = T), 
                                                      bcells_m = median(bcells, na.rm = T))
df_cl_l_2 <- df_cl_l_2 %>% group_by(slide) %>% mutate(cd8_m = median(cd8, na.rm = T), 
                                                      cd4_m = median(cd4, na.rm = T), 
                                                      foxp3_m = median(foxp3, na.rm = T), 
                                                      cd20_m = median(cd20, na.rm = T), 
                                                      cd20cxcr5_m = median(cd20cxcr5, na.rm = T), 
                                                      cd79bCoexp_m = median(cd79bCoexp, na.rm = T), 
                                                      cd8_foxp3_m = median(cd8_foxp3, na.rm = T), 
                                                      cd4_foxp3_m = median(cd4_foxp3, na.rm = T), 
                                                      bcells_m = median(bcells, na.rm = T))

cell_class <- c('cd8_foxp3', 'cd4_foxp3', 'foxp3', 'cd8', 'bcells')

assign_high_low <- function(df, cell_class){
  for (c in cell_class){
    df[paste0(c, "_class")] <- ifelse(df[[paste0(c, "_m")]] >= median(df[[paste0(c, "_m")]], na.rm = T), "High", "Low")
    df[paste0(c, "_class")] <-  factor(df[[paste0(c, "_class")]], levels = c("High", "Low"))
  }
  return(df)
}


df_cl_l_1 <- assign_high_low(df_cl_l_1, cell_class)


df_cl_l_2 <- assign_high_low(df_cl_l_2, cell_class)

df_cl_l_3 <- rbind(df_cl_l_1, df_cl_l_2)

write.csv(df_cl_l_3, file.path(save_dir, "df_cl_l_3_slide_median.csv"))

#cd8/foxp3 high low
cell_class <- c("cd4", "cd20", "cd20cxcr5", "cd79bCoexp", "bcells")
# cell_class <- c("cd20cxcr5")
cell_label <- c("% of CD4+FOXP3- T cells", "% of CD20+CXCR5- B cells", "% of CD20+CXCR5+ B cells", "% of CD79b+ B cells",  "% of B cells")
# cell_label <- c("% of CD20+CXCR5+ B cells")
save_path <- file.path(save_dir, "cd8_foxp3_low_high_cl_l_median")
dir.create(save_path)
for (i in 1:length(cell_class)){
  print(cell_class[i])
  pdf(file.path(save_path, paste0(cell_class[i], ".pdf")))
  
  p <- ggboxplot(subset(df_cl_l_3, !is.na(cd8_foxp3_class)), x = 'cd8_foxp3_class', y = cell_class[i], add = 'jitter', color = 'cd8_foxp3_class', palette = c("firebrick4", "tan3"), facet.by = "Hotspots") + 
    stat_compare_means(label = "p.format", label.x = 1.3, paired = F) + #label.y = 0.35
    xlab('') + ylab(cell_label[i]) #+ ylim(c(0,0.4))
    cus_theme
  
  print(p)
  dev.off()
  
}

#cd4/foxp3 high low
cell_class <- c("cd8", "cd20", "cd20cxcr5", "cd79bCoexp", "bcells")
cell_label <- c("% of CD8+ T cells", "% of CD20+CXCR5- B cells", "% of CD20+CXCR5+ B cells", "% of CD79b+ B cells", "% of B cells")
save_path <- file.path(save_dir, "cd4_foxp3_low_high_cl_l_median")
dir.create(save_path)
for (i in 1:length(cell_class)){
  print(cell_class[i])
  pdf(file.path(save_path, paste0(cell_class[i], ".pdf")))
  
  p <- ggboxplot(subset(df_cl_l_3, !is.na(cd4_foxp3_class)), x = 'cd4_foxp3_class', y = cell_class[i], add = 'jitter', color = 'cd4_foxp3_class', palette = c("firebrick4", "tan3"), facet.by = "Hotspots") + 
    stat_compare_means(label = "p.format", label.x = 1.3, paired = F) + #, label.y = 1 
    xlab('') + ylab(cell_label[i]) +
    cus_theme
  
  print(p)
  dev.off()
  
}

#foxp3 perc high low
cell_class <- c("cd8", "cd20", "cd20cxcr5", "cd79bCoexp", "bcells", "cd4")
cell_label <- c("% of CD8+ T cells", "% of CD20+CXCR5- B cells", "% of CD20+CXCR5+ B cells", "% of CD79b+ B cells", "% of B cells", "% of CD4+FOXP3- cells")
save_path <- file.path(save_dir, "foxp3_low_high_cl_l_median")
dir.create(save_path)
for (i in 1:length(cell_class)){
  print(cell_class[i])
  pdf(file.path(save_path, paste0(cell_class[i], ".pdf")))
  
  p <- ggboxplot(subset(df_cl_l_3, !is.na(foxp3_class)), x = 'foxp3_class', y = cell_class[i], add = 'jitter', color = 'foxp3_class', palette = c("firebrick4", "tan3"), facet.by = "Hotspots") + 
    stat_compare_means(label = "p.format", label.x = 1.3, paired = F) + #, label.y = 1
    xlab('') + ylab(cell_label[i]) +
    cus_theme
  
  print(p)
  dev.off()
  
}

#cd8 perc high low
cell_class <- c("foxp3", "cd20", "cd20cxcr5", "cd79bCoexp", "bcells", "cd4")
cell_label <- c("% of CD4+FOXP3+ T cells", "% of CD20+CXCR5- B cells", "% of CD20+CXCR5+ B cells", "% of CD79b+ B cells", "% of B cells", "% of CD4+FOXP3- cells")
save_path <- file.path(save_dir, "cd8_low_high_cl_l_median")
dir.create(save_path)
for (i in 1:length(cell_class)){
  print(cell_class[i])
  pdf(file.path(save_path, paste0(cell_class[i], ".pdf")))
  
  p <- ggboxplot(subset(df_cl_l_3, !is.na(cd8_class)), x = 'cd8_class', y = cell_class[i], add = 'jitter', color = 'cd8_class', palette = c("firebrick4", "tan3"), facet.by = "Hotspots") + 
    stat_compare_means(label = "p.format", label.x = 1.3, paired = F) + #, label.y = 1
    xlab('') + ylab(cell_label[i]) +
    cus_theme
  
  print(p)
  dev.off()
  
}

#bcell perc high low
cell_class <- c("foxp3", "cd20", "cd20cxcr5", "cd79bCoexp", "cd8", "cd4", "cd4_foxp3", "cd8_foxp3")
cell_label <- c("% of CD4+FOXP3+ T cells", "% of CD20+CXCR5- B cells", "% of CD20+CXCR5+ B cells", "% of CD79b+ B cells", "% of CD8+ T cells", "% of CD4+FOXP3- cells",
                "CD4+FOXP3- / CD4+FOXP3+", "CD8+ / CD4+FOXP3+")
save_path <- file.path(save_dir, "bcell_low_high_cl_l_median")
dir.create(save_path)
for (i in 1:length(cell_class)){
  print(cell_class[i])
  pdf(file.path(save_path, paste0(cell_class[i], ".pdf")))
  
  p <- ggboxplot(subset(df_cl_l_3, !is.na(bcells_class)), x = 'bcells_class', y = cell_class[i], add = 'jitter', color = 'bcells_class', palette = c("firebrick4", "tan3"), facet.by = "Hotspots") + 
    stat_compare_means(label = "p.format", label.x = 1.3, paired = F) + 
    xlab('') + ylab(cell_label[i]) +
    cus_theme
  
  print(p)
  dev.off()
  
}


#statis
cell_class <- c("cd8", "cd20", "cd20cxcr5", "cd79bCoexp", "bcells", "cd4")
df <- subset(df_cl_l_3, !is.na(foxp3_class)) %>% group_by(Hotspots, foxp3_class) %>% summarise_at(cell_class, mean)
write.csv(df, file.path(save_dir, "foxp3_low_high_cl_l_median/mean.csv"))

cell_class <- c("cd4", "cd20", "cd20cxcr5", "cd79bCoexp", "bcells")
df <- subset(df_cl_l_3, !is.na(cd8_foxp3_class)) %>% group_by(Hotspots, cd8_foxp3_class) %>% summarise_at(cell_class, mean)
write.csv(df, file.path(save_dir, "cd8_foxp3_low_high_cl_l_median/mean.csv"))

cell_class <- c("cd8", "cd20", "cd20cxcr5", "cd79bCoexp", "bcells")
df <- subset(df_cl_l_3, !is.na(cd4_foxp3_class)) %>% group_by(Hotspots, cd4_foxp3_class) %>% summarise_at(cell_class, mean)
write.csv(df, file.path(save_dir, "cd4_foxp3_low_high_cl_l_median/mean.csv"))

cell_class <- c("foxp3", "cd20", "cd20cxcr5", "cd79bCoexp", "bcells", "cd4")
df <- subset(df_cl_l_3, !is.na(cd8_class)) %>% group_by(Hotspots, cd8_class) %>% summarise_at(cell_class, mean)
write.csv(df, file.path(save_dir, "cd4_foxp3_low_high_cl_l_median/mean.csv"))


cell_class <- c("foxp3", "cd20", "cd20cxcr5", "cd79bCoexp", "cd8", "cd4", "cd4_foxp3", "cd8_foxp3")
df <- subset(df_cl_l_3, !is.na(bcells_class)) %>% group_by(Hotspots, bcells_class) %>% summarise_at(cell_class, mean)
write.csv(df, file.path(save_dir, "bcell_low_high_cl_l_median/mean.csv"))

#outside TLS by density####

df_cell_cd8_foxp3_bcell_cl <- read.csv("/Users/hzhang/Documents/report/lusc_b/scripts/ITLR/results/clump_cl/df_cell_cd8_foxp3_bcell_by_density.csv")
df_cell_cd8_foxp3_bcell_l <- read.csv("/Users/hzhang/Documents/report/lusc_b/scripts/ITLR/results/clump_l/df_cell_cd8_foxp3_bcell_by_density.csv")

df <- df_cell_cd8_foxp3_bcell_l %>% group_by(slide) %>% summarise(n = n())
mean(df$n)
sd(df$n)

df <- subset(df_cell_cd8_foxp3_bcell_cl, Region == 'TLS') %>% group_by(slide) %>% summarise(n = n())
mean(df$n)
sd(df$n)



df_cell_cd8_foxp3_bcell_cl$Hotspots <- 'Intratumoral IH'
df_cell_cd8_foxp3_bcell_l$Hotspots <- 'Peritumoral IH'

select_col <- c('slide', 'Hotspots', 'Region', 'clump_id', 'cd8', 'cd4', 'foxp3', 'cd20', 'cd20cxcr5', 'cd79bCoexp', 'Adjacent', 'Distal', 'Adjacent_cd8', 'Distal_cd8', 'Adjacent_foxp3', 'Distal_foxp3', 'sum_cells', 'bcells')

df_cl_l <- rbind(df_cell_cd8_foxp3_bcell_l[, select_col], df_cell_cd8_foxp3_bcell_cl[, select_col])


df_cl_l <- subset(df_cl_l, Region == "UD")


df_cl_l[,c('cd8', 'cd4', 'foxp3', 'cd20', 'cd20cxcr5', 'cd79bCoexp')] <- df_cl_l[,c('cd8', 'cd4', 'foxp3', 'cd20', 'cd20cxcr5', 'cd79bCoexp')]/ (50*50) * 10^6

# df_cl_l <- df_cl_l %>% group_by(slide, Hotspots, clump_id) %>% summarise(cd8 = sum(cd8, na.rm = TRUE), cd4 = sum(cd4, na.rm = TRUE),
#                                                                          foxp3 = sum(foxp3, na.rm=TRUE), cd20 = sum(cd20, na.rm=TRUE), 
#                                                                          cd20cxcr5 = sum(cd20cxcr5, na.rm = TRUE), cd79bCoexp = sum(cd79bCoexp, na.rm=TRUE),
#                                                                          Adjacent = mean(Adjacent, na.rm = T), Distal = mean(Distal, na.rm = T),
#                                                                          Adjacent_cd8 = mean(Adjacent_cd8), 
#                                                                          Distal_cd8 = mean(Distal_cd8),
#                                                                          Adjacent_foxp3 = mean(Adjacent_foxp3),
#                                                                          Distal_foxp3 = mean(Distal_foxp3), 
#                                                                          sum_cells = sum(sum_cells, na.rm = TRUE),
#                                                                          area = n())



df_cl_l$cd8_foxp3 <- df_cl_l$cd8 / df_cl_l$foxp3
df_cl_l$cd4_foxp3 <- df_cl_l$cd4 / df_cl_l$foxp3

#cutoff: c/cl median
#cutoff: all hs median
df_cl_l_1 <- subset(df_cl_l, Hotspots == 'Intratumoral IH')
df_cl_l_2 <- subset(df_cl_l, Hotspots == 'Peritumoral IH')

# 
ggscatter(df_cl_l_2, x = 'cd8', y = 'foxp3', add = "reg.line",
          conf.int = TRUE,
          add.params = list(color = "gray40",
                            fill = "lightgray"),
          xlab = 'CD8+ T cells / CD4+FOXP3+ T cells', ylab= "% of B cells")+
  stat_cor(method = "pearson") +
  theme(plot.title=element_text(size=14),
        axis.text=element_text(size=7),
        axis.title=element_text(size=16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))
# 
# 
# 


df_cl_l_1$cd8_foxp3_class <- ifelse(df_cl_l_1$cd8_foxp3 >= median(df_cl_l_1$cd8_foxp3, na.rm = T), "High", "Low")
df_cl_l_2$cd8_foxp3_class <- ifelse(df_cl_l_2$cd8_foxp3 >= median(df_cl_l_2$cd8_foxp3, na.rm = T), "High", "Low")

df_cl_l_1$cd4_foxp3_class <- ifelse(df_cl_l_1$cd4_foxp3 >= median(df_cl_l_1$cd4_foxp3, na.rm = T), "High", "Low")
df_cl_l_2$cd4_foxp3_class <- ifelse(df_cl_l_2$cd4_foxp3 >= median(df_cl_l_2$cd4_foxp3, na.rm = T), "High", "Low")

df_cl_l_1$foxp3_class <- ifelse(df_cl_l_1$foxp3 >= median(df_cl_l_1$foxp3, na.rm = T), "High", "Low")
df_cl_l_2$foxp3_class <- ifelse(df_cl_l_2$foxp3 >= median(df_cl_l_2$foxp3, na.rm = T), "High", "Low")

df_cl_l_1$cd8_class <- ifelse(df_cl_l_1$cd8 >= median(df_cl_l_1$cd8, na.rm = T), "High", "Low")
df_cl_l_2$cd8_class <- ifelse(df_cl_l_2$cd8 >= median(df_cl_l_2$cd8, na.rm = T), "High", "Low")

df_cl_l_1$bcell_class <- ifelse(df_cl_l_1$bcells >= median(df_cl_l_1$bcells, na.rm = T), "High", "Low")
df_cl_l_2$bcell_class <- ifelse(df_cl_l_2$bcells >= median(df_cl_l_2$bcells, na.rm = T), "High", "Low")

df_cl_l_3 <- rbind(df_cl_l_1, df_cl_l_2)

df_cl_l_3$cd8_foxp3_class <- factor(df_cl_l_3$cd8_foxp3_class, levels = c("High", "Low"))
df_cl_l_3$cd4_foxp3_class <- factor(df_cl_l_3$cd4_foxp3_class, levels = c("High", "Low"))
df_cl_l_3$foxp3_class <- factor(df_cl_l_3$foxp3_class, levels = c("High", "Low"))
df_cl_l_3$cd8_class <- factor(df_cl_l_3$cd8_class, levels = c("High", "Low"))
df_cl_l_3$bcell_class <- factor(df_cl_l_3$bcell_class, levels = c("High", "Low"))

#write.csv(df_cl_l_3, "/Users/hzhang/Documents/report/lusc_b/plots_v4/fig4_cl_l_by_clump/df_cl_l_3_UD_only.csv")

#cd8/foxp3 high low
cell_class <- c("cd4", "cd20", "cd20cxcr5", "cd79bCoexp", "bcells")
cell_label <- c("CD4+FOXP3- T cells/mm^2", "CD20+CXCR5- B cells/mm^2", "CD20+CXCR5+ B cells/mm^2", "CD79b+ B cells/mm^2",  "B cells/mm^2")

for (i in 1:length(cell_class)){
  print(cell_class[i])
  pdf(file.path("/Users/hzhang/Documents/report/lusc_b/plots_v4/fig4_cl_l_by_clump/UD_only/by_density/cd8_foxp3_low_high_cl_l_median/", paste0(cell_class[i], ".pdf")))
  
  p <- ggboxplot(subset(df_cl_l_3, !is.na(cd8_foxp3_class)), x = 'cd8_foxp3_class', y = cell_class[i], add = 'jitter', color = 'cd8_foxp3_class', palette = c("firebrick4", "tan3"), facet.by = "Hotspots") + 
    stat_compare_means(label = "p.format", label.x = 1.3, paired = F) + 
    xlab('') + ylab(cell_label[i]) +
    cus_theme
  
  print(p)
  dev.off()
  
}

#cd4/foxp3 high low
cell_class <- c("cd8", "cd20", "cd20cxcr5", "cd79bCoexp", "bcells")
cell_label <- c("CD8+ T cells/mm^2", "CD20+CXCR5- B cells/mm^2", "CD20+CXCR5+ B cells/mm^2", "CD79b+ B cells/mm^2", "B cells/mm^2")

for (i in 1:length(cell_class)){
  print(cell_class[i])
  pdf(file.path("/Users/hzhang/Documents/report/lusc_b/plots_v4/fig4_cl_l_by_clump/UD_only/by_density/cd4_foxp3_low_high_cl_l_median/", paste0(cell_class[i], ".pdf")))
  
  p <- ggboxplot(subset(df_cl_l_3, !is.na(cd4_foxp3_class)), x = 'cd4_foxp3_class', y = cell_class[i], add = 'jitter', color = 'cd4_foxp3_class', palette = c("firebrick4", "tan3"), facet.by = "Hotspots") + 
    stat_compare_means(label = "p.format", label.x = 1.3, paired = F) + 
    xlab('') + ylab(cell_label[i]) +
    cus_theme
  
  print(p)
  dev.off()
  
}

#foxp3 perc high low
cell_class <- c("cd8", "cd20", "cd20cxcr5", "cd79bCoexp", "bcells", "cd4")
cell_label <- c("CD8+ T cells/mm^2", "CD20+CXCR5- B cells/mm^2", "CD20+CXCR5+ B cells/mm^2", "CD79b+ B cells/mm^2", "B cells/mm^2", "CD4+FOXP3- cells/mm^2")

for (i in 1:length(cell_class)){
  print(cell_class[i])
  pdf(file.path("/Users/hzhang/Documents/report/lusc_b/plots_v4/fig4_cl_l_by_clump/UD_only/by_density/foxp3_low_high_cl_l_median/", paste0(cell_class[i], ".pdf")))
  
  p <- ggboxplot(subset(df_cl_l_3, !is.na(foxp3_class)), x = 'foxp3_class', y = cell_class[i], add = 'jitter', color = 'foxp3_class', palette = c("firebrick4", "tan3"), facet.by = "Hotspots") + 
    stat_compare_means(label = "p.format", label.x = 1.3, paired = F) + 
    xlab('') + ylab(cell_label[i]) +
    cus_theme
  
  print(p)
  dev.off()
  
}

#cd8 perc high low
cell_class <- c("foxp3", "cd20", "cd20cxcr5", "cd79bCoexp", "bcells", "cd4")
cell_label <- c("CD4+FOXP3+ T cells/mm^2", "CD20+CXCR5- B cells/mm^2", "CD20+CXCR5+ B cells/mm^2", "CD79b+ B cells/mm^2", "B cells/mm^2", "CD4+FOXP3- cells/mm^2")

for (i in 1:length(cell_class)){
  print(cell_class[i])
  pdf(file.path("/Users/hzhang/Documents/report/lusc_b/plots_v4/fig4_cl_l_by_clump/UD_only/by_density/cd8_low_high_cl_l_median/", paste0(cell_class[i], ".pdf")))
  
  p <- ggboxplot(subset(df_cl_l_3, !is.na(cd8_class)), x = 'cd8_class', y = cell_class[i], add = 'jitter', color = 'cd8_class', palette = c("firebrick4", "tan3"), facet.by = "Hotspots") + 
    stat_compare_means(label = "p.format", label.x = 1.3, paired = F) + 
    xlab('') + ylab(cell_label[i]) +
    cus_theme
  
  print(p)
  dev.off()
  
}

#bcell perc high low
cell_class <- c("foxp3", "cd20", "cd20cxcr5", "cd79bCoexp", "cd8", "cd4", "cd4_foxp3", "cd8_foxp3")
cell_label <- c("CD4+FOXP3+ T cells/mm^2", "CD20+CXCR5- B cells/mm^2", "CD20+CXCR5+ B cells/mm^2", "CD79b+ B cells/mm^2", "CD8+ T cells/mm^2", "CD4+FOXP3- cells/mm^2",
                "CD4+FOXP3- / CD4+FOXP3+", "CD8+ / CD4+FOXP3+")

for (i in 1:length(cell_class)){
  print(cell_class[i])
  pdf(file.path("/Users/hzhang/Documents/report/lusc_b/plots_v4/fig4_cl_l_by_clump/UD_only/by_density/bcell_low_high_cl_l_median/", paste0(cell_class[i], ".pdf")))
  
  p <- ggboxplot(subset(df_cl_l_3, !is.na(bcell_class)), x = 'bcell_class', y = cell_class[i], add = 'jitter', color = 'bcell_class', palette = c("firebrick4", "tan3"), facet.by = "Hotspots") + 
    stat_compare_means(label = "p.format", label.x = 1.3, paired = F) + 
    xlab('') + ylab(cell_label[i]) +
    cus_theme
  
  print(p)
  dev.off()
  
}


#statis
cell_class <- c("cd8", "cd20", "cd20cxcr5", "cd79bCoexp", "bcells", "cd4")
df <- subset(df_cl_l_3, !is.na(foxp3_class)) %>% group_by(Hotspots, foxp3_class) %>% summarise_at(cell_class, mean)
write.csv(df, "/Users/hzhang/Documents/report/lusc_b/plots_v4/fig4_cl_l_by_clump/UD_only/by_density/foxp3_low_high_cl_l_median/mean.csv")

cell_class <- c("cd4", "cd20", "cd20cxcr5", "cd79bCoexp", "bcells")
df <- subset(df_cl_l_3, !is.na(cd8_foxp3_class)) %>% group_by(Hotspots, cd8_foxp3_class) %>% summarise_at(cell_class, mean)
write.csv(df, "/Users/hzhang/Documents/report/lusc_b/plots_v4/fig4_cl_l_by_clump/UD_only/by_density/cd8_foxp3_low_high_cl_l_median/mean.csv")

cell_class <- c("cd8", "cd20", "cd20cxcr5", "cd79bCoexp", "bcells")
df <- subset(df_cl_l_3, !is.na(cd4_foxp3_class)) %>% group_by(Hotspots, cd4_foxp3_class) %>% summarise_at(cell_class, mean)
write.csv(df, "/Users/hzhang/Documents/report/lusc_b/plots_v4/fig4_cl_l_by_clump/UD_only/by_density/cd4_foxp3_low_high_cl_l_median/mean.csv")

cell_class <- c("foxp3", "cd20", "cd20cxcr5", "cd79bCoexp", "bcells", "cd4")
df <- subset(df_cl_l_3, !is.na(cd8_class)) %>% group_by(Hotspots, cd8_class) %>% summarise_at(cell_class, mean)
write.csv(df, "/Users/hzhang/Documents/report/lusc_b/plots_v4/fig4_cl_l_by_clump/UD_only/by_density/cd4_foxp3_low_high_cl_l_median/mean.csv")

cell_class <- c("foxp3", "cd20", "cd20cxcr5", "cd79bCoexp", "cd8", "cd4", "cd4_foxp3", "cd8_foxp3")
df <- subset(df_cl_l_3, !is.na(bcell_class)) %>% group_by(Hotspots, bcell_class) %>% summarise_at(cell_class, mean)
write.csv(df, "/Users/hzhang/Documents/report/lusc_b/plots_v4/fig4_cl_l_by_clump/UD_only/by_density/bcell_low_high_cl_l_median/mean.csv")

#compare the slide_level mean, median#====
df_slide_level_mean <- read.csv("/Users/hzhang/Documents/project/sum/final/plot/clump_l_cl/UD_only_median_slide_level/slide_mean/df_cl_l_3_slide_mean.csv")
df_slide_level_median <- read.csv("/Users/hzhang/Documents/project/sum/final/plot/clump_l_cl/UD_only_median_slide_level/slide_median/df_cl_l_3_slide_median.csv")

df_1 <- df_slide_level_mean %>% group_by(slide, Hotspots) %>% summarise(n_mean = unique(foxp3_class))

df_2 <- df_slide_level_median %>% group_by(slide, Hotspots) %>% summarise(n_median = unique(foxp3_class))

merge(df_1, df_2, by = c('slide', 'Hotspots'))
# slide        Hotspots     n_mean  n_median
# 51464 Intratumoral IH   High     High
# 51464  Peritumoral IH    Low     High
# 51746 Intratumoral IH    Low      Low
# 51746  Peritumoral IH   High      Low
# so 51746 high foxp3 in peritumoral IH, 51464 low foxp3 in peritumoral IH
#======================#

#nearby cd8/foxp3 ####
df_cell_cd8_foxp3_bcell_cl <- read.csv("/Users/hzhang/Documents/report/lusc_b/scripts/ITLR/results/clump_cl/df_cell_cd8_foxp3_bcell_v2.csv")
df_cell_cd8_foxp3_bcell_l <- read.csv("/Users/hzhang/Documents/report/lusc_b/scripts/ITLR/results/clump_l/df_cell_cd8_foxp3_bcell_v1.csv")

df <- df_cell_cd8_foxp3_bcell_l %>% group_by(slide) %>% summarise(n = n())
mean(df$n)
sd(df$n)

df_cell_cd8_foxp3_bcell_cl$Hotspots <- 'Intratumoral IH'
df_cell_cd8_foxp3_bcell_l$Hotspots <- 'Peritumoral IH'

select_col <- c('slide', 'Hotspots', 'Region', 'clump_id', 'cd8', 'cd4', 'foxp3', 'cd20', 'cd20cxcr5', 'cd79bCoexp', 'Adjacent', 'Distal', 'Adjacent_cd8', 'Distal_cd8', 'Adjacent_foxp3', 'Distal_foxp3', 'sum_cells', 'bcells')

df_cl_l <- rbind(df_cell_cd8_foxp3_bcell_l[, select_col], df_cell_cd8_foxp3_bcell_cl[, select_col])

mean(table(df_cl_l$Hotspots, df_cl_l$slide)[2,])

df_cl_l[,c('cd8', 'cd4', 'foxp3', 'cd20', 'cd20cxcr5', 'cd79bCoexp')] <- df_cl_l[,c('cd8', 'cd4', 'foxp3', 'cd20', 'cd20cxcr5', 'cd79bCoexp')] * df_cl_l$sum_cells
df_cl_l <- df_cl_l %>% group_by(slide, Hotspots, clump_id) %>% summarise(cd8 = sum(cd8, na.rm = TRUE), cd4 = sum(cd4, na.rm = TRUE),
                                                                         foxp3 = sum(foxp3, na.rm=TRUE), cd20 = sum(cd20, na.rm=TRUE), 
                                                                         cd20cxcr5 = sum(cd20cxcr5, na.rm = TRUE), cd79bCoexp = sum(cd79bCoexp, na.rm=TRUE),
                                                                         Adjacent = mean(Adjacent, na.rm = T), Distal = mean(Distal, na.rm = T),
                                                                         Adjacent_cd8 = mean(Adjacent_cd8), 
                                                                         Distal_cd8 = mean(Distal_cd8),
                                                                         Adjacent_foxp3 = mean(Adjacent_foxp3),
                                                                         Distal_foxp3 = mean(Distal_foxp3), 
                                                                         sum_cells = sum(sum_cells, na.rm = TRUE))


df_cl_l$bcells <- rowSums(df_cl_l[,c('cd20', 'cd20cxcr5', 'cd79bCoexp')]) / df_cl_l$sum_cells

df_cl_l[,c('cd8', 'cd4', 'foxp3', 'cd20', 'cd20cxcr5', 'cd79bCoexp')] <- df_cl_l[,c('cd8', 'cd4', 'foxp3', 'cd20', 'cd20cxcr5', 'cd79bCoexp')] / df_cl_l$sum_cells


df_cl_l$cd8_foxp3 <- df_cl_l$cd8 / df_cl_l$foxp3
df_cl_l$cd4_foxp3 <- df_cl_l$cd4 / df_cl_l$foxp3

#df_cl_l[is.na(df_cl_l)] <- 0

#cutoff: c/cl median
#cutoff: all hs median
df_cl_l_1 <- subset(df_cl_l, Hotspots == 'Intratumoral IH')
df_cl_l_2 <- subset(df_cl_l, Hotspots == 'Peritumoral IH')

df_cl_l_1$adjacent_class <- ifelse(df_cl_l_1$Adjacent >= median(df_cl_l_1$Adjacent, na.rm = T), "High", "Low")
df_cl_l_2$adjacent_class <- ifelse(df_cl_l_2$Adjacent >= median(df_cl_l_2$Adjacent, na.rm = T), "High", "Low")

df_cl_l_1$bcell_class <- ifelse(df_cl_l_1$bcells >= median(df_cl_l_1$bcells, na.rm = T), "High", "Low")
df_cl_l_2$bcell_class <- ifelse(df_cl_l_2$bcells >= median(df_cl_l_2$bcells, na.rm = T), "High", "Low")

df_cl_l_3 <- rbind(df_cl_l_1, df_cl_l_2)

df_cl_l_3$adjacent_class <- factor(df_cl_l_3$adjacent_class, levels = c("High", "Low"))
df_cl_l_3$bcell_class <- factor(df_cl_l_3$bcell_class, levels = c("High", "Low"))

#intra tumor cd8_foxp3 high low
cell_class <- c("cd8", "cd20", "cd20cxcr5", "cd79bCoexp", "bcells", "cd4", "foxp3")
cell_label <- c("% of CD8+ T cells", "% of CD20+CXCR5- B cells", "% of CD20+CXCR5+ B cells", "% of CD79b+ B cells", "% of B cells", "% of CD4+FOXP3- cells", "% of CD4+FOXP3+ cells")

for (i in 1:length(cell_class)){
  print(cell_class[i])
  pdf(file.path("/Users/hzhang/Documents/report/lusc_b/plots_v4/fig4_cl_l_by_clump/intra_cd8_foxp3_low_high_cl_l_median/", paste0(cell_class[i], ".pdf")))
  
  p <- ggboxplot(subset(df_cl_l_3, !is.na(adjacent_class)), x = 'adjacent_class', y = cell_class[i], add = 'jitter', color = 'adjacent_class', palette = c("firebrick4", "tan3"), facet.by = "Hotspots") + 
    stat_compare_means(label = "p.format", label.x = 1.3, paired = F, label.y = 1) + 
    xlab('') + ylab(cell_label[i]) +
    cus_theme
  
  print(p)
  dev.off()
  
}

#bcell perc high low
cell_class <- c("Adjacent")
cell_label <- c("CD8+ / CD4+FOXP3+ within tumor nests")

for (i in 1:length(cell_class)){
  print(cell_class[i])
  pdf(file.path("/Users/hzhang/Documents/report/lusc_b/plots_v4/fig4_cl_l_by_clump/bcell_low_high_cl_l_median/", paste0(cell_class[i], ".pdf")))
  
  p <- ggboxplot(subset(df_cl_l_3, !is.na(bcell_class)), x = 'bcell_class', y = cell_class[i], add = 'jitter', color = 'bcell_class', palette = c("firebrick4", "tan3"), facet.by = "Hotspots") + 
    stat_compare_means(label = "p.format", label.x = 1.3, paired = F) + 
    xlab('') + ylab(cell_label[i]) +
    cus_theme
  
  print(p)
  dev.off()
  
}