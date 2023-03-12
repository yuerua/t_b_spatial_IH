##############################################################################=
# Identify disconnected IH using clump in raster
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


#change region labels====
Outputdir <- '/Users/hzhang/Documents/project/sum/TLS_annotation/hs_he_region/'
save_dir <- '/Users/hzhang/Documents/project/sum/final/TLS_annotation/CellPos/CellPos_region_and_hs/'
t_folder <- "/Users/hzhang/Documents/project/sum/final/TLS_annotation/CellPos/T/"
b_folder <- '/Users/hzhang/Documents/project/sum/final/TLS_annotation/CellPos/B/'
he_folder <- '/Users/hzhang/Documents/project/sum/final/TLS_annotation/CellPos/HE/'
ffs_b <- dir(b_folder, pattern = '.csv$', full.names = TRUE)
ffs_t <-  dir(t_folder, pattern = '.csv$', full.names = TRUE)
ffs_he <-  dir(he_folder, pattern = '.csv$', full.names = TRUE)
s = 10
dir.create(save_dir)

ffs <- gsub("[^[:digit:]]", "",  gsub(".*/", "",  ffs_he))

hs_match <- function(lym_data, data_wg){
  lym_data$x.s <- ceiling(lym_data$x / s) * s
  lym_data$y.s <- ceiling(lym_data$y / s) * s
  
  lym_hs <- merge(lym_data, data_wg, by.x = c('x.s', 'y.s'), by.y = c('x', 'y'), all.x = TRUE, all.y=FALSE,sort = TRUE)
  
  
  lym_hs$hs <- ifelse(!is.na(lym_hs$G.l) & !is.na(lym_hs$G.c) & lym_hs$G.l >= 3.886 & lym_hs$G.c >= 3.886, 'cl',
                      ifelse(!is.na(lym_hs$G.c) & lym_hs$G.c >= 3.886, 'c',
                             ifelse(!is.na(lym_hs$G.l) & lym_hs$G.l >= 3.886, 'l', 'else')))
  return(lym_hs)
}




for (i in 1:length(ffs)){
  print(ffs[i])
  df_slide <- data.frame()
  hs_df <- read.csv(paste0(Outputdir, ffs[i], '_cl.csv'))
  
  df_t <- read.csv(ffs_t[i])
  levels(df_t$region)[which(levels(df_t$region) == "")] <- "UD"
  # df_t <- subset(df_t, region == 'None')
  if (nrow(df_t)>0) {
    # df_t$region = NULL
    df_t <- hs_match(df_t, hs_df)
    df_t$slide <- ffs[i]
    # df_t$region = 'None'
    df_slide <- rbind(df_slide, df_t)}
  
  df_b <- read.csv(ffs_b[i])
  levels(df_b$region)[which(levels(df_b$region) == "")] <- "UD"
  # df_b <- subset(df_b, region == 'None')
  if (nrow(df_b)>0) {
    df_b <- hs_match(df_b, hs_df)
    df_b$slide <- ffs[i]
    # df_b$region = 'None'
    df_slide <- rbind(df_slide, df_b)}
  
  df_he <- read.csv(ffs_he[i])
  levels(df_he$region)[which(levels(df_he$region) == "")] <- "UD"
  # df_he <- subset(df_he, region == 'None')
  if (nrow(df_he)>0) {
    df_he <- hs_match(df_he, hs_df)
    df_he$slide <- ffs[i]
    # df_he$region = 'None'
    df_slide <- rbind(df_slide, df_he)}
  
  # df_all <- rbind(df_t, df_b)
  # hs_bt_all <- rbind(hs_bt_all, df_all)
  write.csv(df_slide, file.path(save_dir, paste0(ffs[i], '.csv')))
}

#=========#
#Identify individual IH====
setwd("/Users/hzhang/Documents/project/sum/final/src")
maskDir <- "/Users/hzhang/Documents/project/sum/CellPos_all/HE_r/mask_r/"
cellPosDir <- "/Users/hzhang/Documents/project/sum/final/TLS_annotation/CellPos/CellPos_region_and_hs/"

Outputdir <- "/Users/hzhang/Documents/project/sum/final/clump_l_cl"

res <- 228*2*16/1000
dist_th <- 250

dir.create(Outputdir, showWarnings = F)

cell_type_list <- c('cd8', 'foxp3')


cellPos_ext <- "_Cellpos_add_hs.csv"

cus_theme <-theme(plot.title=element_text(size=16),
                  axis.text=element_text(size=14),
                  axis.title=element_text(size=16),
                  legend.title = element_text(size = 14),
                  legend.text = element_text(size = 14))


#========#
get_near_dist <- function(df_a, df_b){
  #df_a: hs dataframe #df_b: cell pos
  #change k to number of neighbours 
  closest <- RANN::nn2(data = df_a[,c('x.s','y.s')], query = df_b[,c('x','y')], k = 1)
  
  dist_df <- df_b
  dist_df$c_x <- df_a[closest$nn.idx,]$x.s
  dist_df$c_y <- df_a[closest$nn.idx,]$y.s
  dist_df$clump_id <- df_a[closest$nn.idx,]$clump_id
  dist_df$dist <- closest$nn.dists[,1]
  
  return(dist_df)
}


ffs <- dir(cellPosDir, pattern = '.csv$')
files = c()
for (i in 1:length(ffs)) {
  files.vec <- c('/',cellPosDir,'/', ffs[i])
  files[i] <- paste(files.vec, collapse = "")
}
ffs_ori <- ffs


for (i in 1:length(ffs_ori)) {
  # i = 3
  print(paste0('Processing ', i, ' of ', length(ffs_ori)))
  rdata.file = read.csv(files[i])
  rdata.file$X <- NULL
  # CellPos <- rdata.file[, c('class', 'x', 'y')]
  
  df_hs <- subset(rdata.file, hs %in% c('cl', 'l'))
  
  file_name <- gsub(".csv", "", ffs_ori[i])
  

  Mask <- readImage(paste0(maskDir,'/', file_name, '_HE.ndpi_Mask.jpg'))
  s <- 10
  wp <- nrow(Mask)
  hp <- ncol(Mask)
  
  df_hs_sum <- df_hs %>% group_by(x.s, y.s, hs) %>% summarise(sum_n = n())
  
  l.zm.c <- matrix(nrow=wp/s, ncol=hp/s)
  
  for (i in 1:nrow(df_hs_sum)) {
    l.zm.c[df_hs_sum$x.s[i]/s,df_hs_sum$y.s[i]/s] = 1
  }
  # l.zm.c <- l.zm.c[,ncol(l.zm.c):1]
  l.zm.c <- raster(l.zm.c)
  
  # wid1 <- wid/scalefactor
  # hei1 <- hei/scalefactor
  
  # x.i <- 1*(1:(nrow(l.zm.c)))
  # y.i <- 1*(1:(ncol(l.zm.c)))
  
  rc <- clump(l.zm.c, directions = 8)
  
  clump_id <- getValues(rc)
  l.zm.c[] <- clump_id
  
  clump_df <- setNames(melt(as.matrix(l.zm.c)), c('x.s','y.s','clump_id'))
  clump_df$x.s <- clump_df$x.s * s
  clump_df$y.s <- clump_df$y.s * s
  df_hs_sum <- merge(df_hs_sum, clump_df, by = c('x.s', 'y.s'), all.x = T)
  
  # plot
  # clump_id
  # xy <- xyFromCell(rc,1:ncell(rc))
  # df <- data.frame(xy, clump_id, is_clump = rc[] %in% freq(rc, useNA = 'no')[,1])
  # df[df$is_clump == T, ]
  # 
  # dfm <- ddply(df[df$is_clump == T, ], .(clump_id), summarise, xm = mean(x), ym = mean(y))
  # 
  # 
  # par(bg=NA)
  # plot (rc)
  # # text(dfm[, 2:3], labels = dfm$clump_id, cex.lab=1,  cex.main=1)
  # 
  # dev.copy(png,'/Users/hzhang/Documents/project/sum/final/plot/clump_l_cl/51464_clump_text.png')
  # dev.off()
  
  write.csv(df_hs_sum, file.path(Outputdir, paste0(file_name, '_hs_sum_clump_df.csv')))
}
# dist_df

#===========================================================#
#sum all slides====
dist_all_sum_all_slides <- data.frame()


for (i in 1:length(ffs_ori)) {
  print(paste0('Processing ', i, ' of ', length(ffs_ori)))
  file_name <- gsub(".csv", "", ffs_ori[i])
  dist_all_sum <- read.csv(file.path(Outputdir, paste0(file_name, '_hs_sum_clump_df.csv')))
  dist_all_sum$X <- NULL
  dist_all_sum$slide <- file_name
  dist_all_sum <- dist_all_sum %>% group_by(hs, clump_id) %>% mutate(clump_area = n())
  
  if (i==1){
    dist_all_sum_all_slides <- dist_all_sum
  }else{
    dist_all_sum_all_slides <- rbind(dist_all_sum_all_slides, dist_all_sum)
  }
  
}

dir.create(file.path(Outputdir, "sum"))
write.csv(dist_all_sum_all_slides, file.path(Outputdir,"sum", "dist_all_sum_all_slides.csv"))
#======================#

#correlation with b cell percentages
# load('/Users/hzhang/Documents/report/lusc_b/scripts/data_2.RData')

#merge clump_id with cell no.====
dist_all_sum_all_slides <- read.csv("/Users/hzhang/Documents/project/sum/final/src/prepare_final_script/data/validation/new/dist_all_sum_all_slides.csv")
dist_all_sum_all_slides$X <- NULL

sum_all_no_join <- read.csv("/Users/hzhang/Documents/project/sum/final/src/prepare_final_script/data/validation/new/sum_all_add_tls_new.csv")
sum_all_no_join$X <- NULL

colnames(sum_all_no_join)[which(colnames(sum_all_no_join)=='hs')] <- 'Hotspots'
levels(sum_all_no_join$Hotspots) <- c('Cancer', 'Cancer-immune', 'Else', 'Immune')

colnames(sum_all_no_join)[which(colnames(sum_all_no_join) == 'region')] = 'Region'
levels(sum_all_no_join$Region) <- c('LAG', 'UD', 'TLS')
sum_all_no_join$Region <- factor(sum_all_no_join$Region, levels=c('TLS','LAG','UD'))

levels(dist_all_sum_all_slides$hs) <- c('Cancer-immune', 'Immune')
colnames(dist_all_sum_all_slides)[which(colnames(dist_all_sum_all_slides)=='hs')] <- 'Hotspots'
colnames(dist_all_sum_all_slides)[which(colnames(dist_all_sum_all_slides)=='sum_n')] <- 'clump_area'


sum_all_no_join_clump <- merge(sum_all_no_join, dist_all_sum_all_slides, by = c('x.s', 'y.s', 'slide', 'Hotspots'), all.x = T)


write.csv(sum_all_no_join_clump, "/Users/hzhang/Documents/project/sum/final/src/prepare_final_script/data/validation/new/sum_all_no_join_clump_new.csv")


# sum_all_no_join_clump <- merge(sum_all_no_join_clump, df_cell_cd8_foxp3[, c('clump_id', 'slide', 'Adjacent', 'Distal')], all.x = T)