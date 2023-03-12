##############################################################################=
# Count T, B cells in regions and IH
##############################################################################=
library(plyr)
library(dplyr)
library(reshape2)
library(tidyverse)

# detach("package:plyr", unload=TRUE)
# detach("package:dplyr", unload=TRUE)

#Count B, T cells in IH====
hs_match <- function(lym_data, data_wg){
  levels(lym_data$region)[which(levels(lym_data$region) == "")] <- "None"
  levels(lym_data$region_idx)[which(levels(lym_data$region_idx) == "")] <- "None"
  levels(lym_data$join_idx)[which(levels(lym_data$join_idx) == "")] <- "None"
  
  data_wg$X <- NULL
  
  lym_data$region_idx <- gsub("_Da.*", "", lym_data$join_idx)
  
  lym_data$x.s <- ceiling(lym_data$x / s) * s
  lym_data$y.s <- ceiling(lym_data$y / s) * s
  
  lym_hs <- merge(lym_data, data_wg, by.x = c('x.s', 'y.s'), by.y = c('x', 'y'), all.x = TRUE, all.y=FALSE,sort = TRUE)
  
  cell_count <- aggregate(x ~ class + x.s + y.s + region + region_idx + join_idx, data = lym_hs, length)
  colnames(cell_count)[which(colnames(cell_count) == 'x')] <- 'cell.count'
  
  cell_count <- unique(merge(cell_count,lym_hs[, !colnames(lym_hs) %in% c('x','y')], by = c('x.s', 'y.s','class', 'region', 'region_idx', 'join_idx'),sort = TRUE))

  cell_count$hs <- ifelse(!is.na(cell_count$G.l) & !is.na(cell_count$G.c) & cell_count$G.l >= 3.886 & cell_count$G.c >= 3.886, 'cl',
                          ifelse(!is.na(cell_count$G.c) & cell_count$G.c >= 3.886, 'c',
                                 ifelse(!is.na(cell_count$G.l) & cell_count$G.l >= 3.886, 'l', 'else')))
  return(cell_count)
}

Outputdir <- '/Users/hzhang/Documents/project/sum/TLS_annotation/hs_he_region/'
save_dir <- '/Users/hzhang/Documents/project/sum/final/TLS_annotation/cellcounts_final/'
t_folder <- "/Users/hzhang/Documents/project/sum/final/TLS_annotation/CellPos/T/"
b_folder <- '/Users/hzhang/Documents/project/sum/final/TLS_annotation/CellPos/B/'
he_folder <- '/Users/hzhang/Documents/project/sum/final/TLS_annotation/CellPos/HE/'
ffs_b <- dir(b_folder, pattern = '.csv$', full.names = TRUE)
ffs_t <-  dir(t_folder, pattern = '.csv$', full.names = TRUE)
ffs_he <-  dir(he_folder, pattern = '.csv$', full.names = TRUE)
s = 10

ffs <- gsub("[^[:digit:]]", "",  gsub(".*/", "",  ffs_he))

dir.create(save_dir)


for (i in 1:length(ffs)){
  print(ffs[i])
  hs_df <- read.csv(paste0(Outputdir, ffs[i], '_cl.csv'))

  df_t <- read.csv(ffs_t[i])

  if (nrow(df_t)>0) {

  df_t <- hs_match(df_t, hs_df)
  df_t$slide <- ffs[i]

  write.csv(df_t, paste0(save_dir, ffs[i], '_hscellcount_t.csv'))}

  df_b <- read.csv(ffs_b[i])

  if (nrow(df_b)>0) {
  df_b <- hs_match(df_b, hs_df)
  df_b$slide <- ffs[i]

  write.csv(df_b, paste0(save_dir, ffs[i], '_hscellcount_b.csv'))}
  
  df_he <- read.csv(ffs_he[i])

  if (nrow(df_he)>0) {
  df_he <- hs_match(df_he, hs_df)
  df_he$slide <- ffs[i]

  write.csv(df_he, paste0(save_dir, ffs[i], '_hscellcount_he.csv'))}
  
}

#######################1
# Join by regions====

sort_df <- function(df_t, df_b){
  df_bt <- rbind(df_t, df_b)
  df_bt$X.1<-NULL
  df_bt$X <- NULL
  df_bt$cell.count.c<- NULL
  df_bt$cell.count.l<- NULL
  
  td <- df_bt %>% group_by(x.s,y.s,hs, region, slide, region_idx, join_idx) %>%
    spread(class, sum(cell.count))
  hs_bt_all <- td

  return(hs_bt_all)
}

get_cell_count_all <- function(df){
  df$X.1 <- NULL
  df$X <- NULL

  # df_1 <- df %>% group_by(x.s, y.s, slide, region, region_idx, join_idx) %>% summarise(cell.count.all = sum(cell.count))
  # df_2 <- merge(df_1, df[c('x.s','y.s','hs','slide','region')], by = c('x.s','y.s','slide'), all.x = T)
  # df_2 <- unique(df_2)

  df_2 <- df %>% group_by(x.s, y.s, slide, region, region_idx, join_idx) %>% mutate(cell.count.all = sum(cell.count))
  
  
  df_l <- subset(df, class =='l')
  df_l$cell.count.l = df_l$cell.count
  df_all <- merge(df_2, df_l[c('x.s','y.s','hs','slide','region','region_idx','join_idx', 'cell.count.l')], all.x=T, all.y=T)
  
  df_c <- subset(df, class =='t')
  df_c$cell.count.c = df_c$cell.count
  df_all <- merge(df_all, df_c[c('x.s','y.s','hs','slide','region','region_idx','join_idx', 'cell.count.c')], all.x=T, all.y=T)
  
  df_all[,c("class", "cell.count")] <- NULL
  df_all <- unique(df_all)
  df_all <- df_all %>% group_by(x.s,y.s,hs,slide,region,region_idx,join_idx) %>% summarise(cell.count.l= sum(cell.count.l, na.rm = T),
                                                                                                         cell.count.c = sum(cell.count.c, na.rm = T),
                                                                                                         cell.count.all = sum(cell.count.all, na.rm = T))
  # 
  # df_all[is.na(df_all$cell.count.l), 'cell.count.l'] = 0
  # df_all[is.na(df_all$cell.count.c), 'cell.count.c'] = 0
  
  return(df_all)
}


generate_sum_all <- function(data_dir, ffs){
  sum_all <- data.frame()
  for (f in ffs){
    df_tf <- file.path(data_dir, paste0(f, '_hscellcount_t.csv'))
    df_bf <- file.path(data_dir, paste0(f, '_hscellcount_b.csv'))
    df_hef <- file.path(data_dir, paste0(f, '_hscellcount_he.csv'))
    
    if(file.exists(df_hef)){
      print(f)
      df_t <- read.csv(df_tf)
      df_b <- read.csv(df_bf)
      df_he <- read.csv(df_hef)
      
      hs_bt_all <- sort_df(df_t, df_b)
      
      he_all <- get_cell_count_all(df_he)
      
      sum_all_slide <- merge(hs_bt_all, he_all, all.x = T, all.y = T)
      
      #uniform dataframe
      slide_contain <- data.frame(matrix(ncol = 22, nrow = nrow(sum_all_slide)))
      colnames(slide_contain) <- c("x.s","y.s","hs","region","slide","region_idx","join_idx", "G.c","G.l", 
                                   "cd4","cd8" ,"foxp3","uc","cd20","cd20cxcr5","cd79bCoexp",
                                   "cxcr5","hem","p40","cell.count.all","cell.count.l","cell.count.c")
      for (i in colnames(sum_all_slide)){
        slide_contain[,i] = sum_all_slide[,i]
      }
      
      sum_all <- rbind(sum_all, slide_contain)
    }
  }
  
  for (i in c(10:22)){
    sum_all[is.na(sum_all[,i]),i] <- 0 
  }
  dir.create(file.path(data_dir, "sum"))
  write.csv(sum_all, file = file.path(data_dir,"sum", "sum_all_add_tls_new.csv"))
  return(sum_all)
}


data_dir <- '/Users/hzhang/Documents/project/sum/final/TLS_annotation/cellcounts_final/'
sum_all <- generate_sum_all(data_dir, ffs)

#Sort dataset====
# sum_all <- read.csv("/Users/hzhang/Documents/project/sum/final/TLS_annotation/cellcounts/sum/sum_all_add_tls.csv")
sum_all <- read.csv("/Users/hzhang/Documents/project/sum/final/TLS_annotation/cellcounts_final/sum/sum_all_add_tls_new.csv")

subset(sum_all, region == "None")
head(sum_all)
df <- sum_all %>% group_by(join_idx, slide, region) %>% summarise(n())
region_count <- data.frame(table(df$slide, df$region))

#change region of sum_all without join_idx to else
sum_all[sum_all[,"join_idx"]=="" | is.na(sum_all[,"join_idx"]), "region"] <- "None"

sum_all_no_join <- sum_all
#use tls/lag regardless of whether there's join_idx or not
sum_all_no_join <- sum_all_no_join %>% group_by(x.s, y.s, hs, slide, region) %>% summarise(cd4 = mean(cd4),
                                                                                           cd8 = mean(cd8),
                                                                                           foxp3 = mean(foxp3),
                                                                                           uc = mean(uc),
                                                                                           cd20 = mean(cd20),
                                                                                           cd20cxcr5 = mean(cd20cxcr5),
                                                                                           cd79bCoexp = mean(cd79bCoexp),
                                                                                           cxcr5 = mean(cxcr5),
                                                                                           hem = mean(hem),
                                                                                           p40  = mean(p40),
                                                                                           cell.count.all = mean(cell.count.all),
                                                                                           cell.count.l= mean(cell.count.l),
                                                                                           cell.count.c = mean(cell.count.c))

colnames(sum_all_no_join)[which(colnames(sum_all_no_join) == 'region')] = 'Region'
colnames(sum_all_no_join)[which(colnames(sum_all_no_join) == 'hs')] = 'Hotspots'
levels(sum_all_no_join$Region) <- c('LAG', 'UD', 'TLS')
sum_all_no_join$Region <- factor(sum_all_no_join$Region, levels=c('TLS','LAG','UD'))
levels(sum_all_no_join$Hotspots) <- c('Cancer', 'Cancer-immune', 'Else', 'Immune')

# write.csv(sum_all_no_join, "/Users/hzhang/Documents/project/sum/final/TLS_annotation/cellcounts/sum/sum_all_no_join_1.csv")
