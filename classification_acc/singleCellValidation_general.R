#setwd('X:/Dropbox (ICR)/yuanlab/Projects/lung/tracerx')
#setwd('/Users/kjabbar/Dropbox (ICR)/yuanlab/Projects/lung/tracerx')
setwd('/Users/hzhang/Documents/project/sum/T_accuracy')
source('Kfun.R')
library(dplyr)
library(caret)

  #10pixels = 4.54um
#15pixels = 6.84um
mydist <- function(row){
  dists <- (row[["x"]] - df2$x)^2 + (row[["y"]]- df2$y)^2
  return(cbind(df2[which.min(dists),], distance = min(dists)))
}

#R:\tracerx\tracerx100\HE_region_specific\results\classification\20180124\csv; regional DATA
#combine
gt = '/Volumes/proj4/HZ/backup/desktop/project/lung_reg/code/polyscope_annotation/cell_labels_tile/Tcell/celllabels' #da.csv, name have to be the same as dl
dl = '/Volumes/proj4/HZ/backup/desktop/project/lung_reg/code/polyscope_annotation/cell_labels_tile/Tcell/celllabels_dl' #labels/slidenames/da.csv
dl_corrected = '/Volumes/proj4/HZ/backup/desktop/project/lung_reg/code/polyscope_annotation/cell_labels_tile/Tcell/celllabels_dl_corrected/cv1'
files = dir(path=gt, pattern = '*.ndpi')

csvs <- vector("list", length(files))
names(csvs) <- files

csv_name <- list.files(paste0(dl, '/',files[1]), pattern = '*.csv')
csv_dl <- read.csv(paste0(dl, '/',files[1], '/', csv_name[1]))
csv_corrected <- read.csv(paste0(dl_corrected, '/',files[1], '/', csv_name[1]))
csv_merge <- merge(csv_dl, csv_corrected, by = c("V2","V3"), all.x = TRUE)
levels(csv_merge$V1.y) <- c(levels(csv_merge$V1.y), "None")
csv_merge[is.na(csv_merge)] <- "None"
csv_merge$V1 <- ifelse(csv_merge$V1.x=="uc" | csv_merge$V1.y=="uc", "uc",
                       ifelse(csv_merge$V1.x == "cd8", "cd8", 
                              ifelse(csv_merge$V1.y == "foxp3", "foxp3",
                                     ifelse(csv_merge$V1.y=="cd4", "cd4", "None"))))
csv_merge$V1 <- as.factor(csv_merge$V1)


for (i in 1:length(files)){
  
  filePathGT = file.path(gt, files[i])
  filePathDL = file.path(dl, files[i])
  print(files[i])
  
  Das = dir(path=filePathGT, pattern = '*.csv')
  
  csvs[[i]] <- vector("list", length(Das))
  names(csvs[[i]]) <- Das
  
  for (t in 1:length(Das)){
    
    DaGT <- read.csv(file.path(filePathGT, Das[t]))
    DaDL <- read.csv(file.path(filePathDL, Das[t]))
    print(Das[t])
    
    DaGT = DaGT[, c(2:3, 1)]
    DaDL = DaDL[, c(2:3, 1)]
    names(DaGT) <- c("x", "y", "class")
    names(DaDL) <- c("x", "y", "class")
    
    #remove the new cyan/penumocytes in overall 
    #DaGT = DaGT[! DaGT$class=="p",]
    #DaGT = DaGT[! DaGT$class=="o",]
    #DaGT = DaGT[! DaGT$class=="f",]
    
    GTcells = nrow(DaGT)
    
    if (GTcells>0){
      df2 = DaDL
      df1 = DaGT
      DaR <- cbind(df1, do.call(rbind, lapply(1:nrow(df1), function(x) mydist(df1[x,]))))
      
      names(DaR) <- c("xgt", "ygt", "classgt", "xdl", "ydl", "classdl", "distance")
      DaR$detect_gt <- "True"
      DaR$detect_dl <-ifelse(DaR$distance<=25, "True", "False") ##CHECK!
      DaR$sample = files[i]
      csvs[[i]][[t]] <- DaR
    }
    
  }
}
#comb = flatten(csvs)
#comb2 <- bind_rows(comb, .id = "file_name")

#comb2 <- bind_rows(csvs[[1]][[1]], csvs[[2]][[1]], csvs[[2]][[2]])
comb2 <- bind_rows(csvs[[1]][[1]], csvs[[1]][[2]], csvs[[2]][[1]], csvs[[2]][[2]])

comb2$classgt <- as.factor(as.character(comb2$classgt))
comb2$classdl <- as.factor(as.character(comb2$classdl))
comb2$detect_gt <- factor(comb2$detect_gt, levels = c("False", "True"))
comb2$detect_dl <- factor(comb2$detect_dl, levels = c("False", "True"))

#levels(comb2$classgt) <- c("l", "t", "o", "f")
#levels(comb2$classdl) <- c("f", "l", "o", "t")

detach("package:plyr", unload=TRUE)
detach("package:dplyr", unload=TRUE)
library(plyr)
library(dplyr)
    
comb2 %>% group_by(xdl, ydl) %>% summarise(n = n())
comb2 %>% group_by(xgt, ygt) %>% summarise(n = n())

#comb2_positive <- subset(comb2, !(classgt %in% c("hem", "cxcr5")))
levels(comb2$classgt)[4] <- 'uc'
comb2_positive <- subset(comb2, !(classgt %in% c("uc")))

comb2_detected_only <- subset(comb2, detect_dl == "True")
cfm <- confusionMatrix(comb2_detected_only$classdl, comb2_detected_only$classgt)
cfm_detect <- confusionMatrix(comb2$detect_dl, comb2$detect_gt)
cfm_detect_pos <- confusionMatrix(comb2_positive$detect_dl, comb2_positive$detect_gt)

#saveRDS(cfm, '/home/adminhzhang/Documents/project/lung_reg/sum/B_class_accuracy/cfm_20pixel.rds')
saveRDS(cfm, '/home/adminhzhang/Documents/project/lung_reg/sum/cfm_Tcell.rds')
saveRDS(cfm_detect, '/home/adminhzhang/Documents/project/lung_reg/sum/cfm_Tcell_detect.rds')
