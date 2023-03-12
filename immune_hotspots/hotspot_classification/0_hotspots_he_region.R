##############################################################################=
#Categorise 50x50 um^2 grids into intratumoral IH (Cancer-immune, ci), peritumoral IH (Immune, i), 
#Cancer hotspots (Cancer, c), and non-hotspot area (Else)
##############################################################################=
library(plyr)
library(dplyr)
library(ggpubr)
library(EBImage)
library(spdep)
library(reshape2)
library(tidyverse)

#hotspot categorisation====
setwd('/Users/hzhang/Documents/project/sum/TLS_annotation/code')

maskDir <- "/Users/hzhang/Documents/project/sum/CellPos_all/HE_r/mask_r/"
he_folder <- '/Users/hzhang/Documents/project/sum/TLS_annotation/CellPos/HE_region/'
Outputdir <- '/Users/hzhang/Documents/project/sum/TLS_annotation/hs_he_region/'
#hs_match_dir <- '/home/adminhzhang/Documents/project/lung_reg/sum/hotspot_bt_DRDIN_foxp3correct/hotspot_count/'
dir.create(Outputdir)

ffs_he <- dir(he_folder, pattern = '.csv$', full.names = TRUE)

ffs <- gsub("[^[:digit:]]", "",  gsub(".*/", "",  ffs_he))

s <- 10
n <- 4
min.d <- 0
max.d <- (n * s * sqrt(2))

for (f in 1:length(ffs_he)){
  k = f
  print(paste0('Processing ', f, ' of ', length(ffs_he)))
  
  rdata.file_he <- read.csv(ffs_he[f])
  #rdata.file_he$subclass <- rdata.file_he$class
  #rdata.file_he <- subset(rdata.file_he, subclass =='l')
  #rdata.file_he$class <- 'l'
  
  # rdata.file_b <- read.csv(ffs_b[f])
  # rdata.file_b$subclass <- rdata.file_b$class
  # 
  # rdata.file_c <- subset(rdata.file_b, subclass == 'p40')
  # rdata.file_c$class <- 'c'
  
  CellPos.xy <- rdata.file_he
  
  slide_name <- ffs[f]
  
  Mask <- readImage(paste0(maskDir, ffs[f], '_HE.ndpi_Mask.jpg'))
  w <- nrow(Mask)
  h <- ncol(Mask)
  
  w.s <- w / s
  h.s <- h / s
  
  c.m <- matrix(nrow = w, ncol = h) #b=cancer, t=lymphocyte
  l.m <- matrix(nrow = w, ncol = h)
  cellpos.c <- subset(CellPos.xy, CellPos.xy$class == 't')
  cellpos.l <- subset(CellPos.xy, CellPos.xy$class == 'l')
  
  #?cell x larger than w?
  cellpos.c <- subset(cellpos.c, x <= w&y<=h)
  cellpos.l <- subset(cellpos.l, x <= w&y<=h)
  
  
  for (i in 1:nrow(cellpos.c)) {
    c.m[cellpos.c$x[i], cellpos.c$y[i]] <- 1
  }
  
  nadata <- which(is.na(c.m) == T)
  c.m[nadata] <- 0
  
  
  for (i in 1:nrow(cellpos.l)) {
    l.m[cellpos.l$x[i], cellpos.l$y[i]] <- 1
  }
  
  nadata <- which(is.na(l.m) == T)
  l.m[nadata] <- 0
  
  c.m.s <- matrix(nrow = w.s, ncol = h.s)
  
  for (i in 1:w.s) {
    for (j in 1:h.s) {
      c.m.s[i, j] <-
        sum(c.m[((i * s) - (s - 1)):(i * s), ((j * s) - (s - 1)):(j * s)])
    }
  }
  
  
  l.m.s <- matrix(nrow = w.s, ncol = h.s)
  
  for (i in 1:w.s) {
    for (j in 1:h.s) {
      l.m.s[i, j] <-
        sum(l.m[((i * s) - (s - 1)):(i * s), ((j * s) - (s - 1)):(j * s)])
    }
  }
  
  #Warning!!! Next step takes a long time!(edited)
  Mask_df <- setNames(melt(Mask), c('x','y','value'))
  Mask_df$row <- ceiling(Mask_df$x / s)
  Mask_df$col <- ceiling(Mask_df$y / s)
  data <- Mask_df %>% group_by(row, col) %>% summarise(value = sum(value)) %>% as.matrix()
  data <- data[data[,3]>= ((s * s) / 2),] %>% as.data.frame()
  data <- subset(data, row <= floor(w.s) & col <= floor(h.s))
  data$value <- NULL
  data <- data[order(data$col),]
  # mask.s <- matrix(nrow = w.s, ncol = h.s)
  # for (i in 1:w.s) {
  #   for (j in 1:h.s) {
  #     mask.s[i, j] <- sum(Mask[((i * s) - (s - 1)):(i * s), ((j * s) - (s - 1)):(j *
  #                                                                                  s)])
  #   }
  # }
  # 
  # #Why? Grayscale each pixel 0:1, threshold = s*s/2, data is the index of sum grids
  # data <- which(mask.s >= ((s * s) / 2), arr.ind = T) 
  
  
  c.counts <- c()
  for (i in 1:nrow(data)) {
    c.counts[i] <- c.m.s[data$row[i], data$col[i]]
  }
  
  l.counts <- c()
  for (i in 1:nrow(data)) {
    l.counts[i] <- l.m.s[data$row[i], data$col[i]]
  }
  
  data <- data * s
  
  data <- cbind(data, c.counts, l.counts)
  colnames(data) <- c('x', 'y','cell.count.c', 'cell.count.l')
  
  #Remove all position without c cells
  # a <- which(data$cell.count.c == 0 )
  # if (length(a) != 0) {
  #   data <- data[-a, ]
  # } else{
  # }
  
  sq.count <- nrow(data)
  
  if (sq.count >= 30) {
    data.copy <- data
    coordinates(data.copy) <- c('x', 'y')
    nlist <- dnearneigh(data.copy, d1 = min.d, d2 = max.d)
    nbs <-
      nb2listw(include.self(nlist),
               style = 'B',
               zero.policy = TRUE)
    
    G.c <- localG(data.copy$cell.count.c, nbs)
    G.l <- localG(data.copy$cell.count.l, nbs)
    
    class(G.c) = 'numeric'
    class(G.l) = 'numeric'
    
    data.wg = cbind(data, G.c, G.l)
    
    if (sq.count >= 1 & sq.count < 50) {
      g <- 1.645
    }
    if (sq.count >= 50 & sq.count < 100) {
      g <- 3.083
    }
    if (sq.count >= 100 & sq.count < 1000) {
      g <- 3.289
    }
    if (sq.count >= 1000) {
      g <- 3.886
    }
    
    hotcells.c <- subset(data.wg, data.wg$G.c >= g)
    hotcells.l <- subset(data.wg, data.wg$G.l >= g)
    
    c.hotspot.count <- nrow(hotcells.c)
    l.hotspot.count <- nrow(hotcells.l)
    
    store.rdata <- paste0(c(ffs[f]), '_cl.rdata')
    save(
      data.wg,
      h,
      w,
      s,
      min.d,
      max.d,
      g,
      hotcells.c,
      hotcells.l,
      sq.count,
      nlist,
      file = paste0(Outputdir, store.rdata)
    )
    print('Done')
    
    data.wg$hs <- ifelse(!is.na(data.wg$G.c) & !is.na(data.wg$G.l) & data.wg$G.c >= g & data.wg$G.l >=g, 'cl',
                         ifelse(!is.na(data.wg$G.c) & data.wg$G.c >= g, 'c', 
                                ifelse(!is.na(data.wg$G.l) & data.wg$G.l >= g, 'l', 'else')))
    data.wg$slide <- slide_name
    fn <- paste0(c(ffs[f]), '_cl.csv')
    write.csv(data.wg, paste0(Outputdir, fn))
  }
}

######################################1
#summary IH statistics====
# get all processed hotspot files (these were processed using getHotspotsLocal.R)
folder <- "/Users/hzhang/Documents/project/sum/TLS_annotation/hs_he_region/"
library(survival)

s <- 10
n <- 4
min.d <- 0
max.d <- (n * s * sqrt(2))
library(spdep)
ffs <- dir(folder, pattern = '.rdata$')
files = c()


for (i in 1:length(ffs)) {
  files.vec <- c(folder, ffs[i])
  files[i] <- paste(files.vec, collapse = "")
}

#pre variables

id <- c()
for (i in 1:length(files)){
  #fix this stupid filename
  id[i] <- ffs[i]
  id[i] <- gsub("[^[:digit:]]", "",  gsub(".*/", "",  ffs_he[i]))
}

hotspots<-data.frame(file_name=character(), squares=character(),cancer.hs=character(),lympho.hs=character(),no.of.clhs.overlap=character(),
                     cancer.fractional.overlap=character(),immune.fractional.overlap=character(),normalised.overlap=character(),
                     hCImuPer=character(), hCPer=character(), hImuPer=character(), nonhotPer=character())

overlap.data <- data.frame()
overlap.count <- c()
fractional.cancer.overlap <- c()
fractional.lympho.overlap <- c()
norm.overlap <- c()

for (f in 1:length(files)) {
  print(paste0('Processing ', f, ' of ', length(files)))
  rdata.file <-  files[f]
  load(rdata.file)
  
  
  #calculate the fractions using Sidra's getisord_ptc
  cancer.hs <- nrow(hotcells.c)
  lympho.hs <- nrow(hotcells.l)
  squares <- sq.count
  cl.hs <- subset(data.wg, data.wg$G.c >= 3.886 & data.wg$G.l >= 3.886)
  overlap.count[f] <- nrow(cl.hs)
  fractional.cancer.overlap[f] <- overlap.count[f]/nrow(hotcells.c)
  fractional.lympho.overlap[f] <- overlap.count[f]/nrow(hotcells.l)
  norm.overlap[f] <- overlap.count[f]/sq.count
  
  #new stuff
  # calculate per for each hotspot + the empty habitat (nonhot)
  nonhot.hs <- subset(data.wg, data.wg$G.c < 3.886 & data.wg$G.l < 3.886)
  
  # also add the fci by subtracting from both fi/fc 
  hCImuPer <- c()
  hCPer <- c()
  hImuPer <- c()
  nonhotPer <- c()
  
  hCImuPer[f] = overlap.count[f]/(nrow(hotcells.c)+nrow(hotcells.l)+nrow(nonhot.hs)+overlap.count[f])
  hCPer[f] = nrow(hotcells.c)/(nrow(hotcells.c)+nrow(hotcells.l)+nrow(nonhot.hs)+overlap.count[f])
  hImuPer[f] = nrow(hotcells.l)/(nrow(hotcells.c)+nrow(hotcells.l)+nrow(nonhot.hs)+overlap.count[f])
  nonhotPer[f] = nrow(nonhot.hs)/(nrow(hotcells.c)+nrow(hotcells.l)+nrow(nonhot.hs)+overlap.count[f])
  
  hotspots <- rbind (hotspots, data.frame(file_name=id[[f]],squares=squares,cancer.hs=cancer.hs,
                                          lympho.hs=lympho.hs,no.of.clhs.overlap=overlap.count[[f]],
                                          cancer.fractional.overlap=fractional.cancer.overlap[[f]],
                                          immune.fractional.overlap=fractional.lympho.overlap[[f]],
                                          normalised.overlap=norm.overlap[[f]], hCImuPer=hCImuPer[[f]], hCPer=hCPer[[f]],
                                          hImuPer=hImuPer[[f]], nonhotPer=nonhotPer[[f]]))
}


#rename
hotspots$file_name <- as.character(hotspots$file_name)

colnames(hotspots)[which(names(hotspots) == "cancer.fractional.overlap")] <- "fc"
colnames(hotspots)[which(names(hotspots) == "immune.fractional.overlap")] <- "fi"
colnames(hotspots)[which(names(hotspots) == "normalised.overlap")] <- "fci"

write.csv(hotspots, file = paste0(folder, "hotspotsSummary.csv"), quote = FALSE)

###################1
#plot IH masks====
maskDir <- "/Users/hzhang/Documents/project/sum/CellPos_all/HE_r/mask_r/"
folder <- "/Users/hzhang/Documents/project/sum/TLS_annotation/CellPos/HE_region/"
Outputdir <- "/Users/hzhang/Documents/project/sum/TLS_annotation/hs_he_region/"
s <- 10
n <- 4
min.d <- 0
max.d <- (n * s * sqrt(2))
library(spdep)
library(EBImage)
ffs <- dir(Outputdir, pattern = '.rdata$')
files = c()
for (i in 1:length(ffs)) {
  files.vec <- c(Outputdir, ffs[i])
  files[i] <- paste(files.vec, collapse = "")
}

#filesN <- list.files(folder, pattern = "\\.csv$")
for (k in 1:length(files)) {
  
  print(paste0('Processing ', k, ' of ', length(files)))
  load(files[k])
  ffs[k] <- gsub("[^[:digit:]]", "",  gsub(".*/", "",  ffs_he[k]))
  
  # also load the masks
  Mask <- readImage(paste0(maskDir, ffs[k], '_HE.ndpi_Mask.jpg'))
  Mask[]
  
  
  s <- 10
  n <- 4
  min.d <- 0
  max.d <- (n * s * sqrt(2))
  
  
  # wid <- 19919      #full res image width
  # hei <- 14390      #full res image height  
  # scalefactor <- 10
  imageformat <- "pdf"
  
  wp <- nrow(Mask)
  hp <- ncol(Mask)
  # mask_r = which(Mask > 0, arr.ind = TRUE)
  # w_r = min(mask_r[,1]) - max(mask_r[,1])
  # h_r = min(mask_r[,2]) -  max(mask_r[,2])
  
  data.naomit <- na.omit(data.wg)
  l.zm.c <- matrix(nrow=w/s, ncol=h/s)
  for (i in 1:nrow(data.naomit)) {
    l.zm.c[data.naomit$x[i]/s,data.naomit$y[i]/s] = data.naomit$G.c[i]
  }
  l.zm.c <- l.zm.c[,ncol(l.zm.c):1]
  
  
  # wid1 <- wid/scalefactor
  # hei1 <- hei/scalefactor
  
  x.i <- 1*(1:(nrow(l.zm.c)))
  y.i <- 1*(1:(ncol(l.zm.c)))
  
  
  # plot lympho hotspots map
  data.naomit <- na.omit(data.wg)
  l.zm.l <- matrix(nrow=w/s, ncol=h/s)
  for (i in 1:nrow(data.naomit)) {
    l.zm.l[data.naomit$x[i]/s,data.naomit$y[i]/s] = data.naomit$G.l[i]
  }
  l.zm.l <- l.zm.l[,ncol(l.zm.l):1]
  
  x.i <- 1*(1:(nrow(l.zm.l)))
  y.i <- 1*(1:(ncol(l.zm.l)))
  
  # plot overlap hotspots map
  data.naomit <- na.omit(data.wg)
  add.gl.gc <- c()
  
  for (f in 1:nrow(data.wg)){
    if(data.wg$G.c[f] >= 3.886 & data.wg$G.l[f] >= 3.886){add.gl.gc[f] <- 2}
    else if (data.wg$G.c[f] != 0 & data.wg$G.l[f] != 0){add.gl.gc[f] <- 1}
    else {add.gl.gc[f] <- 0}
  }
  overlap <- data.wg[,c('x','y')]
  overlap <- cbind(overlap,add.gl.gc)
  
  l.zm.o <- matrix(nrow=w/s, ncol=h/s)
  for (i in 1:nrow(data.naomit)) {
    l.zm.o[data.naomit$x[i]/s,data.naomit$y[i]/s] = overlap$add.gl.gc[i]
  }
  l.zm.o <- l.zm.o[,ncol(l.zm.o):1]
  
  
  x.i <- 1*(1:(nrow(l.zm.o)))
  y.i <- 1*(1:(ncol(l.zm.o)))
  
  dev.new()
  
  #change w and h
  png(file=paste0(Outputdir, ffs[k], ".png"), width = wp, height = hp)
  # a <- image(x.i, y.i,l.zm.o, xlab="",ylab="",bty ="n",xaxt='n', yaxt='n',col=c('white','rosybrown1'), breaks=c(-0.1,0.1,2.2), useRaster = TRUE)
  # a <- image(x.i, y.i,l.zm.c, xlab="",ylab="",bty ="n",xaxt='n', yaxt='n',col=c('red3'), breaks=c(3.886, max(data.naomit$G.c)+10), useRaster = TRUE, add = TRUE)
  # a <- image(x.i, y.i,l.zm.l, xlab="",ylab="",bty ="n",xaxt='n', yaxt='n',
  #            col=c('navyblue'), 
  #            breaks=c(3.886,max(data.naomit$G.l)+10), useRaster = TRUE, add = TRUE)
  
  a <- image(x.i, y.i,l.zm.o, xlab="",ylab="",bty ="n",xaxt='n', yaxt='n',col=c('white','#E4C1C1'), breaks=c(-0.1,0.1,2.2), useRaster = TRUE)
  a <- image(x.i, y.i,l.zm.c, xlab="",ylab="",bty ="n",xaxt='n', yaxt='n',col=c('#C1534E'), breaks=c(3.886, max(data.naomit$G.c)+10), useRaster = TRUE, add = TRUE)
  a <- image(x.i, y.i,l.zm.l, xlab="",ylab="",bty ="n",xaxt='n', yaxt='n',
             col=c('#636AA8'), 
             breaks=c(3.886,max(data.naomit$G.l)+10), useRaster = TRUE, add = TRUE)
  image(x.i, y.i,l.zm.o, xlab="",ylab="",bty ="n",xaxt='n', yaxt='n',col=c('#A0C588'), breaks=c(1.2,2.2), useRaster = TRUE, add = TRUE)
  
  #dev.copy(png,paste0(Outputdir, ffs[f], ".png"), width = wp, height = hp)
  dev.off()
  
}

#resize
# img_t = readImage('/Users/hzhang/Documents/project/sum/CellPos_all/HE_r/hs_he_r/50922.png')
# s=10
# w_r = (dim(img_t)[1] - (floor(dim(img_t)[1]/s) * s)) / (10*s) + 1
# h_r =  (dim(img_t)[2] - (floor(dim(img_t)[2]/s) * s)) / (10*s) + 1
# img_t1 = resize(img_t, w = dim(img_t)[1] * w_r, h = dim(img_t)[2] * h_r)
# a = matrix(nrow = 11/3, ncol = 11/3)

#######################1
#match B, T cell pos to IH====

hs_match <- function(lym_data, data_wg){
  
  lym_data$x.s <- ceiling(lym_data$x / s) * s
  lym_data$y.s <- ceiling(lym_data$y / s) * s
  
  lym_hs <- merge(lym_data, data_wg, by.x = c('x.s', 'y.s'), by.y = c('x', 'y'), all.x = TRUE, all.y=TRUE,sort = TRUE)
  
  cell_count <- aggregate(x ~ class + x.s + y.s, data = lym_hs, length)
  colnames(cell_count)[which(colnames(cell_count) == 'x')] <- 'cell.count'
  #all <- merge(lym_hs, cell_count, by = c('x.s', 'y.s','class'), all.x = TRUE, sort = TRUE)
  # lym_hs$X <- NULL
  # library(tidyverse)
  # td =  lym_hs%>%
  #   group_by(x.s,y.s,class) %>%
  #   mutate(cell_count = n())
  # td$x <- NULL
  # td$y <- NULL
  # td <- unique(td)
  # td_1 = td %>%
  #   group_by(x.s,y.s) %>%
  #   spread(class, sum(cell_count))
  # td[is.na(td)] = 0
  
  cell_count <- unique(merge(cell_count,lym_hs[, !colnames(lym_hs) %in% c('x','y')], by = c('x.s', 'y.s','class'),sort = TRUE))
  
  cell_count$hs <- ifelse(!is.na(cell_count$G.l) & !is.na(cell_count$G.c) & cell_count$G.l >= 3.886 & cell_count$G.c >= 3.886, 'cl',
                          ifelse(!is.na(cell_count$G.c) & cell_count$G.c >= 3.886, 'c',
                                 ifelse(!is.na(cell_count$G.l) & cell_count$G.l >= 3.886, 'l', 'else')))
  return(cell_count)
}


Outputdir <- '/Users/hzhang/Documents/project/sum/CellPos_all/HE_r/hs_he_r/'
t_folder <- "/Users/hzhang/Documents/project/sum/CellPos_all/T/"
b_folder <- '/Users/hzhang/Documents/project/sum/CellPos_all/B/'
ffs_b <- dir(b_folder, pattern = '.csv$', full.names = TRUE)
ffs_t <-  dir(t_folder, pattern = '.csv$', full.names = TRUE)
hs_bt_all <- data.frame()

for (i in 1:length(ffs)){
  print(ffs[i])
  hs_df <- read.csv(paste0(Outputdir, ffs[i], '_cl.csv'))
  df_t <- read.csv(ffs_t[i])
  df_b <- read.csv(ffs_b[i])
  
  df_t <- hs_match(df_t, hs_df)
  df_t$slide <- ffs[i]
  
  df_b <- hs_match(df_b, hs_df)
  df_b$slide <- ffs[i]
  
  write.csv(df_t, paste0(Outputdir, ffs[i], '_hscellcount_t.csv'))
  write.csv(df_b, paste0(Outputdir, ffs[i], '_hscellcount_b.csv'))
  
  df_all <- rbind(df_t, df_b)
  hs_bt_all <- rbind(hs_bt_all, df_all)
}

write.csv(hs_bt_all, '/Users/hzhang/Documents/project/sum/CellPos_all/HE_r/hotspots_r/hs_bt_all.csv')


td <- hs_bt_all %>% group_by(x.s,y.s,slide) %>%
  spread(class, sum(cell.count))

write.csv(td, '/Users/hzhang/Documents/project/sum/CellPos_all/HE_r/hotspots_r/hs_bt_all_spread.csv')

sum_all <- td
for (i in c(4,5,10:19)){
  print(i)
  sum_all[is.na(sum_all[,i]),i] <- 0 
}

write.csv(sum_all, '/Users/hzhang/Documents/project/sum/CellPos_all/HE_r/hotspots_r/sum_all.csv')

#Check if IH is consistent
sum_all_2 <- read.csv('/Users/hzhang/Documents/project/sum/hotspot_bt_DRDIN_foxp3correct/csvs/sum_all_2.csv')
sum_all_2$X.1 <- NULL


sum_all_3 <- merge(sum_all[,c('x.s','y.s','slide','hs','cell.count.c','cell.count.l')], sum_all_2, by = c('x.s','y.s','slide'),all.x = T, all.y = T, sort=T)

write.csv(sum_all_3, '/Users/hzhang/Documents/project/sum/CellPos_all/HE_r/hotspots_r/sum_all_correct.csv')
