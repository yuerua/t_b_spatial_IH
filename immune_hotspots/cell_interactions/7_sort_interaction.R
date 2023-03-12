##############################################################################=
# Sort cell interaction statistics in individual IH 
##############################################################################=

library(plyr)
library(dplyr)

input_dir <- '/Users/hzhang/Documents/project/sum/final/graph/clump_l_cl_TLS_LAG_th'
save_dir <-  '/Users/hzhang/Documents/project/sum/final/graph/clump_l_cl_TLS_LAG_th/sort'

dir.create(save_dir)
ffs <- dir(input_dir, pattern = '.csv$')

df_all <- data.frame()
for (f in ffs){
  print(f)
  df <- read.csv(file.path(input_dir, f))#
  df_all <- rbind(df_all, df)
}

write.csv(df_all, file.path(save_dir, "all_slide_interaction.csv"))

head(df_all)

#interaction features####
df_all$X <- NULL

df_inter <- df_all[, 1:(which(colnames(df_all) == "cd20")-1)]

#exclude p40
df_inter_lym <- df_inter[,!grepl("p40", colnames(df_inter))]
sum <- rowSums(df_inter_lym)
df_inter_lym <- df_inter_lym/sum
df_inter_lym$sum <- sum
df_inter_lym <- cbind(df_inter_lym, df_all[,c('hs', 'clump_id', 'slide', 'clump_area')])

#predict cl/l
df_inter_lym[is.na(df_inter_lym)] <- 0

df_inter_lym$hs_dig <- ifelse(df_inter_lym$hs == 'cl', 1, 0)

colnames(df_inter_lym)
df <- df_inter_lym[, c(1:21, 27)]
logit  <- glm(hs_dig ~ ., data = df, family = "binomial")

sum_df <- summary(logit)

sum_df_coef <- data.frame(sum_df$coefficients)
sum_df_coef$p_adj <- p.adjust(sum_df_coef$Pr...z.., method = 'BH')


saveRDS(sum_df, file.path(save_dir, "inter_sum_df.RData"))
write.csv(df_inter_lym, file.path(save_dir, "df_inter_lym.csv"))
write.csv(sum_df_coef, file.path(save_dir, "inter_sum_df_coef.csv"))

#compute probability
# prob=predict(logit,type=c("response"))
# 
# df$prob=prob
# library(pROC)
# g <- roc(hs_dig ~ prob, data = df)
# plot(g)   


#per cell types inter####
# df_inter_cd8 <- df_inter[,grepl("cd79bCoexp", colnames(df_inter))]
# sum <- rowSums(df_inter_cd8)
# df_inter_cd8 <- df_inter_cd8 / sum
# df_inter_cd8$sum <- sum
# df_inter_cd8 <- cbind(df_inter_cd8, df_all[,c('hs', 'clump_id', 'slide', 'clump_area')])
# 
# df_inter_cd8[is.na(df_inter_cd8)] <- 0
# 
# df_inter_cd8$hs_dig <- ifelse(df_inter_cd8$hs == 'cl', 1, 0)
# 
# colnames(df_inter_cd8)
# df <- df_inter_cd8[, c(1:7, 13)]
# 
# df <- df[,!grepl("p40", colnames(df))]
# logit  <- glm(hs_dig ~ ., data = df, family = "binomial")
# 
# 
# sum_df <- summary(logit)
# 
# sum_df_coef <- data.frame(sum_df$coefficients)
# sum_df_coef$p_adj <- p.adjust(sum_df_coef$Pr...z.., method = 'BH')
# 
# saveRDS(sum_df, file.path(save_dir, "cd79b_inter_sum_df.RData"))
# write.csv(sum_df_coef, file.path(save_dir, "cd79b_inter_sum_df_coef.csv"))



#cell proportion features####

df_prop <- df_all[, c(which(colnames(df_all) == "cd20"):which(colnames(df_all) == "foxp3"))]

sum<- rowSums(df_prop)
df_prop <- df_prop / sum
df_prop$sum <- sum
df_prop <- cbind(df_prop, df_all[,c('hs', 'clump_id', 'slide', 'clump_area')])

df_prop[is.na(df_prop)] <- 0
df_prop$hs_dig <- ifelse(df_prop$hs == 'cl', 1, 0)
colnames(df_prop)

df <- df_prop[, c(1:6, 12)]
df <- df[,!grepl("p40", colnames(df))]

logit  <- glm(hs_dig ~ ., data = df, family = "binomial")

sum_df <- summary(logit)

sum_df_coef <- data.frame(sum_df$coefficients)
sum_df_coef$p_adj <- p.adjust(sum_df_coef$Pr...z.., method = 'BH')

saveRDS(sum_df, file.path(save_dir, "cell_prop_sum_df.RData"))
write.csv(df_prop, file.path(save_dir, "df_prop.csv"))
write.csv(sum_df_coef, file.path(save_dir, "cell_prop_sum_df_coef.csv"))

#cell density features####

df_den <- df_all[, c(which(colnames(df_all) == "cd20"):which(colnames(df_all) == "foxp3"))]

sum<- df_all$clump_area
df_den <- df_den / sum
df_den$sum_clump_area <- sum
df_den <- cbind(df_den, df_all[,c('hs', 'clump_id', 'slide', 'clump_area')])

df_den[is.na(df_den)] <- 0
df_den$hs_dig <- ifelse(df_den$hs == 'cl', 1, 0)
colnames(df_den)

df <- df_den[, c(1:6, 12)]
df <- df[,!grepl("p40", colnames(df))]

logit  <- glm(hs_dig ~ ., data = df, family = "binomial")

sum_df <- summary(logit)

sum_df_coef <- data.frame(sum_df$coefficients)
sum_df_coef$p_adj <- p.adjust(sum_df_coef$Pr...z.., method = 'BH')

saveRDS(sum_df, file.path(save_dir, "den_sum_df.RData"))
write.csv(sum_df_coef, file.path(save_dir, "den_sum_df_coef.csv"))

#df_inter_prop####
df_inter_prop <- cbind(df_inter_lym[,1:(which(colnames(df_inter_lym) == 'sum') - 1)], df_prop[,1:6])
df_inter_prop <- cbind(df_inter_prop, df_all[,c('hs', 'clump_id', 'slide', 'clump_area')])
df_inter_prop[is.na(df_inter_prop)] <- 0
df_inter_prop$hs_dig <- ifelse(df_inter_prop$hs == 'cl', 1, 0)

write.csv(df_inter_prop, file.path(save_dir, "df_inter_prop.csv"))

#cell interaction and proportion features####

df_inter_den <- cbind(df_inter_lym[,1:(which(colnames(df_inter_lym) == 'sum') - 1)], df_den[,1:6])
# sum<- rowSums(df_den)
# df_den <- df_den / sum
# df_den$sum <- sum
df_inter_den <- cbind(df_inter_den, df_all[,c('hs', 'clump_id', 'slide', 'clump_area')])

df_inter_den[is.na(df_inter_den)] <- 0
df_inter_den$hs_dig <- ifelse(df_inter_den$hs == 'cl', 1, 0)
colnames(df_inter_den)

# write.csv(df_inter_den, file.path(save_dir, "df_inter_den.csv"))

df <- df_inter_den[, c(1:27, 32)]
df <- df[,!grepl("p40", colnames(df))]
head(df)
logit  <- glm(hs_dig ~ ., data = df, family = "binomial")

sum_df <- summary(logit)

sum_df_coef <- data.frame(sum_df$coefficients)
sum_df_coef$p_adj <- p.adjust(sum_df_coef$Pr...z.., method = 'BH')

saveRDS(sum_df, file.path(save_dir, "inter_den_sum_df.RData"))
write.csv(sum_df_coef, file.path(save_dir, "inter_den_sum_df_coef.csv"))

#========================#
#LDA
#
head(df_inter_den)

library(MASS)
library(tidyverse)
library(caret)

df <- df_inter_den[,c(1:27, 36)]
model <- lda(hs_dig~., data = df)

#glm net
#try
install.packages("glmnet", repos = "https://cran.us.r-project.org")
library(glmnet)
data(BinomialExample)
x <- BinomialExample$x
y <- BinomialExample$y

fit <- glmnet(x, y, family = "binomial")
cvfit <- cv.glmnet(x, y, family = "binomial", type.measure = "class")
plot(cvfit)
coef(cvfit)
#==========#
#inter + cell densities
#nothing has non-zero coefficient
# write.csv(df_inter_den, "/Users/hzhang/Documents/report/lusc_b/scripts/graph/results/clump_l_cl_without_TLS_LAG/sort/df_inter_den.csv")

df <- df_inter_den[,c(1:27, 32)]

x = as.matrix(df[,1:27])
y = df$hs_dig

fit <- glmnet(x, y, family = "binomial")
cvfit <- cv.glmnet(x, y, family = "binomial", type.measure = "class")

coef(cvfit)

#==========#
#inter + cell densities
#nothing has non-zero coefficient

df_inter_prop <- cbind(df_inter_lym[,1:(which(colnames(df_inter_lym) == 'sum') - 1)], df_prop[,1:6])
df_inter_prop <- cbind(df_inter_prop, df_all[,c('hs', 'clump_id', 'slide', 'clump_area')])
df_inter_prop[is.na(df_inter_prop)] <- 0
df_inter_prop$hs_dig <- ifelse(df_inter_prop$hs == 'cl', 1, 0)

write.csv(df_inter_prop, file.path(save_dir, "df_inter_prop.csv"))

df <- df_inter_prop[,c(1:27, 32)]

colnames(df_inter_prop)
x = as.matrix(df[,1:27])
y = df$hs_dig

fit <- glmnet(x, y, family = "binomial")
cvfit <- cv.glmnet(x, y, family = "binomial", type.measure = "class")

coef(cvfit)


