par(mfrow = c(1,3))
par(pty="s")

setwd("~/R/input_files/OBJ_METHOD_FULL_ANALYSIS_060721")
library("readxl")
library("writexl")
library(parallel)
library(cutpointr)
#####################Functions#######################
df1 <- data.frame(0,0,0,0,0,0,0,0,0,0,0,0)
colnames(df1) <- c("AUC", "optimal_cut", "acc", "sensitivity", 
                "specificity", "TP", "FN", "FP", "TN", "SP", "NON_SP", "TOTAL")

measure_opt <- function(data, col_name1, col_name2, new, previous) {
  
  opt_cut <- cutpointr(data[[col_name1]], 
                       data[[col_name2]], 
                       pos_class = "SP",
                       neg_class = "NON-SP", 
                       method = maximize_metric, 
                       metric = youden,
                       direction = "<=")
  
  npos <- data[data[[col_name2]] == "SP", ] #total SP
  nneg <- data[data[[col_name2]] == "NON-SP", ] #total NON-SP
  
  tp <- opt_cut$sensitivity *  nrow(npos) 
  tn <- opt_cut$specificity *  nrow(nneg)
  fp <- nrow(nneg) - tn
  fn <- nrow(npos) - tp
  sp <- tp + fn
  non_sp <- fp + tn 
  total <- sp + non_sp
  
  library(tidyverse)
  new <- previous %>% 
    add_row(
      AUC = opt_cut$AUC,
      optimal_cut = opt_cut$optimal_cutpoint,
      acc = opt_cut$acc,
      sensitivity = opt_cut$sensitivity, 
      specificity = opt_cut$specificity,
      TP = tp, FN = fn, FP = fp, TN = tn, 
      SP = sp, NON_SP = non_sp, 
      TOTAL = total
    )
  
  #summary(opt_cut)
}


####################FunctionsEND######################



Ecoli_data <- read_excel("SPi_scores_all_Ecoli_SPs_060721.xlsx", sheet = "SPi_score_all_proteins_possible", na = "NA")

#uncomment the following command if inner membrane proteins are to be excluded from our analysis
#Ecoli_data <- Ecoli_data[Ecoli_data$`Localization (STEPdb)` != "B", ]
##################


#uncomment the following command if proteins >1000 are to be eliminated from our analysis. 
Ecoli_data <- Ecoli_data[Ecoli_data$Uniprot_Full_Len <= 1000, ]
###########

#uncomment the following command if proteins of CM quality > 0.75 are considered
#Ecoli_data <- Ecoli_data[Ecoli_data$trRosetta_CM_quality > 0.75, ]
#########


Ecoli_data$Global_index <- 1:nrow(Ecoli_data)
#cyto

N_index <- grep('N', Ecoli_data$`Localization (STEPdb)`)
A_index <- grep('A', Ecoli_data$`Localization (STEPdb)`)
F1_index <- grep('F1', Ecoli_data$`Localization (STEPdb)`)
B_index <- grep('B', Ecoli_data$`Localization (STEPdb)`)
r_index <- grep('r', Ecoli_data$`Localization (STEPdb)`)

Cyto_index <- unique(c(N_index, A_index, F1_index, B_index, r_index))

#lipo
E_index <- grep('E', Ecoli_data$`Localization (STEPdb)`)
I_index <- grep('I', Ecoli_data$`Localization (STEPdb)`)

# EI_index <- unique(c(E_index, I_index))
Lipo_index <- unique(c(E_index, I_index))


#rest
F2_index <- grep('F2', Ecoli_data$`Localization (STEPdb)`)
G_index <- grep('G', Ecoli_data$`Localization (STEPdb)`)
F3_index <- grep('F3', Ecoli_data$`Localization (STEPdb)`)
H_index <- grep('H', Ecoli_data$`Localization (STEPdb)`)
F4_index <- grep('F4', Ecoli_data$`Localization (STEPdb)`)
X_index <- grep('X', Ecoli_data$`Localization (STEPdb)`)


#NS corresponds to non-specific
Peri_NS_index <- unique(c(F2_index, G_index, F3_index, H_index, F4_index, X_index))

#Dataframes:
Cyto_df <- Ecoli_data[Cyto_index, ]

Lipo_df <- Ecoli_data[Lipo_index, ]

tmp1 <- as.data.frame(Ecoli_data[Peri_NS_index, ncol(Ecoli_data)])
tmp2 <- as.data.frame(Ecoli_data[Lipo_index, ncol(Ecoli_data)])

Peri_NS_index_minus_Lipo_index <- which((tmp1[,1] %in% tmp2[,1]) == FALSE)

Peri_NS_df <- Ecoli_data[tmp1[Peri_NS_index_minus_Lipo_index,1], ]


tmp3 <- as.data.frame(tmp1[Peri_NS_index_minus_Lipo_index,1], )
tmp4 <- as.data.frame(Ecoli_data[Cyto_index, ncol(Ecoli_data)])

Peri_minus_Cyto_index <- which((tmp3[,1] %in% tmp4[,1]) == FALSE)

Peri_df <- Ecoli_data[tmp3[Peri_minus_Cyto_index,1], ] 

#Identifying SPs 

Cyto_SPI <- Cyto_df[Cyto_df$SigP5_Prediction == "SP(Sec/SPI)", ]
Peri_SPI <- Peri_df[Peri_df$SigP5_Prediction == "SP(Sec/SPI)", ]

Cyto_SPII <- Cyto_df[Cyto_df$SigP5_Prediction == "LIPO(Sec/SPII)", ]
Lipo_SPII <- Lipo_df[Lipo_df$SigP5_Prediction == "LIPO(Sec/SPII)", ]
Peri_SPII <- Peri_df[Peri_df$SigP5_Prediction == "LIPO(Sec/SPII)", ]

Cyto_TAT <- Cyto_df[Cyto_df$SigP5_Prediction == "TAT(Tat/SPI)", ]
Peri_TAT <- Peri_df[Peri_df$SigP5_Prediction == "TAT(Tat/SPI)", ]


#ROCS of each SP type
#SPI
Peri_SPI$ORIGINAL <- "SP" 
Cyto_SPI$ORIGINAL <- "NON-SP"

SPI_ALL <- rbind(Peri_SPI, Cyto_SPI)
df2 <- measure_opt(SPI_ALL, "rnum2_top5", "ORIGINAL", df2, df1)


#ROC curves
roc_obj <- pROC::roc(response = SPI_ALL$ORIGINAL, 
                     predictor = SPI_ALL$rnum2_top5, 
                     levels = c("SP", "NON-SP"))
roc_df <- data.frame(FPR=1 - roc_obj$specificities, 
                     TPR=roc_obj$sensitivities)

plot(roc_df, type ="l", lwd = 4, 
     col = "skyblue1", 
     main = "Sec-SPI",
     cex.main=2, 
     cex.lab=2, 
     cex.axis=1.5)

abline(coef = c(0,1), lwd=4)



#SPII
Lipo_SPII$ORIGINAL <- "SP" 
Cyto_SPII$ORIGINAL <- "NON-SP"

SPII_ALL <- rbind(Lipo_SPII, Cyto_SPII)
df3 <- measure_opt(SPII_ALL, "rnum2_top5", "ORIGINAL", df3, df2)

roc_obj <- pROC::roc(response = SPII_ALL$ORIGINAL, 
                     predictor = SPII_ALL$rnum2_top5, 
                     levels = c("SP", "NON-SP"))
roc_df <- data.frame(FPR=1 - roc_obj$specificities, 
                     TPR=roc_obj$sensitivities)

plot(roc_df, type ="l", lwd = 4, 
     col = "rosybrown1", 
     main = "Sec-SPII",
     cex.main=2, 
     cex.lab=2, 
     cex.axis=1.5)

abline(coef = c(0,1), lwd=4)


#TAT
Peri_TAT$ORIGINAL <- "SP" 
Cyto_TAT$ORIGINAL <- "NON-SP"

TAT_ALL <- rbind(Peri_TAT, Cyto_TAT)
df4 <- measure_opt(TAT_ALL, "rnum2_top5", "ORIGINAL", df4, df3)

roc_obj <- pROC::roc(response = TAT_ALL$ORIGINAL, 
                     predictor = TAT_ALL$rnum2_top5, 
                     levels = c("SP", "NON-SP"))
roc_df <- data.frame(FPR=1 - roc_obj$specificities, 
                     TPR=roc_obj$sensitivities)

plot(roc_df, type ="l", lwd = 4, 
     col = "green", 
     main = "Tat-SPI",
     cex.main=2, 
     cex.lab=2, 
     cex.axis=1.5)

abline(coef = c(0,1), lwd=4)

#####  Making confusion matirx along with our predictions

##Input - SignalP5 intial data
#True negatives of SignalP5
Cyto_TN <- nrow(Cyto_df[Cyto_df$SigP5_Prediction == "OTHER", ])

#False positives of SignalP5
Cyto_SPI_FP <- nrow(Cyto_df[Cyto_df$SigP5_Prediction == "SP(Sec/SPI)", ])
Cyto_SPII_FP <- nrow(Cyto_df[Cyto_df$SigP5_Prediction == "LIPO(Sec/SPII)", ]) 
Peri_SPII_FP <- nrow(Peri_df[Peri_df$SigP5_Prediction == "LIPO(Sec/SPII)", ])
CYto_TAT_FP <- nrow(Cyto_df[Cyto_df$SigP5_Prediction == "TAT(Tat/SPI)", ]) 
Lipo_SPI_FP <- nrow(Lipo_df[Lipo_df$SigP5_Prediction == "SP(Sec/SPI)", ])
Lipo_TAT_FP <- nrow(Lipo_df[Lipo_df$SigP5_Prediction == "TAT(Tat/SPI)", ])


#True positives of SignalP5
Peri_SPI_TP <- nrow(Peri_df[Peri_df$SigP5_Prediction == "SP(Sec/SPI)", ])
Lipo_SPII_TP <- nrow(Lipo_df[Lipo_df$SigP5_Prediction == "LIPO(Sec/SPII)", ])
Peri_TAT_TP <- nrow(Peri_df[Peri_df$SigP5_Prediction == "TAT(Tat/SPI)", ])

#False negatives of SignalP5
Lipo_SPII_FN <- nrow(Lipo_df[Lipo_df$SigP5_Prediction == "OTHER", ])
Peri_FN <- nrow(Peri_df[Peri_df$SigP5_Prediction == "OTHER", ])

######## Intial INPUT from SignalP5#########################

Cyto_STEPdb <- c(Cyto_TN, Cyto_SPI_FP, Cyto_SPII_FP, CYto_TAT_FP)
Lipo_STEPdb <- c(Lipo_SPII_FN, Lipo_SPI_FP, Lipo_SPII_TP, Lipo_TAT_FP)
Peri_STEPdb <- c(Peri_FN, Peri_SPI_TP, Peri_SPII_FP, Peri_TAT_TP)

Confusion_mat_initial <- rbind(Cyto_STEPdb, Lipo_STEPdb, Peri_STEPdb)
colnames(Confusion_mat_initial) <- c("OTHER", "SPI", "SPII", "SPI")

############## Applying our contact-map appraoch on the SignalP5 matrix#############

CM_pred <- as.matrix(df4[-c(1), ])

#SPI
CM_SPI_TP <- CM_pred[1,6]
CM_SPI_FN <- CM_pred[1,7]
CM_SPI_FP <- CM_pred[1,8]
CM_SPI_TN <- CM_pred[1,9]

#SPII
CM_SPII_TP <- CM_pred[2,6]
CM_SPII_FN <- CM_pred[2,7]
CM_SPII_FP <- CM_pred[2,8]
CM_SPII_TN <- CM_pred[2,9]

#TAT
CM_TAT_TP <- CM_pred[3,6]
CM_TAT_FN <- CM_pred[3,7]
CM_TAT_FP <- CM_pred[3,8]
CM_TAT_TN <- CM_pred[3,9]

#Confusion matrix from the CM approach
Cyto_CM <- c(Cyto_TN + CM_SPI_FP + CM_SPII_FP + CM_TAT_TP, 
             Cyto_SPI_FP - CM_SPI_TN, 
             Cyto_SPII_FP - CM_SPII_TN, 
             CYto_TAT_FP - CM_TAT_TN)

Lipo_CM <- c(Lipo_SPII_FN + CM_SPII_FN, 
             Lipo_SPI_FP, 
             Lipo_SPII_TP - CM_SPII_FN, 
             Lipo_TAT_FP)

Peri_CM <- c(Peri_FN + CM_SPI_FN + CM_TAT_FN,
             Peri_SPI_TP - CM_SPI_FN, 
             Peri_SPII_FP, 
             Peri_TAT_TP - CM_TAT_FN)

Confusion_mat_CM <- rbind(Cyto_CM, Lipo_CM, Peri_CM)
colnames(Confusion_mat_CM) <- c("OTHER", "SPI", "SPII", "SPI")

