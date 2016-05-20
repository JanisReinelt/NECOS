## set workingdirectory 
setwd('/home/raid1/reinelt/scripts/NECOS/r/')

## load packages
library(ggplot2)
library(Hmisc)
library(gridExtra)
library(cowplot)
library(zoo)
library(reshape2)

############
############

## Load workspace
load("/afs/cbs.mpg.de/projects/nro132_neural-consequences-stress/results/R/workspace/ILM_summary_workspace.RData")

## participants to exclude
#drop_id <- NA
drop_id <- c("NECOS008","NECOS038", "NECOS048", "NECOS059" )

## get group information from meta data --> master_sheet_google & create group vector
master_sheet_google <- read.csv(file="/scr/nil1/reinelt/NeCoS/master_sheet_google.csv", sep=",", header=T)
stress <- master_sheet_google[master_sheet_google$Stress==1,]$NECOS_ID
control  <- master_sheet_google[master_sheet_google$Stress==0,]$NECOS_ID

## subset df, choose measures, add names below!!!
endocrine_df <- subset(ILM_blood_imputed, 
                       select = c("NECOS_ID", "group", "timepoint", 
                                  "Serum_Copeptin", "Plasma_ACTH", 
                                  "Serum_Cortisol", "Speichel_Cortisol",
                                  "STAI"),
                       subset = !is.element(NECOS_ID, drop_id))
##########################################
### names of measures that are choosen ###
###      need to be added HERE!!!!     ###
##########################################
colnames(endocrine_df)[4:ncol(endocrine_df)] <- c("Copeptin","ACTH","Ser_Cort","Sal_Cort", "STAI")
number_measures <- length(colnames(endocrine_df)[4:ncol(endocrine_df)])

## add real time column
rt_timepoint <-   c(20, 135.77 ,173.971,190.425,
                                 206.946,227.784,241.994,250.954,
                                 270.482,286.474,306.717, 322.585, 338.283, 355.212) 
endocrine_df <- cbind(endocrine_df[,1:3], rt_timepoint , endocrine_df[,4:ncol(endocrine_df)])

str(endocrine_df)

## create empty auc_df 
auc_df <- data.frame(NECOS_ID = sort(factor(rep(unique(endocrine_df$NECOS_ID), 13))),
                     group = NA,
                     time_interval = factor(rep(c('delta_1_2','delta_2_3','delta_3_4',
                                              'delta_4_5','delta_5_6','delta_6_7',
                                              'delta_7_8','delta_8_9','delta_9_10',
                                              'delta_10_11','delta_11_12','delta_12_13',
                                              'delta_13_14'), length(unique(endocrine_df$NECOS_ID)))),
                     Copeptin = NA,
                     ACTH = NA,
                     Ser_Cort = NA,
                     Sal_Cort = NA)

## add group column to auc_df 
auc_df$group <- ifelse(is.element(auc_df$NECOS_ID, stress), 1, 
                ifelse(is.element(auc_df$NECOS_ID, control), 0, NA))
auc_df$group <- factor(auc_df$group, levels = c(0,1), labels = c("control", "stress"))


## create ID vector for df splitting
IDs <-  unique(endocrine_df$NECOS_ID)

## create func to add consecutive rows of a given vector
add_consec_rows <- function(y){rollapply(y, 2, sum)} 

## create 
measures <- as.factor(colnames(endocrine_df[5:ncol(endocrine_df)]))

## calculate time deltas from endocrine_df splitted by NECOS_ID
time_row_diffs <- sapply(IDs, function(x){diff(endocrine_df[grep(x, endocrine_df$NECOS_ID), "rt_timepoint"])})

## iterate through measures, calculate AUC, and add it to auc_df
for(imeasure in measures){
    print(imeasure)

    ## apply add_consec_rows to endocrine_df splitted by NECOS_ID
    measure_row_sums <- sapply(IDs, function(x){add_consec_rows(endocrine_df[grep(x, endocrine_df$NECOS_ID), imeasure])})
    
    ## calculate auc with trapezoid formula (((m1+m2)*(t2-t1))/2)
    auc_measure <- as.data.frame(measure_row_sums * time_row_diffs)/2
    
    #colnames(auc_copep) <- IDs
    
    auc_measure$time_interval <- as.factor(c('delta_1_2','delta_2_3','delta_3_4',
                                           'delta_4_5','delta_5_6','delta_6_7',
                                            'delta_7_8','delta_8_9','delta_9_10',
                                            'delta_10_11','delta_11_12','delta_12_13',
                                            'delta_13_14'))
    
    auc_measure_l <- melt(auc_measure, id.vars=c("time_interval"))
    
    auc_df[,imeasure] <-  auc_measure_l$value
  
}

## rename colums of auc_df
colnames(auc_df)[4:ncol(auc_df)] <-  c(paste(measures, '_AUC', sep=''))

colnames(auc_df[4:7])  <-  c(paste(measures, '_AUC', sep=''))

colnames(auc_df)

## create auc_base 
auc_base <- data.frame(NECOS_ID = IDs,
                     group = NA)
## add group column to auc_df 
auc_base$group <- ifelse(is.element(auc_base$NECOS_ID, stress), 1, 
                       ifelse(is.element(auc_base$NECOS_ID, control), 0, NA))
auc_base$group <- factor(auc_base$group, levels = c(0,1), labels = c("control", "stress"))

## create empty auc_g df
measures_auc_g <- as.factor(c(paste(colnames(auc_df[4:ncol(auc_df)]), '_g', sep=''),
                              paste(colnames(auc_df[4:ncol(auc_df)]), '_g_1-5', sep=''),
                              paste(colnames(auc_df[4:ncol(auc_df)]), '_g_5-14', sep='')))
auc_g <- auc_base[,1:2]
for(i in measures_auc_g){
    auc_g[,i]<- NA
}

## calculate AUCg (area under the curve in respect to ground, including all intervals)
for(imeasure in measures_auc_g[1:number_measures]){
    print(imeasure)
    # AUCg 
    par <- sapply(IDs, function(x){sum(auc_df[grep(x, auc_df$NECOS_ID), as.character(strsplit(imeasure, '_g'))])})
    auc_g[,imeasure] <- par
}

## calculate AUCg from time point 1 - 5 
t_1_5 <- as.factor(c("delta_1_2", "delta_2_3", "delta_3_4", "delta_4_5"))

for(imeasure in measures_auc_g[(number_measures+1): (2*number_measures)]){
    print(imeasure)
     # AUCg 1-5
     temp_auc_df <- subset(auc_df, subset = is.element(time_interval, t_1_5))
     par <- sapply(IDs, function(x){sum(temp_auc_df[grep(x, temp_auc_df$NECOS_ID), as.character(strsplit(imeasure, '_g_1-5'))])})
     auc_g[,imeasure] <- par
}

## calculate AUCg from time point 6 - 14 
for(imeasure in measures_auc_g[((2*number_measures)+1) : length(measures_auc_g)]){
    print(imeasure)
    #print(strsplit(imeasure, '_6-14'))
    #print(gsub('_6-14', '_1-5',imeasure))
    auc_g[,imeasure] <- auc_g[,as.character(strsplit(imeasure, '_5-14'))] - auc_g[,as.character(gsub('_5-14', '_1-5',imeasure))]
}




## create empty auc_i df
measures_auc_i <- as.factor(paste(colnames(auc_df[4:ncol(auc_df)]), '_i', sep=''))
auc_i <- auc_base[,1:2]
for(i in measures_auc_i){
    auc_i[,i]<- NA
}

## rt_delta (real time points t14 - t5 --> to calculate rectangel)
rt_delta <- unique(endocrine_df$rt_timepoint[endocrine_df$timepoint == 14] 
                    - endocrine_df$rt_timepoint[endocrine_df$timepoint == 5])



for(imeasure in colnames(endocrine_df[5:ncol(endocrine_df)])){
    rect <- subset(endocrine_df, select = imeasure, subset = timepoint == 5) * rt_delta)

}




############test which timepoints differ significantly from baseline for area under the curve calculation##############
arithmean <- function(x,y){(x+y)/2}
all_w$cortbaseline <- all_w$`T5_Salivary Cortisol`
stress_cortbaseline <- all_w$cortbaseline[all_w$stress==1]


#Test until which timepoint in TSST group values differ from baseline for AUCi

wilcox.test(all$value[all$Timepoint=="T6"&all$stress==1&all$measure=="Salivary Cortisol"],stress_cortbaseline,paired=T)

wilcox.test(all$value[all$Timepoint=="T7"&all$stress==1&all$measure=="Salivary Cortisol"],stress_cortbaseline,paired=T)

wilcox.test(all$value[all$Timepoint=="T8"&all$stress==1&all$measure=="Salivary Cortisol"],stress_cortbaseline,paired=T)

wilcox.test(all$value[all$Timepoint=="T9"&all$stress==1&all$measure=="Salivary Cortisol"],stress_cortbaseline,paired=T)

wilcox.test(all$value[all$Timepoint=="T10"&all$stress==1&all$measure=="Salivary Cortisol"],stress_cortbaseline,paired=T)

wilcox.test(all$value[all$Timepoint=="T11"&all$stress==1&all$measure=="Salivary Cortisol"],stress_cortbaseline,paired=T)

wilcox.test(all$value[all$Timepoint=="T12"&all$stress==1&all$measure=="Salivary Cortisol"],stress_cortbaseline,paired=T)

wilcox.test(all$value[all$Timepoint=="T13"&all$stress==1&all$measure=="Salivary Cortisol"],stress_cortbaseline,paired=T)

wilcox.test(all$value[all$Timepoint=="T14"&all$stress==1&all$measure=="Salivary Cortisol"],stress_cortbaseline,paired=T)
##cortisol is significantly diff. from baseline in stressgroup until T13

#########################

##remove timepoints where cort in stress group is not significantly different from baseline /CIs do not overlap#
aucidf_c <- aucdf_c
aucidf_c <- aucidf_c[,-c(1:5,14)]
###sum area under the curve for points of interest
aucidf_c$AUCp <- apply(aucidf_c,1,sum)
#correct for baseline
#row 6 = T5
#row 15 = T14
basecorr <- function(x,t14,t5){(x*(t14-t5))}
for (k in 3:69)
{
  aucidf_c[k,"base"]<- mapply(basecorr,cort_w[6,k],cort_w[15,"Time"],cort_w[6,"Time"])
  
}

#shift rows
aucidf_c[1:67,"base"] <- aucidf_c[3:69,"base"]
#remove additional rows
aucidf_c <- aucidf_c[-(68:69),]

aucidf_c$auci <- aucidf_c$AUCp-aucidf_c$base
aucdf_c$stress <- all_w$stress
aucidf_c$stress <- all_w$stress
aucidf_c$condition <- as.character(aucidf_c$stress)
aucidf_c$condition[aucidf_c$condition==1] <- "TSST"
aucidf_c$condition[aucidf_c$condition==0] <- "Placebo-TSST"
aucidf_c$ID <- cort_w$ID


#statistical analysis of AUCg & AUCi


wilcox.test(aucdf_c$AUCg~aucdf_c$stress,paired=F)


###...and AUCi
wilcox.test(aucidf_c$auci~aucidf_c$stress,paired=F)

aucbp <- ggplot(aucidf_c,aes(x=condition,y=AUCg,group=condition))+
  geom_boxplot(colour="black",aes(fill=condition),alpha=0.5,outlier.colour = temp[4],outlier.size = 1.5)+
  ylim(c(0,4200))+
  geom_point(colour="gray32",size=1)+
  scale_fill_manual(values=temp[c(4,5)])+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))+
  ggtitle("Area under the Curve with respect to ground (AUCg)Saliva Cortisol")+
  ylab("AUCg nmol/l * min")+
  annotate(geom="segment", x = 1.2, xend = 1.8, y = 3400,yend=3400,colour="black",alpha=0.7)+
  annotate(geom="text",x=1.5,y=3500,label="**",size=8)

aucbp


aucibp <- ggplot(aucidf_c,aes(x=condition,y=auci,group=condition))+
  geom_boxplot(colour="black",aes(fill=condition,outlier.colour = condition),alpha=0.5)+
  ylim(c(0,4500))+
  geom_point(colour="gray32",size=1)+
  scale_fill_manual(values=temp[c(4,5)])+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))+
  ggtitle("Area under the Curve with respect to increase (AUCi) Saliva Cortisol")+
  ylab("AUCi nmol/l * min")+
  annotate(geom="segment", x = 1.2, xend = 1.8, y = 3300,yend=3300,colour="black",alpha=0.7)+
  annotate(geom="text",x=1.5,y=3500,label="***",size=8)

aucibp


