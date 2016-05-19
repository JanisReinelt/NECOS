##set workingdirectory 
setwd('/home/raid1/reinelt/scripts/NECOS/r/')

##loading packages
library(ggplot2)
library(Hmisc)
library(gridExtra)
library(cowplot)

############
############

### Load workspace
load("/afs/cbs.mpg.de/projects/nro132_neural-consequences-stress/results/R/workspace/ILM_summary_workspace.RData")

#participants to exclude
#drop_id <- NA
drop_id <- c("NECOS008","NECOS018","NECOS038", "NECOS048", "NECOS059" )

### subset df 
endocrine_df <- subset(ILM_blood_imputed, 
                       select = c("NECOS_ID", "group", "timepoint", 
                                                     "Serum_Copeptin", "Plasma_ACTH", 
                                                     "Serum_Cortisol", "Speichel_Cortisol"),
                       subset = !is.element(NECOS_ID, drop_id))

endocrine_df$rt_timepoint <-   c(20, 135.77 ,173.971,190.425,
                                 206.946,227.784,241.994,250.954,
                                 270.482,286.474,306.717, 322.585, 338.283, 355.212) 

str(endocrine_df)

#ID stress Timepoint value             measure
#1 NECOS001      1        T0    12 VAS 'stressfulness'
#2 NECOS002      1        T0    23 VAS 'stressfulness'
#3 NECOS003      1        T0    19 VAS 'stressfulness'
#4 NECOS004      1        T0    63 VAS 'stressfulness'
#5 NECOS005      0        T0     0 VAS 'stressfulness'
#6 NECOS006      0        T0    21 VAS 'stressfulness'
#  ................................................
#2000 NECOS060      0       T14 1.319 Salivary Cortisol
#2001 NECOS061      0       T14 1.730 Salivary Cortisol
#2002 NECOS056      1       T14 7.008 Salivary Cortisol
#2003 NECOS055      0       T14 3.449 Salivary Cortisol
#2004 NECOS057      0       T14 1.294 Salivary Cortisol
#2005 NECOS062      1       T14 0.767 Salivary Cortisol
#2006 NECOS063      0       T14 1.087 Salivary Cortisol
#2007 NECOS064      0       T14 1.363 Salivary Cortisol
#2008 NECOS065      1       T14 1.388 Salivary Cortisol
#2009 NECOS066      1       T14 1.683 Salivary Cortisol
#2010 NECOS067      0       T14 0.946 Salivary Cortisol


##create AUC function measure1+measure2/betweentime/2
auc <- function(m1,m2,t1,t2){((m1+m2)*(t2-t1))/2}


#create empty df
aucdf_c <- df <- read.table(text = "",col.names = paste("interval",1:16,sep=""))
x <- 1:16
for (k in 3:69)
{
  for ( i in 1: nrow(cort_w))
  { x[i]<- mapply(auc, cort_w[i,k],cort_w[(i+1),k],cort_w[i,"Time"],cort_w[(i+1),"Time"])
  
  }
  aucdf_c[k,] <- x
}

aucdf_c <- aucdf_c[-(1:2),-(15:16)]
#interval1 ~ T0 - T1, etc.

##calculate AUCg
aucdf_c$AUCg <- apply(aucdf_c,1,sum)

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


