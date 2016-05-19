##set workingdirectory 
setwd('/home/raid1/reinelt/github/NECOS/r/')

##loading packages
library(ggplot2)
library(Hmisc)
library(gridExtra)
library(cowplot)

############
############
####            set summary file and specify parameter to plot

#read in 
master.frame_resting.state <- read.csv(file="/scr/nil1/reinelt/NeCoS/R_data/master.frame_resting.state_movement.csv")
master.frame_resting.state$scan_id <- factor(master.frame_resting.state$scan_id)

#Mandy`s outlier function
add.outlier <- function(p,labvar = as.character(p$mapping$y)){
    df <- data.frame(y = with(p$data,eval(p$mapping$y)),
                     x = with(p$data,eval(p$mapping$x)))
    
    df.l <- split(df,df$x)
    
    mm <- Reduce(rbind, lapply(df.l,FUN = function(df){
        data.frame(y = df$y[df$y <= (quantile(df$y)[2] - 1.5 * IQR(df$y)) | df$y >= (quantile(df$y)[4] + 1.5 * IQR(df$y))],
                   x = df$x[df$y <= (quantile(df$y)[2] - 1.5 * IQR(df$y)) | df$y >= (quantile(df$y)[4] + 1.5 * IQR(df$y))]
        )})
    )
    
    
    mm$x <- factor(mm$x,levels=sort(as.numeric(as.character(unique(p$data[,as.character(p$mapping$x)])))),
                   labels = levels(p$data[,as.character(p$mapping$x)])
    )
    
    names(mm) <- c(as.character(p$mapping$y),as.character(p$mapping$x))
    mm <- merge(p$data[,c(names(mm),labvar)],mm)
    
    p + geom_text(data=mm,
                  aes_string(label=labvar),
                  vjust = -0.5)
}

#parameter to plot
parameter <- "Mean.FD"

#participants to exclude
#drop_id <- NA
drop_id <- c("NECOS008", "NECOS018","NECOS038","NECOS048","NECOS059")

#############
#prepare dfs
master.frame_drop <- subset(master.frame_resting.state, select = c("NECOS_ID","scan_id", "group", parameter), subset = !is.element(NECOS_ID, drop_id))
colnames(master.frame_drop)[4] <- "par"

df_stress   <- subset(master.frame_drop, select = c("NECOS_ID","scan_id", "group", "par"), subset= master.frame_resting.state$group == "stress")
#colnames(df_stress)[4] <- "par"

df_control  <- subset(master.frame_drop, select = c("NECOS_ID","scan_id", "group", "par"), subset= master.frame_resting.state$group == "control")
#colnames(df_control)[4] <- "par"

##boxplots
#by group, with outliers
boxplot_stress   <- ggplot(df_stress, aes(x= scan_id, y= par))+
                            background_grid(major = "xy", minor = "none", size.major = 1)+
                            geom_boxplot()+
                            #scale_y_continuous(limits = c(0, max(stress_df$Serum_Copeptin, na.rm=T)+max(summary_copep_imputed$sd, na.rm=T)))+
                            #geom_text(aes(label=ifelse(fbox>2*stress_summ$sd, NECOS_ID,""), hjust=1.1))+
                            labs(title= paste(parameter, "stress group", sep=" "))+
                            ylab(parameter)+
                            xlab("scan_id")+
                            stat_summary(fun.y=mean, colour="darkred", geom="point", 
                                        shape=18, size=3,show.legend = FALSE)+
                            scale_y_continuous(limits=c(0, max(master.frame_drop$par)))

add.outlier(boxplot_stress, "NECOS_ID")


if(parameter=="sd_FD"){
    df_control <- subset(df_control, subset= df_control$NECOS_ID != "NECOS008" ) 
} 
boxplot_control  <- ggplot(df_control, aes(x= scan_id, y= par))+
                            background_grid(major = "xy", minor = "none", size.major = 1)+
                            geom_boxplot()+
                            #scale_y_continuous(limits = c(0, max(stress_df$Serum_Copeptin, na.rm=T)+max(summary_copep_imputed$sd, na.rm=T)))+
                            #geom_text(aes(label=ifelse(fbox>2*stress_summ$sd, NECOS_ID,""), hjust=1.1))+
                            labs(title= paste(parameter, "control group", sep=" "))+
                            ylab(parameter)+
                            xlab("scan_id")+
                            stat_summary(fun.y=mean, colour="darkred", geom="point", 
                                        shape=18, size=3,show.legend = FALSE)+
                            scale_y_continuous(limits=c(0, max(master.frame_drop$par)))

add.outlier(boxplot_control, "NECOS_ID")

##violine plot
#by group, with outliers
violine_plot_stress <- ggplot(df_stress, aes(x=scan_id, y=par))+
                            geom_violin()+
                            #geom_dotplot(binaxis='y', stackdir='center', size=10 ,dotsize=1)+
                            geom_jitter(shape=16, position=position_jitter(0.1))+
                            #geom_boxplot(width=0.1)+
                            labs(title= paste(parameter, "stress group", sep=" "))+
                            ylab(parameter)+
                            xlab("scan_id")+
                            stat_summary(fun.y=mean, colour="darkred", geom="point", 
                                         shape=18, size=8 ,show.legend = TRUE)+
                            stat_summary(fun.y=median, colour="darkred", geom="point", 
                                         shape=10, size=8 ,show.legend = TRUE)+
                            scale_y_continuous(limits=c(0, max(master.frame_drop$par)))

add.outlier(violine_plot_stress,"NECOS_ID")


if(parameter=="sd_FD"){
    df_control <- subset(df_control, subset= df_control$NECOS_ID != "NECOS008" ) 
} 
violine_plot_control <-ggplot(df_control, aes(x=scan_id, y=par))+
                            geom_violin()+
                            #geom_dotplot(binaxis='y', stackdir='center', size=10 ,dotsize=1)+
                            geom_jitter(shape=16, position=position_jitter(0.1))+
                            labs(title= paste(parameter, "control group", sep=" "))+
                            ylab(parameter)+
                            xlab("scan_id")+
                            stat_summary(fun.y=mean, colour="darkred", geom="point", 
                                         shape=18, size=8 ,show_guide = FALSE)+
                            stat_summary(fun.y=median, colour="darkred", geom="point", 
                                         shape=10, size=8 ,show_guide = FALSE)+
                            scale_y_continuous(limits=c(0, max(master.frame_drop$par)))

add.outlier(violine_plot_control,"NECOS_ID")


#violine plots
#by group, with outliers

violine_jitter <- ggplot(master.frame_drop, aes(x=scan_id, y=par, fill=group))+
                    #theme_gray()+
                    background_grid(major = "y", minor = "none", size.major = .3 ,colour.major = "black")+
                    scale_fill_manual("Group",
                                        values=c(  "green","red"),
                                        labels=c("control, n=", "stressor, n="))+
                    geom_violin(position=position_dodge(.75))+
                    geom_dotplot(binaxis='y', stackdir='center', binwidth = .01 ,position=position_dodge(.75), dotsize=.5)+
                    labs(title= parameter)+
                    ylab(parameter)+
                    xlab("scan_id")+
                    scale_y_continuous(limits=c(0, max(master.frame_drop$par)))
                
violine_box <- ggplot(master.frame_drop, aes(x=scan_id, y=par, fill=group))+
                #theme_gray()+
                background_grid(major = "y", minor = "none", size.major = .3 ,colour.major = "black")+
                scale_fill_manual("Group",
                                  values=c(  "green","red"),
                                  labels=c("control, n=", "stressor, n="))+
                geom_violin(position=position_dodge(.75))+
                geom_boxplot(width=0.2, position=position_dodge(.75))+
                stat_summary(fun.y=mean, colour="darkred", geom="point", 
                             shape=18, size=2, position=position_dodge(.75), show_guide = TRUE)+
                labs(title= parameter)+
                ylab(parameter)+
                xlab("scan_id")+
                scale_y_continuous(limits=c(0, max(master.frame_drop$par)))


plot_grid(violine_jitter, violine_box, labels=c("A", "B"), ncol = 1, nrow = 2)
#############


ggplot(master.frame_resting.state[master.frame_resting.state$NECOS_ID=="NECOS059",], aes(sd_FD))+
    geom_histogram()


ggplot(master.frame_resting.state, aes(x=n_fd_over_1., y=sd_FD))+
    geom_point()+
    geom_smooth(method=lm)

cor.test(master.frame_resting.state$n_fd_over_1., master.frame_resting.state$sd_FD)


ggplot(master.frame_resting.state, aes(x=n_fd_over_1., y=Mean.FD))+
    geom_point()+
    geom_smooth(method=lm)

cor.test(master.frame_resting.state$n_fd_over_1., master.frame_resting.state$Mean.FD)

ggplot(master.frame_resting.state, aes(x=sd_FD, y=Mean.FD))+
    geom_point()+
    geom_smooth(method=lm)

cor.test(master.frame_resting.state$sd_FD, master.frame_resting.state$Mean.FD)





add.outlier(violine_box,"NECOS_ID")





df_meanFD_stress$NECOS_ID[df_meanFD_stress$scan_id==5 & df_meanFD_stress$Mean.FD > 0.2]



################
    
master_bis_s <- subset(master_summary, rest_num < 3) 
master_bis_s$rest_num <- factor(master_bis_s$rest_num)
str(master_bis_s)

master_bis_s_anov = ezANOVA(data = master_bis_s, dv = Mean.FD, wid = NECOS_ID, 
                    within = .(rest_num), between = group,
                    type=3)
master_bis_s_anov
#####################
master_nach_s <- subset(master_summary, rest_num >= 3) 
master_nach_s$rest_num <- factor(master_nach_s$rest_num)
str(master_nach_s)

master_nach_s = ezANOVA(data = master_nach_s, dv = Mean.FD, wid = NECOS_ID, 
                            within = .(rest_num), between = group,
                            type=3)
master_nach_s
#################
master_2_6 <- subset(master_summary, rest_num > 1) 
master_2_6$rest_num <- factor(master_2_6$rest_num)
str(master_2_6)

master_2_6 = ezANOVA(data = master_2_6, dv = Mean.FD, wid = NECOS_ID, 
                        within = .(rest_num), between = group,
                        type=3)
master_2_6
###################
master_1_5 <- subset(master_summary, rest_num < 6) 
master_1_5$rest_num <- factor(master_1_5$rest_num)
str(master_1_5)

master_1_5 = ezANOVA(data = master_1_5, dv = Mean.FD, wid = NECOS_ID, 
                     within = .(rest_num), between = group,
                     type=3)
master_1_5
#######################
master_2_3 <- subset(master_summary, rest_num > 1 & rest_num < 4 ) 
master_2_3$rest_num <- factor(master_2_3$rest_num)
str(master_2_3)

master_2_3 = ezANOVA(data = master_2_3, dv = Mean.FD, wid = NECOS_ID, 
                     within = .(rest_num), between = group,
                     type=3)
master_2_3



t.test(stress_summary$Mean.FD ~ stress_summary$rest_num[stress_summary$rest_num==2 && stress_summary$rest_num==2])

Wilcox_test_meanFD <- as.list(c(1:6))
t_test_meanFD <- as.list(c(1:6))
for(i in 1:6){
    Wilcox_test_meanFD[[i]] <- wilcox.test(Mean.FD[rest_num == i]~ group[rest_num == i], data = master_summary)
    t_test_meanFD[[i]] <- t.test(Mean.FD[rest_num == i]~ group[rest_num == i], data = master_summary)
}


master_summary_rest1 <-  subset(master_summary, subset= rest_num==1)
master_summary_rest2 <-  subset(master_summary, subset= rest_num==2)
master_summary_rest3 <-  subset(master_summary, subset= rest_num==3)
master_summary_rest4 <-  subset(master_summary, subset= rest_num==4)
master_summary_rest5 <-  subset(master_summary, subset= rest_num==5)
master_summary_rest6 <-  subset(master_summary, subset= rest_num==6)


ggplot(master_summary_rest6, aes(x=Mean.FD, fill=group)) + geom_density(alpha=.3)

boxplot( master_summary$Mean.FD[master_summary$group==1] ~ master_summary$rest_num[master_summary$group==1], col="red")
boxplot( master_summary$Mean.FD[master_summary$group==0] ~ master_summary$rest_num[master_summary$group==0], col="green")

boxplot(Mean.FD ~ rest_num * group, data=master_summary, notch=TRUE,
        col=(c("gold","darkgreen")),
        main="mean_FD", xlab="rest and group") 

lines(master_summary$Mean.FD[master_summary$group==0], col="green")
lines(scale(serc), col="green")



Wilcox_test_meanFD

