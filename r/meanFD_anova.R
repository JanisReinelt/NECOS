setwd('/home/raid1/reinelt/scripts/NECOS/r/')


library(ez)
library(ggplot2)

master.frame_resting.state <- read.csv(file="/scr/nil1/reinelt/NeCoS/R_data/master.frame_resting.state_movement.csv")
master.frame_resting.state$scan_id <- factor(master.frame_resting.state$scan_id)

master_sheet_google <- read.csv(file="/scr/nil1/reinelt/NeCoS/master_sheet_google.csv", sep=",", header=T)
stress <- master_sheet_google[master_sheet_google$Stress==1,]$NECOS_ID
control  <- master_sheet_google[master_sheet_google$Stress==0,]$NECOS_ID


#parameter to plot
parameter <- "Mean.FD"

#participants to exclude
#drop_id <- NA
drop_id <- c("NECOS008","NECOS018","NECOS038", "NECOS048", "NECOS059" )

#############
#prepare dfs
master.frame_drop <- subset(master.frame_resting.state, select = c("NECOS_ID","scan_id", "group", parameter), subset = !is.element(NECOS_ID, drop_id))
master.frame_drop$NECOS_ID <- factor(master.frame_drop$NECOS_ID)

colnames(master.frame_drop)[4] <- "par"
######

######
meanFD_rest_groups <- data.frame(NECOS_ID = sort(factor(rep(unique(master.frame_drop$NECOS_ID), 3))),
                                 rs_block = factor(rep(c('block_1','block_2','block_3'),62)),
                                 meanFD = NA)

meanFD_rest_groups$group <- ifelse(is.element(meanFD_rest_groups$NECOS_ID, stress), 1, 
                                           ifelse(is.element(meanFD_rest_groups$NECOS_ID, control), 0, NA))
meanFD_rest_groups$group <- factor(meanFD_rest_groups$group, levels = c(0,1), labels = c("control", "stress"))

ids <- unique(master.frame_drop$NECOS_ID)
length(ids)

for (i in ids){
    #rest_group1
    meanFD_rest_groups$meanFD[meanFD_rest_groups$NECOS_ID == i 
                              & meanFD_rest_groups$rs_block == 'block_1' ]  <-   mean(master.frame_drop$par[master.frame_drop$NECOS_ID == i 
                                                                                                                & master.frame_drop$scan_id == c(1,2)])
    #rest_group2
    meanFD_rest_groups$meanFD[meanFD_rest_groups$NECOS_ID == i 
                              & meanFD_rest_groups$rs_block == 'block_2' ]  <-   mean(master.frame_drop$par[master.frame_drop$NECOS_ID == i 
                                                                                                                & master.frame_drop$scan_id == c(3,4)])
    
    #rest_group3
    meanFD_rest_groups$meanFD[meanFD_rest_groups$NECOS_ID == i 
                              & meanFD_rest_groups$rs_block == 'block_3' ]  <-   mean(master.frame_drop$par[master.frame_drop$NECOS_ID == i 
                                                                                                                & master.frame_drop$scan_id == c(5,6)])
    
}
################################

###################
###################
##### Plot ########
###################
###################

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

meanFD_boxplot <- ggplot(meanFD_rest_group1_2, aes(x=rs_block, y=meanFD, colour = group))+
    geom_boxplot()


str(master.frame_drop)
str(meanFD_rest_groups)

violine_box <- ggplot(meanFD_rest_groups, aes(x=rs_block, y=meanFD, fill=group))+
    #theme_gray()+
    background_grid(major = "y", minor = "none", size.major = .3 ,colour.major = "black")+
    scale_fill_manual("Group",
                      values=c(  "green","red"),
                      labels=c("control", "stressor"))+
    geom_violin(position=position_dodge(.85))+
    geom_boxplot(width=0.2, position=position_dodge(.85))+
    stat_summary(fun.y=mean, colour="darkred", geom="point", 
                 shape=18, size=2, position=position_dodge(.85))+
    #labs(title= mean FD)+
    ylab("mean FD")+
    xlab("scan blocks")+
    scale_y_continuous(limits=c(0, max(master.frame_drop$par)))


###############
###############

############################
###### analyse motion ######
############################
meanFD_anov = ezANOVA(data = meanFD_rest_groups, dv = meanFD, wid = NECOS_ID, 
                      within = .(rs_block), between = group,
                      type=3)
meanFD_anov
str(meanFD_rest_groups)

t.test(meanFD_rest_groups$meanFD[meanFD_rest_groups$rs_block == 'block_1' &  meanFD_rest_groups$group == 'stress'],
           meanFD_rest_groups$meanFD[meanFD_rest_groups$rs_block == 'block_1' &  meanFD_rest_groups$group == 'control'])

t.test(meanFD_rest_groups$meanFD[meanFD_rest_groups$rs_block == 'block_2' &  meanFD_rest_groups$group == 'stress'],
       meanFD_rest_groups$meanFD[meanFD_rest_groups$rs_block == 'block_2' &  meanFD_rest_groups$group == 'control'])

t.test(meanFD_rest_groups$meanFD[meanFD_rest_groups$rs_block == 'block_3' &  meanFD_rest_groups$group == 'stress'],
       meanFD_rest_groups$meanFD[meanFD_rest_groups$rs_block == 'block_3' &  meanFD_rest_groups$group == 'control'])


ggplot(meanFD_rest_groups, aes(x=rs_block, y=meanFD, colour = group))+
    geom_boxplot()
################################

meanFD_rest_group1_2 <- subset(meanFD_rest_groups, subset = rs_block != 'block_3')

meanFD_group1_2_anov = ezANOVA(data = meanFD_rest_group1_2, dv = meanFD, wid = NECOS_ID, 
                      within = .(rs_block), between = group,
                      type=3)
meanFD_group1_2_anov


################################

meanFD_rest_groups_wide <- reshape(meanFD_rest_groups, v.names = 'meanFD', idvar = 'NECOS_ID', timevar = 'rs_block', direction = 'wide')
meanFD_rest_groups_wide$group <- factor(meanFD_rest_groups_wide$group) 

#str(meanFD_rest_groups_wide)

meanFD_rest_groups_wide$delta_rest_gr2_gr1 <-  meanFD_rest_groups_wide$meanFD.block_2 - meanFD_rest_groups_wide$meanFD.block_1 
#meanFD_rest_groups_wide$demeaned_del_rest_gr2_gr1 <- meanFD_rest_groups_wide$delta_rest_gr2_gr1 - mean(meanFD_rest_groups_wide$delta_rest_gr2_gr1)

df_stress   <- subset(meanFD_rest_groups,  subset = meanFD_rest_groups$group == "stress")
df_stress_wide <- reshape(df_stress, v.names = 'meanFD', idvar = 'NECOS_ID', timevar = 'rs_block', direction = 'wide')

df_stress_wide$demeaned_rest_gr1 <- df_stress_wide$meanFD.block_1 - mean(meanFD_rest_groups_wide$meanFD.block_1)
df_stress_wide$demeaned_rest_gr2 <- df_stress_wide$meanFD.block_2 - mean(meanFD_rest_groups_wide$meanFD.block_2)
df_stress_wide$demeaned_rest_gr3 <- df_stress_wide$meanFD.block_3 - mean(meanFD_rest_groups_wide$meanFD.block_3)

df_stress_wide$delta_rest_gr2_gr1 <-  df_stress_wide$meanFD.block_2 - df_stress_wide$meanFD.block_1 
df_stress_wide$demeaned_del_rest_gr2_gr1 <- df_stress_wide$delta_rest_gr2_gr1 - mean(meanFD_rest_groups_wide$delta_rest_gr2_gr1)


df_control  <- subset(meanFD_rest_groups,  subset= meanFD_rest_groups$group == "control")
df_control_wide <- reshape(df_control, v.names = 'meanFD', idvar = 'NECOS_ID', timevar = 'rs_block', direction = 'wide')

df_control_wide$demeaned_rest_gr1 <- df_control_wide$meanFD.block_1 - mean(meanFD_rest_groups_wide$meanFD.block_1)
df_control_wide$demeaned_rest_gr2 <- df_control_wide$meanFD.block_2 - mean(meanFD_rest_groups_wide$meanFD.block_2)
df_control_wide$demeaned_rest_gr3 <- df_control_wide$meanFD.block_3 - mean(meanFD_rest_groups_wide$meanFD.block_3)

df_control_wide$delta_rest_gr2_gr1 <-  df_control_wide$meanFD.block_2 - df_control_wide$meanFD.block_1 
df_control_wide$demeaned_del_rest_gr2_gr1 <- df_control_wide$delta_rest_gr2_gr1 - mean(meanFD_rest_groups_wide$delta_rest_gr2_gr1)


t.test(df_control_wide$demeaned_del_rest_gr2_gr1,
       df_stress_wide$demeaned_del_rest_gr2_gr1)

write.csv(df_stress_wide, file= '/scr/nil3/reinelt/NECOS/2nd_level_analysis/meanFD_stress_group.csv')
write.csv(df_control_wide, file= '/scr/nil3/reinelt/NECOS/2nd_level_analysis/meanFD_control_group.csv')

#df_stress   <- subset(meanFD_rest_groups, select = c("NECOS_ID","rs_block", "group", "meanFD"), subset= meanFD_rest_groups$group == "stress")
#colnames(df_stress)[4] <- "par"

df_control  <- subset(meanFD_rest_groups, select = c("NECOS_ID","rs_block", "group", "meanFD"), subset= meanFD_rest_groups$group == "control")
#colnames(df_control)[4] <- "par"

df_stress_wide <- reshape(df_stress, v.names = 'meanFD', idvar = 'NECOS_ID', timevar = 'rs_block', direction = 'wide')

