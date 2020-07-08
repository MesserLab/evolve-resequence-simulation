## This script is used to create ROC curve comparisons among different trait architecture or experimental design scenarios.
## First define the type of comparisons that you want to make. Then, under the following if statement, define the appropriate variables for type of your comparison. Also define the variables outside of the if statement. Startin from line 75, also find the appropriate if statement and follow instructions there. 
## This step is usually very quick, so an interactive R session will be more suitable.
## You might find it easier to write your own R code for your specific purpose to plot the ROC curve. 


rm(list=ls(all=TRUE))

Type <- "ComputationalMethod" # can be "QTL_architecture", "WindowSize", "ComputationalMethod", "Replication", "Direction", "Experimental_design"

if (Type == "Direction"){
  NExpRep <- c(1)
  Directions <- c("Plus", "Minus")
} else if (Type == "Replication"){
  NsExpRep <- c(1,2,5)
  Direction <- c("Minus")
} else {
  NExpRep <- 1
  Directions <- c("Minus"
  )
}
WindowSize <- 0
# ComputationalMethod <- "D" # comment out if doing method comparison
filename <- "NQTL10_FewerSNPs" # When Type=="Replication", "Direction", "ComputationalMethod"
title <- paste0(Type, "_", filename) # When Type=="Replication", "Direction", "ComputationalMethod"
OutPath="/fs/cbsubscb10/storage/rl683/TemporalScan/" 
SampleSize=50 # make sure that the simulation had SampleSize LARGER or EQUAL to this
NSimRep=100 # make sure that the simulation had NSimRep LARGER or EQUAL to this

library(tidyverse)

if (Type %in% c("ComputationalMethod")){
  for (Direction in Directions){
    for (sim_rep_id in 1:NSimRep){
      setwd(paste0(OutPath, "Simulations/", filename, "/SimRep", sim_rep_id, "/", "ExpRep", Direction, NExpRep, "/")) ## set directory to experimental replication
      ROC <- read.table("ROC_TwoTimepoints_NGen5.txt", sep = ",",  header = T, stringsAsFactors = F)
      temp <- read.table("ROC_AllTimepoints_NGen5.txt", sep = ",",  header = T, stringsAsFactors = F)
      temp2 <- read.table("ROC_AllTimepoints_NGen5_TransformedD.txt", sep = ",",  header = T, stringsAsFactors = F) %>%
        mutate(Computational_method = "TransformedD")
      ROC <- rbind(ROC, temp, temp2)
      ROC[,5]<-paste0(ROC[,5], "_", ROC[,6])
      ROC <- ROC[,-6]
      colnames(ROC)<- c("Proportion_of_genetic_variance_in_gen_1", 
                        "Proportion_of_QTL_detected_weighted", 
                        "False_positive_rate",
                        "Window_size",
                        "Computational_method")
      if(sim_rep_id==1){
        ROC_temp <- ROC[,1:3]
      } else {
        ROC_temp <- ((sim_rep_id-1)*ROC_temp + ROC[,1:3])/sim_rep_id
      }
    }
    ROC_temp<-ROC_temp[c(1:100, 501:600, 101:200, 301:400, 201:300, 401:500),]
    ROC<-ROC[c(1:100, 501:600, 101:200, 301:400, 201:300, 401:500),]
    color <- sapply(ROC[,5], function(x) {which(unique(ROC[,5])==x)})
    labels <- c("D (2 Timepoints)",
                "Transformed D (2 Timepoints)", 
                "WFABC (2 Timepoints)", 
                "WFABC (5 Timepoints)",
                "ApproxWF (2 Timepoints)",
                "ApproxWF (5 Timepoints)")
    Computational_method <-ROC[,5]
    ROC_final <- cbind(ROC_temp, Computational_method, color)
    rm(Computational_method)
    
    
    p1 <- ggplot(ROC_final, aes(False_positive_rate,Proportion_of_genetic_variance_in_gen_1, color = factor(color, labels = labels))) + 
      geom_line(size = 2, alpha = 0.8, linetype="dotted") + 
      geom_line(data=ROC_final[which(ROC_final[,2]>0 | ROC_final[,3]>0),], aes(False_positive_rate, Proportion_of_genetic_variance_in_gen_1, color = factor(color, labels = labels)), size = 2, alpha = 1) + 
      scale_x_continuous(breaks = seq(0,0.1,0.02)) + 
      scale_y_continuous(limits = c(0,1)) + 
      coord_cartesian(xlim = c(0,0.1)) +
      xlab("False positive rate") +
      ylab("Proportion of genetic variance detected") +
      theme_bw() +
      theme(legend.position="none") +
      theme(axis.title = element_text(size = 30)) +
      theme(axis.text=element_text(size=30)) +
      theme(text = element_text(size=30)) +
      ggtitle(NULL)
    
    
    p2 <- ggplot(ROC_final, aes(False_positive_rate,Proportion_of_QTL_detected_weighted, color = factor(color, labels = labels))) + 
      geom_line(size = 2, alpha = 0.8, linetype="dotted") + 
      geom_line(data=ROC_final[which(ROC_final[,2]>0 | ROC_final[,3]>0),], aes(False_positive_rate, Proportion_of_QTL_detected_weighted, color = factor(color, labels = labels)), size = 2, alpha = 1) + 
      scale_x_continuous(breaks = seq(0,0.1,0.02)) +
      scale_y_continuous(limits = c(0,1)) + 
      coord_cartesian(xlim = c(0,0.1)) +
      xlab("False positive rate") +
      ylab("Proportion of QTLs detected") +
      theme_bw() +
      theme(legend.position="none") +
      theme(axis.title = element_text(size = 30)) +
      theme(axis.text=element_text(size=30)) +
      theme(text = element_text(size=30)) +
      ggtitle(NULL)
    
    p3 <- p1 +
      theme(legend.position = c(0.7, 0.2)) +
      theme(legend.text=element_text(size=30)) +
      theme(legend.background = element_rect(fill = "transparent", color="black")) +
      theme(legend.title = element_blank()) +
      theme(axis.title = element_text(size = 30)) +
      theme(axis.text=element_text(size=30))
    
    p4 <- p2 +
      theme(legend.position = c(0.65, 0.83)) +
      theme(legend.text=element_text(size=30)) +
      theme(legend.background = element_rect(fill = "transparent", color="black")) +
      theme(legend.title = element_blank()) +
      theme(axis.title = element_text(size = 30)) +
      theme(axis.text=element_text(size=30))
    
    setwd(paste0(OutPath, "Figures/ROC/", Type, "/"))
    
    png(paste0(title[1], "_", Direction, "_Zoomed_Variance.png"), width = 1500, height = 760, units = "px", pointsize = 20)
    print(p1)
    dev.off()
    
    png(paste0(title[1], "_", Direction, "_Zoomed_Unweighted.png"), width = 1500, height = 760, units = "px", pointsize = 20)
    print(p2)
    dev.off()
    
    png(paste0(title[1], "_", Direction, "_Zoomed_Variance_Legend.png"), width = 800, height = 750, units = "px", pointsize = 20)
    print(p3)
    dev.off()
    
    png(paste0(title[1], "_", Direction, "_Zoomed_Unweighted_Legend.png"), width = 800, height = 750, units = "px", pointsize = 20)
    print(p4)
    dev.off()
    
    png(paste0("~/evolve-resequence-simulation/Figures/figure_2.png"), width = 800, height = 750, units = "px", pointsize = 20)
    print(p4)
    dev.off()

    jpeg(paste0("~/evolve-resequence-simulation/Figures/figure_2.jpeg"), width = 800, height = 750, units = "px", pointsize = 20)
    print(p4)
    dev.off()
    
  }
}

