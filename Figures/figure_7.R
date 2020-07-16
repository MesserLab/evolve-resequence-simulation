## This script is used to create ROC curve comparisons among different trait architecture or experimental design scenarios.
## First define the type of comparisons that you want to make. Then, under the following if statement, define the appropriate variables for type of your comparison. Also define the variables outside of the if statement. Startin from line 75, also find the appropriate if statement and follow instructions there. 
## This step is usually very quick, so an interactive R session will be more suitable.
## You might find it easier to write your own R code for your specific purpose to plot the ROC curve. 

## A and B. ROC plots ##########################
rm(list=ls(all=TRUE))

Type <- "QTL_architecture" # can be "QTL_architecture", "WindowSize", "ComputationalMethod", "Replication", "Direction", "Experimental_design"

if (Type == "Direction"){
  NExpRep <- c(1)
  Directions <- c("Plus", "Minus")
} else if (Type == "Replication"){
  NsExpRep <- c(1,2,5)
  Direction <- c("Minus")
} else {
  NExpRep <- 1
  Directions <- c("Plus"
                  , "Minus"
  )
}
WindowSize <- 0
ComputationalMethod <- "D" # comment out if doing method comparison
filename <- "NQTL10_FewerSNPs" # When Type=="Replication", "Direction", "ComputationalMethod"
title <- paste0(Type, "_", filename) # When Type=="Replication", "Direction", "ComputationalMethod"
OutPath="/fs/cbsubscb10/storage/rl683/TemporalScan/" 
SampleSize=50 # make sure that the simulation had SampleSize LARGER or EQUAL to this
NSimRep=100 # make sure that the simulation had NSimRep LARGER or EQUAL to this

library(ggplot2)

if (Type %in% c("QTL_architecture", "Experimental_design")){ # define filenames, plot title, and plot labels in the begining of the following loop 
  for (i in c(7)){
    titles <- c("Number_of_SNPs"
                , "Number_of_QTLs"
                , "Starting_frequency_10_QTLs"
                , "Starting_frequency_100_QTLs"
                , "Effect_size_distribution_10_QTLs"
                , "Effect_size_distribution_100_QTLs"
                , "Dominance_10_QTLs"
                , "Dominance_100_QTLs"
                , "Position_10_QTLs"
                , "Position_100_QTLs"
                , "Epistasis_scenarios"
                , "ESDistxNGenxSelectedSize"
                , "Number_of_generations_100_QTLs"
                , "Number_of_generations_10_QTLs"
    )
    
    filenames <-   list(c("NQTL10", "NQTL10_FewerSNPs")
                        , c("NQTL2", "NQTL10", "NQTL20", "NQTL100", "NQTL200")
                        , c("NQTL10", "NQTL10_FreqLow5", "NQTL10_FreqHigh5")
                        , c("NQTL100", "NQTL100_FreqLow5", "NQTL100_FreqHigh5")
                        , c("NQTL10", "NQTL10_ESDistE")
                        , c("NQTL100", "NQTL100_ESDistE")
                        , c("NQTL10", "NQTL10_D1",  "NQTL10_D0")
                        , c("NQTL100", "NQTL100_D0", "NQTL100_D1")
                        , c("NQTL10", "NQTL10_Clustered")
                        , c("NQTL100", "NQTL100_Clustered")
                        , c("NQTL10", "NQTL10_EpiSce9", "NQTL10_EpiSce4", "NQTL10_EpiSce8", "NQTL10_EpiSce1", "NQTL10_EpiSce5", "NQTL10_EpiSce10", "NQTL10_EpiSce7", "NQTL10_EpiSce11")
                        , c("NQTL100", "NQTL100_ESDistE", "NQTL100_NGen20SelectedSize200", "NQTL100_NGen20SelectedSize200ESDistE")
                        , c("NQTL100", "NQTL100")
                        , c("NQTL10", "NQTL10")
    )
    
    labels <- list(c("~14000 SNPs \n(Standard)", "~1400 SNPs")
                   , c("2 QTLs", "10 QTLs", "20 QTLs", "100 QTLs",  "200 QTLs")
                   , c("Random (Standard)", "Lower than 5%", "Higher than 5%")
                   , c("Random (Standard)", "Lower than 5%", "Higher than 5%")
                   , c("Equal", "Exponential")
                   , c("Equal (Standard)", "Exponential")
                   , c("Additive", "Dominant", "Recessive")
                   , c("Mutant is codominant", "Mutant is recessive", "Mutant is dominant")
                   , c("Random (Standard)", "Clustered")
                   , c("Random (Standard)", "Clustered")
                   , c("No epistatisis (Standard)", "Synergistic, weak", "Synergistic, strong", "Antaganistic, weak", "Antaganistic, strong", "Sign, weak",  "Sign, strong", "Reciprocal sign, weak", "Reciprocal sign, strong")
                   , c("Fixed, 10 Generations, 10% Selected (Basline)", "Exponential, 10 Generations, 10% Selected", "Fixed, 20 Generations, 20% Selected", "Exponential, 20 Generations, 20% Selected")
                   , c("10 Generations (Standard)", "5 Generations")
                   , c("10 Generations (Standard)", "5 Generations")
    )
    
    title <- titles[i]
    filenames <- filenames[[i]]
    labels <- labels[[i]]
    
    
    n_files <- length(filenames)
    for (Direction in Directions) {
      for (filename in filenames){
        for (k in 1:NSimRep){
          if (Direction == "Plus") {
            DirectionName <- "SingleDirection"
          } else if (Direction == "Minus"){
            DirectionName <- "OppositeDirections"
          }
          setwd(paste0(OutPath, "Simulations/", filename, "/SimRep", k, "/ExpRep", Direction, NExpRep))
          ROC <- read.table("ROC_AllTimepoints_NGen5_TransformedD.txt", sep = ",",  header = T, stringsAsFactors = F)
          if(dim(ROC)[2]==6){
            ROC <- cbind(ROC, ROC[,2])
          }
          colnames(ROC)<- c("Proportion_of_genetic_variance_in_gen_1", 
                            "Proportion_of_QTL_detected_weighted", 
                            "False_positive_rate",
                            "Window_size",
                            "Computational_method",
                            "Numbers_of_timepoints_sampled",
                            "Proportion_of_QTL_detected_unweighted")
          ROC <- ROC[which(ROC$Window_size==WindowSize & ROC$Computational_method==ComputationalMethod),]
          if (k == 1){
            ROC_final <- ROC[,c(1:3,7)]
          } else {
            ROC_final <- ROC[,c(1:3,7)]+ROC_final
          }
        }
        ROC_final <- ROC_final/NSimRep
        ROC_final <- cbind(ROC_final, filename)
        color <- which(filenames == filename)
        ROC_final <- cbind(ROC_final, color) 
        if (!exists("ROC_QTL_archetecture_pooled")){
          ROC_QTL_archetecture_pooled <- ROC_final
        } else {
          ROC_QTL_archetecture_pooled <- rbind(ROC_QTL_archetecture_pooled, ROC_final)
        } 
      }
      library(ggplot2)
      
      p1 <- ggplot(ROC_QTL_archetecture_pooled, aes(False_positive_rate, Proportion_of_genetic_variance_in_gen_1, color = factor(color, labels = labels))) + 
        geom_line(size = 2, alpha = 0.8, linetype="dotted") + 
        geom_line(data=ROC_QTL_archetecture_pooled[which(ROC_QTL_archetecture_pooled[,1]>0 | ROC_QTL_archetecture_pooled[,3]>0),], aes(False_positive_rate, Proportion_of_genetic_variance_in_gen_1, color = factor(color, labels = labels)), size = 2, alpha = 1) + 
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
      p2 <- ggplot(ROC_QTL_archetecture_pooled, aes(False_positive_rate,Proportion_of_QTL_detected_weighted, color = factor(color, labels = labels))) + 
        geom_line(size = 2, alpha = 0.8, linetype="dotted") + 
        geom_line(data=ROC_QTL_archetecture_pooled[which(ROC_QTL_archetecture_pooled[,2]>0 | ROC_QTL_archetecture_pooled[,3]>0),], aes(False_positive_rate, Proportion_of_QTL_detected_weighted, color = factor(color, labels = labels)), size = 2, alpha = 1) + 
        scale_x_continuous(breaks = seq(0,0.1,0.02)) + 
        scale_y_continuous(limits = c(0,1)) + 
        coord_cartesian(xlim = c(0,0.1)) +
        xlab("False positive rate") +
        ylab("Proportion of QTLs detected weighted") +
        theme_bw() +
        theme(legend.position="none") +
        theme(axis.title = element_text(size = 30)) +
        theme(axis.text=element_text(size=30)) +
        theme(text = element_text(size=30)) +
        ggtitle(NULL) 
      p3 <- ggplot(ROC_QTL_archetecture_pooled, aes(False_positive_rate,Proportion_of_QTL_detected_unweighted, color = factor(color, labels = labels))) +
        geom_line(size = 2, alpha = 0.8, linetype="dotted") +
        geom_line(data=ROC_QTL_archetecture_pooled[which(ROC_QTL_archetecture_pooled[,2]>0 | ROC_QTL_archetecture_pooled[,3]>0),], aes(False_positive_rate, Proportion_of_QTL_detected_unweighted, color = factor(color, labels = labels)), size = 2, alpha = 1) +
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
      
      p4 <- p1 +
        theme(legend.position = c(0.75, 0.2)) +
        theme(legend.text=element_text(size=30)) +
        theme(legend.background = element_rect(fill = "transparent", color="black")) +
        theme(legend.title = element_blank()) +
        theme(axis.title = element_text(size = 30)) +
        theme(axis.text=element_text(size=30))
      p5 <- p2 + 
        theme(legend.position = c(0.75, 0.2)) +
        theme(legend.text=element_text(size=30)) +
        theme(legend.background = element_rect(fill = "transparent", color="black")) +
        theme(legend.title = element_blank()) +
        theme(axis.title = element_text(size = 30)) +
        theme(axis.text=element_text(size=30))
      p6 <- p3 +
        theme(legend.position = c(0.75, 0.2)) +
        theme(legend.text=element_text(size=30)) +
        theme(legend.background = element_rect(fill = "transparent", color="black")) +
        theme(legend.title = element_blank()) +
        theme(axis.title = element_text(size = 30)) +
        theme(axis.text=element_text(size=30))
      
      setwd(paste0(OutPath, "Figures/ROC/", title, "/"))
      
      png(paste0(title[1], "_", DirectionName, "_Zoomed_NGen5_Variance.png"), width = 800, height = 750, units = "px", pointsize = 20)
      print(p1)
      dev.off()
      png(paste0(title[1], "_", DirectionName, "_Zoomed_NGen5_Weighted.png"), width = 800, height = 750, units = "px", pointsize = 20)
      print(p2)
      dev.off()
      png(paste0(title[1], "_", DirectionName, "_Zoomed_NGen5_Unweighted.png"), width = 800, height = 750, units = "px", pointsize = 20)
      print(p3)
      dev.off()
      png(paste0(title[1], "_", DirectionName, "_Zoomed_NGen5_Variance_Legend.png"), width = 800, height = 750, units = "px", pointsize = 20)
      print(p4)
      dev.off()
      png(paste0(title[1], "_", DirectionName, "_Zoomed_NGen5_Weighted_Legend.png"), width = 800, height = 750, units = "px", pointsize = 20)
      print(p5)
      dev.off()
      png(paste0(title[1], "_", DirectionName, "_Zoomed_NGen5_Unweighted_Legend.png"), width = 800, height = 750, units = "px", pointsize = 20)
      print(p6)
      dev.off()
      
      rm(ROC_QTL_archetecture_pooled)
    }
  }
} 
panel_a <- p6
panel_b <- p4

## D. Plot D histogram (combining all simulation reps) ####################################################################
setwd("/fs/cbsubscb10/storage/rl683/TemporalScan/")

library(ggplot2)
library(tidyverse)
scenarios=c("NQTL10", "NQTL10_D0", "NQTL10_D1", "NQTL100", "NQTL10_Clustered", "NQTL100_Clustered")
for (j in c(1,2,3)) {
  print(j)
  scenario=scenarios[j]
  for (i in 1:100){
    setwd(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Simulations/", scenario, "/SimRep", i, "/ExpRepPlus1"))
    df1 <- read.table("SampleFreqAllTimepoints.txt", sep = ",", header = T)
    setwd(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Simulations/", scenario, "/SimRep", i, "/ExpRepMinus1"))
    df2 <- read.table("SampleFreqAllTimepoints.txt", sep = ",", header = T)
    df1 <- df1[,c(1, 3, 4, 11)]
    colnames(df1) <- c("PermanentID", "Position", "EffectSize", "FrequencyPlus")
    df2 <- df2[,c(1, 3, 4, 11, 7)]
    colnames(df2) <- c("PermanentID", "Position", "EffectSize", "FrequencyMinus", "StartingFrequency")
    df <- merge(df1, df2, by="PermanentID", all=T)
    #Position<-apply(df[,c(2,5)], 1, function(x){mean(x, na.rm = T)})
    EffectSize<-apply(df[,c(3,6)], 1, function(x){mean(x, na.rm = T)})
    df <- df[,c(1,4,7,8)]
    df <- cbind(df, EffectSize)
    df[is.na(df)]<-0
    df <- mutate(df, D=abs(2*asin(sqrt(FrequencyPlus/100))-2*asin(sqrt(FrequencyMinus/100)))/pi) 
    
    if (i==1){
      gg_frame <- df
    }
    else {
      gg_frame<-rbind(gg_frame, df)
    }
  }

  EffectType <- rep("Neutral", dim(gg_frame)[1])
  EffectType[which(gg_frame$EffectSize>0)] <- "QTL"
  EffectType[which(gg_frame$EffectSize<0)] <- "QTL"
  gg_frame <-cbind(gg_frame, EffectType)
  gg_frame <- mutate(gg_frame, D_binned = cut(D, breaks=seq(0,1,0.1), include.lowest=T))
  
  p4 <- ggplot(data=gg_frame, aes(x=as.numeric(D_binned), colour = EffectType, fill=EffectType, group = EffectType)) +
    geom_bar(aes(y=..prop..), position='dodge') +
    scale_color_manual(values=c("Black", "Black"), label=c(" Neutral loci", " QTLs")) +
    scale_fill_manual(values=c("Grey", "Black"), label=c(" Neutral loci", " QTLs")) +
    xlab("Transformed D") +
    ylab("Proportion") +
    scale_x_continuous(breaks = c(seq(0, 10, by=1)+0.5),
                       labels = c(seq(0, 1, by=0.1))) +
    coord_cartesian(ylim=c(0,0.7)) +
    theme_bw() +
    theme(axis.text = element_text(size=30), axis.title = element_text(size=30)) +
    theme(text = element_text(size=30)) +
    theme(legend.position=c(0.75, 0.87)) +
    theme(legend.text=element_text(size=30)) +
    theme(legend.background = element_rect(fill = "transparent", color="black")) +
    theme(legend.title = element_blank()) +
    theme(panel.grid.minor.x = element_blank())
  
  
  png(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Figures/Misc/D_hist_", scenario, "_all.png"), width = 800, height = 750, units = "px", pointsize = 20)
  print(p4)
  dev.off()
  
  if (j==1) {
    panel_c <- p4
  } else if (j==2) {
    panel_e <- p4
  } else if (j==3) {
    panel_d <- p4
  }
}

## Assemble these into one figure #########################################

library(cowplot)
figure_7_up <- plot_grid(panel_a, panel_b, labels=c("A", "B"), nrow = 1, label_size=30)
figure_7_down <- plot_grid(panel_c, panel_d, panel_e, labels=c("C", "D", "E"), nrow = 1, label_size=30)

figure_7<- plot_grid(figure_7_up, figure_7_down, nrow = 2, rel_heights = c(1, 2/3))

png(paste0("~/evolve-resequence-simulation/Figures/figure_7.png"), width = 800*2, height = 750*(1+2/3), units = "px", pointsize = 20)
print(figure_7)
dev.off()

jpeg(paste0("~/evolve-resequence-simulation/Figures/figure_7.jpeg"), width = 800*2, height = 750*(1+2/3), units = "px", pointsize = 20)
print(figure_7)
dev.off()

set_smaller_font <- function(x){
  x +
    theme(legend.text = element_text(size=20),
          axis.title = element_text(size=20), 
          axis.text = element_text(size=20))
}
panel_a_smaller_font <- set_smaller_font(panel_a) +
  theme(legend.position = c(0.85, 0.2))
panel_b_smaller_font <- set_smaller_font(panel_b) +
  theme(legend.position = c(0.85, 0.2))
panel_c_smaller_font <- set_smaller_font(panel_c)
panel_d_smaller_font <- set_smaller_font(panel_d)
panel_e_smaller_font <- set_smaller_font(panel_e)

figure_7_up_smaller_font <- plot_grid(panel_a_smaller_font, panel_b_smaller_font, labels=c("A", "B"), nrow = 1, label_size=30)
figure_7_down_smaller_font <- plot_grid(panel_c_smaller_font, panel_d_smaller_font, panel_e_smaller_font, labels=c("C", "D", "E"), nrow = 1, label_size=30)

figure_7_smaller_font<- plot_grid(figure_7_up_smaller_font, figure_7_down_smaller_font, nrow = 2, rel_heights = c(1, 2/3))

jpeg(paste0("~/evolve-resequence-simulation/Figures/figure_7_smaller_font.jpeg"), width = 800*2, height = 750*(1+2/3), units = "px", pointsize = 20)
print(figure_7_smaller_font)
dev.off()
