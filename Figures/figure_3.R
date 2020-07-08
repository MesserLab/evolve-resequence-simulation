## This script is used to create ROC curve comparisons among different trait architecture or experimental design scenarios.
## First define the type of comparisons that you want to make. Then, under the following if statement, define the appropriate variables for type of your comparison. Also define the variables outside of the if statement. Startin from line 75, also find the appropriate if statement and follow instructions there. 
## This step is usually very quick, so an interactive R session will be more suitable.
## You might find it easier to write your own R code for your specific purpose to plot the ROC curve. 


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
  for (i in c(2)){
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
                        , c("NQTL10",  "NQTL10_D0", "NQTL10_D1")
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
                   , c("Equal (Standard)", "Exponential")
                   , c("Equal (Standard)", "Exponential")
                   , c("Mutant is codominant \n(Standard)", "Mutant is recessive", "Mutant is dominant")
                   , c("Mutant is codominant \n(Standard)", "Mutant is recessive", "Mutant is dominant")
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
        ylab("Proportion of QTLs detected weighted by effect size") +
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
        theme(legend.position = c(0.85, 0.75)) +
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
      
      png(paste0("~/evolve-resequence-simulation/Figures/figure_3.png"), width = 800, height = 750, units = "px", pointsize = 20)
      print(p6)
      dev.off()

      jpeg(paste0("~/evolve-resequence-simulation/Figures/figure_3.jpeg"), width = 800, height = 750, units = "px", pointsize = 20)
      print(p6)
      dev.off()
      
      rm(ROC_QTL_archetecture_pooled)
    }
  }
} 
