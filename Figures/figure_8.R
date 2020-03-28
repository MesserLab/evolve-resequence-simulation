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
  for (i in c(11)){
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
                   , c("Equal", "Exponential")
                   , c("Equal (Standard)", "Exponential")
                   , c("Mutant is codominant", "Mutant is recessive", "Mutant is dominant")
                   , c("Mutant is codominant", "Mutant is recessive", "Mutant is dominant")
                   , c("Random (Standard)", "Clustered")
                   , c("Random (Standard)", "Clustered")
                   , c("No epistatisis", "Synergistic, weak", "Synergistic, strong", "Antaganistic, weak", "Antaganistic, strong", "Sign, weak",  "Sign, strong", "Reciprocal sign, weak", "Reciprocal sign, strong")
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
        scale_y_continuous(limits = c(0,0.61)) +
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
        theme(legend.position = c(0.7, 0.24)) +
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

## B. Heritability ################
library(ggplot2)
for (i in 11){
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
  
  filenames <-   list(c("NQTL10_NGen5", "NQTL10_FewerSNPs_NGen5")
                      , c("NQTL10_NGen5", "NQTL100_NGen5", "NQTL2_NGen5", "NQTL20_NGen5", "NQTL200_NGen5")
                      , c("NQTL10_NGen5", "NQTL10_FreqLow5_NGen5", "NQTL10_FreqHigh5_NGen5")
                      , c("NQTL100_NGen5", "NQTL100_FreqLow5_NGen5", "NQTL100_FreqHigh5_NGen5")
                      , c("NQTL10_NGen5", "NQTL10_ESDistE_NGen5")
                      , c("NQTL100_NGen5", "NQTL100_ESDistE_NGen5")
                      , c("NQTL10_NGen5",  "NQTL10_D0_NGen5", "NQTL10_D1_NGen5")
                      , c("NQTL100_NGen5", "NQTL100_D0_NGen5", "NQTL100_D1_NGen5")
                      , c("NQTL10_NGen5", "NQTL10_Clustered_NGen5")
                      , c("NQTL100_NGen5", "NQTL100_Clustered_NGen5")
                      , c("NQTL10", "NQTL10_EpiSce9", "NQTL10_EpiSce4", "NQTL10_EpiSce8", "NQTL10_EpiSce1", "NQTL10_EpiSce5", "NQTL10_EpiSce10", "NQTL10_EpiSce7", "NQTL10_EpiSce11")
                      , c("NQTL100", "NQTL100_ESDistE", "NQTL100_NGen20SelectedSize200", "NQTL100_NGen20SelectedSize200ESDistE")
                      , c("NQTL100", "NQTL100_NGen5")
                      , c("NQTL10", "NQTL10_NGen5")
  )
  
  labels <- list(c("~14000 SNPs \n(Standard)", "~1400 SNPs")
                 , c("10 QTL (Standard)", "100 QTL (Standard)", "2 QTL", "20 QTL",  "200 QTL")
                 , c("Random (Standard)", "Lower than 5%", "Higher than 5%")
                 , c("Random (Standard)", "Lower than 5%", "Higher than 5%")
                 , c("Equal (Standard)", "Exponential")
                 , c("Equal (Standard)", "Exponential")
                 , c("Mutant is codominant \n(Standard)", "Mutant is recessive", "Mutant is dominant")
                 , c("Mutant is codominant \n(Standard)", "Mutant is recessive", "Mutant is dominant")
                 , c("Random (Standard)", "Clustered")
                 , c("Random (Standard)", "Clustered")
                 , c("No epistatisis", "Synergistic, weak", "Synergistic, strong", "Antagonistic, weak", "Antagonistic, strong", "Sign, weak",  "Sign, strong", "Reciprocal sign, weak", "Reciprocal sign, strong")
                 , c("Fixed, 10 Generations, 10% Selected (Basline)", "Exponential, 10 Generations, 10% Selected", "Fixed, 20 Generations, 20% Selected", "Exponential, 20 Generations, 20% Selected")
                 , c("10 Generations (Standard)", "5 Generations")
                 , c("10 Generations (Standard)", "5 Generations")
  )
  title <- titles[i]
  filenames <- filenames[[i]]
  labels <- labels[[i]]
  n_files <- length(filenames)
  
  for (filename in filenames){
    print(filename)
    pooled_trait_mean <- matrix(NA, ncol=10, nrow=200)
    pooled_selected_trait_mean <- matrix(NA, ncol=10, nrow=200)
    pooled_trait_variance <- matrix(NA, ncol=10, nrow=200)
    pooled_h_squared <- matrix(NA, ncol=10, nrow=200)
    pooled_number_of_mutant_alleles <- matrix(NA, ncol=10, nrow=200)
    row_number<-1
    for (k in 1:100){
      for (Direction in c("Plus", "Minus")) {
        sign<-Direction=="Minus"
        setwd(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Simulations/", filename, "/SimRep", k, "/ExpRep", Direction, "1/"))
        full_trait <- matrix(NA,ncol = 11, nrow = 1000)
        selected_trait_mean <- vector("numeric", length = 10)
        for (i in 1:11){
          assign(paste0("gen_", i,"_trait"), read.table(paste0("Gen", i, "_Trait.txt"), fill=T))
          full_trait[,i] <- get(paste0("gen_", i,"_trait"))[,1]
        }  
        
        ## get summary stats from all individuals per generation
        trait_mean <- apply(full_trait, 2, mean)
        full_trait_sorted <- apply(full_trait, 2, function(x) {sort(x, decreasing=sign)})
        selected_trait_mean <- apply(full_trait_sorted[901:1000,], 2, mean)
        trait_variance <- apply(full_trait, 2, var)
        
        ## get heritability from each generation
        h_squared <- vector("numeric", length = 10)
        for (j in 1:10) {
          h_squared[j] <- (trait_mean[j+1]-trait_mean[j])/(selected_trait_mean[j]-trait_mean[j])
        }
        ## add this replicate to pooled data
        pooled_trait_mean[row_number,] <- trait_mean[1:10]
        pooled_selected_trait_mean[row_number,] <- selected_trait_mean[1:10]
        pooled_trait_variance[row_number,] <- trait_variance[1:10]
        pooled_h_squared[row_number,] <- h_squared
        
        ## number of QTLs over time
        for (l in 1:10){
          assign(paste0("gen_", l,"_number_of_mutant_alleles"), read.table(paste0("Gen", l, "_NSegregatingQTL.txt"), fill=T))
          pooled_number_of_mutant_alleles[row_number,l] <- get(paste0("gen_", l,"_number_of_mutant_alleles"))[1,1]
        }  
        row_number <- row_number+1
      }
    }
    
    #########################################################################################
    
    pooled_h_squared[which(pooled_trait_variance==0)] <- 1
    pooled_additive_variance <- pooled_h_squared*pooled_trait_variance
    pooled_additive_variance[which(pooled_trait_variance==0)] <- 0
    
    mean_h_squared <- apply(pooled_h_squared, 2, function(x){mean(x, na.rm = T)})
    mean_trait_variance <- apply(pooled_trait_variance, 2, function(x){mean(x, na.rm = T)})
    mean_additive_variance <- apply(pooled_additive_variance, 2, function(x){mean(x, na.rm = T)})
    mean_trait_value <- apply(pooled_trait_mean, 2, function(x){mean(x, na.rm = T)})
    
    
    mean_number_of_mutant_alleles <- apply(pooled_number_of_mutant_alleles, 2, mean)
    color <- c(1,2,3,4,5,6,7,8,9,10)[which(filenames == filename)]
    time_points <- c(1:10)
    mean_number_of_lost_mutant_alleles <- max(mean_number_of_mutant_alleles)-mean_number_of_mutant_alleles
    
    gg_frame <- data.frame(time_points, mean_trait_variance, mean_additive_variance, mean_number_of_mutant_alleles, mean_number_of_lost_mutant_alleles, mean_h_squared, filename, color, mean_trait_value)
    if (filename == filenames[1]){
      gg_frame_all <- gg_frame
    } else {
      gg_frame_all <- rbind(gg_frame_all, gg_frame)
    }
  }
  
  
  p1 <- ggplot(gg_frame_all, aes(time_points,mean_trait_variance, color = factor(color, labels = labels))) + 
    geom_line(size=1.2) + 
    geom_point(size = 3, alpha = 0.8) +
    xlim(1, 5) + 
    ylim(0, max(gg_frame_all$mean_trait_variance)+0.1) + 
    #theme(legend.position="none") +
    theme(axis.title = element_text(size = 30)) +
    theme(axis.text=element_text(size=30)) +
    ggtitle(NULL)
  
  p2 <- ggplot(gg_frame_all, aes(time_points,mean_additive_variance, color = factor(color, labels = labels))) + 
    geom_line(size=1.2) + 
    geom_point(size = 3, alpha = 0.8) +
    xlim(1, 5) + 
    ylim(0, max(gg_frame_all$mean_trait_variance)+0.1) + 
    #theme(legend.position="none") +
    theme(axis.title = element_text(size = 30)) +
    theme(axis.text=element_text(size=30)) +
    ggtitle(NULL)
  
  p3 <- ggplot(gg_frame_all, aes(time_points,mean_h_squared, color = factor(color, labels = labels))) + 
    geom_line(size=1.8) + 
    geom_point(size = 5) +
    theme_bw() +
    ggtitle(NULL) +
    scale_x_continuous(limits = c(1,5), breaks = seq(1,5,1)) + 
    scale_y_continuous(limits = c(0.65,1.02), breaks = seq(0.7,1,0.1)) + 
    xlab("Generation number") +
    ylab("Mean heritability") +
    theme(legend.position="none") +
    theme(axis.title = element_text(size = 30)) +
    theme(axis.text=element_text(size=30)) +
    theme(text = element_text(size=30)) +
    ggtitle(NULL) +
    theme(legend.position = c(0.7, 0.24)) +
    theme(legend.text=element_text(size=30)) +
    theme(legend.background = element_rect(fill = "transparent", color="black")) +
    theme(legend.title = element_blank()) +
    theme(axis.title = element_text(size = 30)) +
    theme(axis.text=element_text(size=30))
  
  p4 <- ggplot(gg_frame_all, aes(time_points, mean_number_of_lost_mutant_alleles, color = factor(color, labels = labels))) + 
    geom_line(size=1.2) + 
    geom_point(size = 3, alpha = 0.8) +
    ## xlim(0, 10) + 
    ylim(0, max(gg_frame_all$mean_number_of_mutant_alleles)+1) + 
    guides(color=guide_legend(title=title)) +
    theme(legend.position = c(0.65, 0.2)) +
    theme(axis.title = element_text(size = 30)) +
    theme(axis.text=element_text(size=30)) +
    theme(legend.text=element_text(size=25)) +
    theme(legend.background = element_rect(fill = "transparent")) +
    theme(text = element_text(size=30)) +
    ggtitle(NULL)
  
  p5 <- ggplot(gg_frame_all, aes(time_points, mean_trait_value, color = factor(color, labels = labels))) + 
    geom_line(size=1.2) + 
    geom_point(size = 3, alpha = 0.8) +
    ## xlim(0, 0.025) + 
    ylim(min(gg_frame_all&mean_trait_value)-1, max(mean_trait_value)+1) + 
    theme(legend.position="none") +
    theme(axis.title = element_text(size = 30)) +
    theme(axis.text=element_text(size=30)) +
    ggtitle(NULL)
  
  
  
  setwd(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Figures/Misc/"))
  png(paste0(title[1], "_heritability.png"), width = 800, height = 750, units = "px", pointsize = 20)
  print(p3)
  dev.off()
}

panel_b <- p3

## Assemble these into one figure #########################################

library(cowplot)
figure_8 <- plot_grid(panel_a, panel_b, labels=c("A", "B"), nrow = 1, label_size=30)

png(paste0("~/evolve-resequence-simulation/Figures/figure_8.png"), width = 800*2, height = 750, units = "px", pointsize = 20)
print(figure_8)
dev.off()

