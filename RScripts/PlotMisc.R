## This script is used to plot miscellaneous figures for the simulation. It now includes plots of mean trait trajectory, trait/genetic variance trajectory, allele frequency trajectory, and distribution of D-value on the chromosome. More plot will be added in the future. 
## This script is meant to run with interactive R sessions.

rm(list=ls(all=TRUE))
library(ggplot2)
library(reshape)

## Plot trait charateristics across all replicates in the standard model ##########################################################################################
## Prepare data
library(tidyverse)
library(ggplot2)
scenarios <-  c("NQTL10", "NQTL100")
for (scenario in scenarios) {
  pooled_trait_mean <- matrix(NA, ncol=10, nrow=200)
  pooled_trait_variance <- matrix(NA, ncol=10, nrow=200)
  for (k in 1:2) {
    Direction <- c("Plus", "Minus")[k]
    for (j in 1:100){
      setwd(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Simulations/", scenario, "/SimRep", j, "/ExpRep", Direction, "1"))
      full_trait <- matrix(NA,ncol = 10, nrow = 1000)
      selected_trait_mean <- vector("numeric", length = 10)
      for (i in 1:10){
        assign(paste0("Gen", i,"_Trait"), read.table(paste0("Gen", i, "_Trait.txt"), fill=T))
        full_trait[,i] <- get(paste0("Gen", i,"_Trait"))[,1]
      }  
      
      ## get summary stats from all individuals per generation
      trait_mean <- apply(full_trait, 2, mean)
      trait_variance <- apply(full_trait, 2, var)
      
      ## enter these into matrices
      pooled_trait_mean[j+(k-1)*100,] <- trait_mean
      pooled_trait_variance[j+(k-1)*100,] <- trait_variance
    }
  }
  
  ## Average trait
  pooled_trait_mean<-as.data.frame(pooled_trait_mean)
  colnames(pooled_trait_mean)<-1:10
  pooled_trait_mean<-cbind(pooled_trait_mean, ReplicateNumber=c(1:200))
  pooled_trait_mean<-cbind(pooled_trait_mean, Direction=c(rep("Plus",100), rep("Minus",100)))
  pooled_trait_mean <- melt(pooled_trait_mean, id.vars=c("ReplicateNumber", "Direction"), value.name="AverageTraitValue", variable.name="Timepoint")
  colnames(pooled_trait_mean)[3:4]<-c("Timepoint", "AverageTraitValue")
  pooled_trait_mean[,3] <- as.numeric(as.character(pooled_trait_mean[,3]))
  p1 <-ggplot(data=pooled_trait_mean[which(pooled_trait_mean[,3]<=5),]) +
    geom_line(aes(x=Timepoint-1, y=AverageTraitValue, group=ReplicateNumber, color=Direction), alpha=0.6) +
    theme_bw() +
    theme(legend.position=c(0.2, 0.9)) +
    xlab("Generation number") +
    ylab("Average phenotype value") +
    theme(axis.title = element_text(size = 30)) +
    theme(axis.text=element_text(size=30)) +
    theme(text = element_text(size=30)) + 
    scale_color_manual(values=c("darkorange2", "purple"), labels = c("\"Low\" line", "\"High\" line")) +
    theme(legend.text=element_text(size=30)) +
    theme(legend.background = element_rect(fill = "transparent")) +
    guides(colour = guide_legend(override.aes = list(size=2))) +
    theme(legend.title = element_blank())
  
  ## Trait Variance
  pooled_trait_variance<-as.data.frame(pooled_trait_variance)
  colnames(pooled_trait_variance)<-1:10
  pooled_trait_variance<-cbind(pooled_trait_variance, ReplicateNumber=c(1:200))
  pooled_trait_variance<-cbind(pooled_trait_variance, Direction=c(rep("Plus",100), rep("Minus",100)))
  pooled_trait_variance <- melt(pooled_trait_variance, id.vars=c("ReplicateNumber", "Direction"), value.name="TraitVariance", variable.name="Timepoint")
  colnames(pooled_trait_variance)[3:4]<-c("Timepoint", "TraitVariance")
  pooled_trait_variance <- mutate(pooled_trait_variance, Timepoint=as.numeric(as.character(Timepoint)))
  p2 <-ggplot(data=filter(pooled_trait_variance, Timepoint<=5)) +
    geom_line(aes(x=Timepoint-1, y=TraitVariance, group=ReplicateNumber, color=Direction), alpha=0.6) +
    theme_bw() +
    theme(legend.position=c(0.8, 0.9)) +
    xlab("Generation number") +
    ylab("Phenotypic/genetic variance") +
    theme(axis.title = element_text(size = 30)) +
    theme(axis.text=element_text(size=30)) +
    theme(text = element_text(size=30)) +
    scale_color_manual(values=c("darkorange2", "purple"), labels = c("\"Low\" line", "\"High\" line")) +
    theme(legend.text=element_text(size=30)) +
    theme(legend.background = element_rect(fill = "transparent")) +
    guides(colour = guide_legend(override.aes = list(size=2))) +
    theme(legend.title = element_blank())
  
  # Plotting
  png(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Figures/Misc/TraitMean_", scenario, ".png"), width = 800, height = 750, units = "px", pointsize = 20)
  print(p1)
  dev.off()
  
  png(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Figures/Misc/TraitVariance_", scenario, ".png"), width = 800, height = 750, units = "px", pointsize = 20)
  print(p2)
  dev.off()
}


## Plot phenotype distribution change in one single experiment #########################################################################################
scenarios <- c("NQTL100", "NQTL10")
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

for (scenario in scenarios) {
  for (j in 1:1){
    for (k in 1:2) {
      Direction <- c("Plus", "Minus")[k]
      setwd(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Simulations/", scenario, "/SimRep", j, "/ExpRep", Direction, "1"))
      for (i in 1:5){
        assign(paste0("Gen", i,"_Trait"), read.table(paste0("Gen", i, "_Trait.txt"), fill=T))
        if (i==1 & k==1 ) {
          full_trait <- cbind(get(paste0("Gen", i,"_Trait"))[,1], i, Direction)
        } else {
          full_trait <- rbind(full_trait, cbind(get(paste0("Gen", i,"_Trait"))[,1], i, Direction))
        }
      }  
    }
    
    full_trait<-data.frame(full_trait)
    full_trait[,1]<-as.numeric(as.character(full_trait[,1]))
    full_trait[,2]<-as.factor(as.integer(as.character(full_trait[,2])))
    colnames(full_trait)<-c("PhenotypeValue", "GenerationNumber", "DirectionofSelection")
    
    p1 <- ggplot(data=full_trait, aes(x=GenerationNumber, y=PhenotypeValue, fill=DirectionofSelection)) +
      #geom_violin(bw = 1.5, alpha=0.8) +
      geom_split_violin(alpha=0.8, adjust=2.8) +
      theme_bw() +
      xlab("Generation number") +
      ylab("Phenotype value") +
      theme(axis.title = element_text(size = 30)) +
      theme(axis.text=element_text(size=30)) +
      theme(text = element_text(size=30)) + 
      scale_x_discrete(labels=0:4) +
      scale_fill_manual(values=c("darkorange2", "purple"), labels = c(" \"Low\" line", " \"High\" line")) +
      theme(legend.position=c(0.2, 0.9)) +
      theme(legend.text=element_text(size=30)) +
      theme(legend.background = element_rect(fill = "transparent")) +
      theme(legend.title = element_blank())
    
    
    
    png(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Figures/Misc/Trait_dist_", scenario, "_", j, ".png"), width = 800, height = 750, units = "px", pointsize = 20)
    print(p1)
    dev.off()
  }
}


## Plot allele frequency trajectory in one single experiment #########################################################################################
library(ggplot2)
library(reshape)

scenarios <- c("NQTL100", "NQTL100_Clustered")
legend_positions <- list(c(0.27, 0.9), c(0.7, 0.9))
for (k in 1:2) {
  scenario <- scenarios[k]
  lengend_position <- legend_positions[[k]]
  for (i in 1){
    setwd(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Simulations/", scenario, "/SimRep", i, "/ExpRepPlus1"))
    df1 <- read.table("SampleFreqAllTimepoints.txt", sep = ",", header = T)
    setwd(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Simulations/", scenario, "/SimRep", i, "/ExpRepMinus1"))
    df2 <- read.table("SampleFreqAllTimepoints.txt", sep = ",", header = T)
    
    colnames(df1)[7:11]<-0:4
    df1<- df1[,c(1, 4, 7:11)]
    df1 <- melt(df1, id.vars=c("ID_perm", "es"), value.name="Frequency", variable.name="Timepoint")
    
    colnames(df2)[7:11]<--4:0
    df2[,7:11]<-df2[,11:7]
    df2<- df2[,c(1, 4, 7:11)]
    df2 <- melt(df2, id.vars=c("ID_perm", "es"), value.name="Frequency", variable.name="Timepoint")
    df <- rbind(df2, df1)
    df[,4]<-df[,4]/100
    colnames(df) <- c("PermanentID", "EffectSize", "Timepoint", "Frequency")
    
    p1 <- mutate(df, EffectSize=as.factor(EffectSize)) %>%
      ggplot()  +
      theme_bw() +
      geom_line(aes(x=Timepoint, y=Frequency, group=PermanentID, colour = EffectSize, alpha=EffectSize, size=EffectSize)) +
      scale_alpha_manual(values=c(1, 0.012, 1), labels = c(" -1 allele", " Neutral allele", " +1 allele")) +
      scale_size_manual(values=c(1, 0.3, 1), labels = c(" -1 allele", " Neutral allele", " +1 allele")) +
      scale_color_manual(values=c("Blue", "Black", "Red"), labels = c(" -1 allele", " Neutral allele", " +1 allele")) +
      guides(colour = guide_legend(override.aes = list(size=c(2,1,2), alpha=c(1, 0.5, 1)))) +
      annotate("label", label="\"High\" line", x = 6.1, y = 0.8, size=10, colour='black', alpha=0.7) +
      annotate("label", label="\"Low\" line", x = 3.9, y = 0.8, size=10, colour='black', alpha=0.7) +
      annotate("segment", x=5, xend=5, y=0, yend=1.5, size=2) +
      annotate("segment", x=5.3, xend=6.9, y=0.75, yend=0.75, size=2, arrow = arrow(length = unit(0.5, "cm")), color ="black") +
      annotate("segment",x=4.7, xend=3.1, y=0.75, yend=0.75, size=2, arrow = arrow(length = unit(0.5, "cm")), color ="black") +
      scale_x_discrete(breaks=c(-4:4), labels=c(seq(4,0,-1), seq(1,4,1)), name="Generation number") +
      scale_y_continuous(breaks=seq(0.00,1.00,0.2), labels=seq(0.00,1.00,0.2), name="Allele frequency") +
      coord_cartesian(ylim = c(0,1)) +
      theme(axis.text = element_text(size=30), axis.title = element_text(size=30)) +
      theme(text = element_text(size=30)) +
      theme(legend.position=lengend_position) +
      theme(legend.text=element_text(size=30)) +
      theme(legend.background = element_rect(fill = "transparent")) +
      theme(legend.title = element_blank())
    
    
    png(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Figures/Misc/FrequencyTrajectory_", scenario, "_", i, ".png"), width = 800, height = 750, units = "px", pointsize = 20)
    print(p1)
    dev.off()
  }
}

## Plot distribution of D on chromosome in one single experiment #######################################################################################

library(ggplot2)
library(reshape)
library(tidyverse)

scenarios <- c("NQTL100", "NQTL100_Clustered")
legend_positions <- list(c(0.1, 0.08), c(0.1, 0.08))
for (k in 1:2) {
  scenario <- scenarios[k]
  lengend_position <- legend_positions[[k]]
  for (i in 1){
    setwd(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Simulations/", scenario, "/SimRep", i, "/ExpRepPlus1"))
    df1 <- read.table("SampleFreqAllTimepoints.txt", sep = ",", header = T)
    setwd(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Simulations/", scenario, "/SimRep", i, "/ExpRepMinus1"))
    df2 <- read.table("SampleFreqAllTimepoints.txt", sep = ",", header = T)
    
    df1 <- df1[,c(1, 3, 4, 11)]
    colnames(df1) <- c("PermanentID", "Position", "EffectSize", "FrequencyPlus")
    df2 <- df2[,c(1, 3, 4, 11, 7)]
    colnames(df2) <- c("PermanentID", "Position", "EffectSize", "FrequencyMinus", "StartingFrequency")
    df <- merge(df1, df2, by="PermanentID", all=T)
    Position<-apply(df[,c(2,5)], 1, function(x){mean(x, na.rm = T)})
    EffectSize<-apply(df[,c(3,6)], 1, function(x){mean(x, na.rm = T)})
    df <- df[,c(1,4,7,8)]
    df <- cbind(df, EffectSize)
    df[is.na(df)]<-0
    df <- cbind(df, D=df$FrequencyPlus-df$FrequencyMinus) 
    
    p1 <- mutate(df, EffectSize=as.factor(EffectSize), Position=Position/10^6) %>%
      ggplot()  +
      theme_bw() +
      geom_point(aes(x=Position, y=D/100, group=PermanentID, colour = EffectSize, alpha=EffectSize, size=EffectSize)) +
      scale_alpha_manual(values=c(1, 0.3, 1), labels = c(" -1 allele", " Neutral allele", " +1 allele")) +
      scale_size_manual(values=c(5, 1, 5), labels = c(" -1 allele", " Neutral allele", " +1 allele")) +
      scale_color_manual(values=c("Blue", "Black", "Red"), labels = c(" -1 allele", " Neutral allele", " +1 allele")) +
      guides(colour = guide_legend(override.aes = list(size=c(2,1,2), alpha=c(1, 0.5, 1)))) +
      scale_x_continuous(name="Position on chromosome (Mbp)") +
      scale_y_continuous(breaks=seq(-1.00,1.00,0.2), labels=seq(-1.00,1.00,0.2), name="D-value with sign") +
      coord_cartesian(ylim = c(-1,1)) +
      theme(axis.text = element_text(size=30), axis.title = element_text(size=30)) +
      theme(text = element_text(size=30)) +
      theme(legend.position=lengend_position) +
      theme(legend.text=element_text(size=30)) +
      theme(legend.background = element_rect(fill = "transparent")) +
      theme(legend.title = element_blank())
    
    
    png(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Figures/Misc/D_vs_position_", scenario, "_", i, ".png"), width = 1600, height = 750, units = "px", pointsize = 20)
    print(p1)
    dev.off()
  }
}


## Plot D ~ Starting Frequency (combining all simulation reps) ####################################################################
library(ggplot2)

scenarios <- c("NQTL100", "NQTL100_D0", "NQTL100_D1", "NQTL10", "NQTL10_D0", "NQTL10_D1", "NQTL10_FreqHigh5", "NQTL10_FreqLow5", "NQTL100_FreqHigh5", "NQTL100_FreqLow5")

for (k in 7:10) {
  scenario <- scenarios[k]
  for (i in 1:100){
    print(i)
    setwd(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Simulations/", scenario, "/SimRep", i, "/ExpRepPlus1"))
    df1 <- read.table("SampleFreqAllTimepoints.txt", sep = ",", header = T)
    setwd(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Simulations/", scenario, "/SimRep", i, "/ExpRepMinus1"))
    df2 <- read.table("SampleFreqAllTimepoints.txt", sep = ",", header = T)
    df1 <- df1[,c(1, 3, 4, 16)]
    colnames(df1) <- c("PermanentID", "Position", "EffectSize", "FrequencyPlus")
    df1 <- df1[which(df1$EffectSize!=0),]
    df2 <- df2[,c(1, 3, 4, 16, 7)]
    colnames(df2) <- c("PermanentID", "Position", "EffectSize", "FrequencyMinus", "StartingFrequency")
    df2 <- df2[which(df2$EffectSize!=0),]
    df <- merge(df1, df2, by="PermanentID", all=T)
    #Position<-apply(df[,c(2,5)], 1, function(x){mean(x, na.rm = T)})
    #EffectSize<-apply(df[,c(3,6)], 1, function(x){mean(x, na.rm = T)})
    df <- df[,c(1,4,7,8)]
    #df <- cbind(df, Position, EffectSize)
    df[is.na(df)]<-0
    df <- cbind(df, D=abs(df$FrequencyPlus-df$FrequencyMinus)) 
    if (i==1){
      gg_frame <- df
    }
    else {
      gg_frame<-rbind(gg_frame, df)
    }
  }
  
  gg_frame[,-1]<-gg_frame[,-1]/100
  
  fit1 <- glm(D~StartingFrequency, data=gg_frame, family=binomial)
  
  p1 <- ggplot(fit1$model, aes_string(x = names(fit1$model)[2], y = names(fit1$model)[1])) + 
    geom_point(alpha=0.05) +
    geom_smooth(method = "glm", method.args = list(family = "binomial"), col="red") +
    # stat_smooth(method = "lm", col = "red") +
    # labs(title = paste("Adj R2 = ",signif(summary(fit1)$adj.r.squared, 5),
    #                    " P =",signif(summary(fit1)$coef[2,4], 5))) +
    scale_y_continuous(limits = c(0,1)) + 
    theme(legend.position="none") +
    xlab("Starting Frequency") +
    ylab("D-value") +
    theme(axis.title = element_text(size = 30)) +
    theme(axis.text=element_text(size=30)) +
    theme(title = element_text(size=20)) +
    theme_bw() +
    theme(axis.text = element_text(size=30), axis.title = element_text(size=30)) +
    theme(text = element_text(size=30)) +
    theme(legend.position=lengend_position) +
    theme(legend.text=element_text(size=30)) +
    theme(legend.background = element_rect(fill = "transparent")) +
    theme(legend.title = element_blank())
  
  setwd("/fs/cbsubscb10/storage/rl683/TemporalScan/Figures/Misc/")
  png(paste0("DvsFreq_", scenario, ".png"), width = 800, height = 750, units = "px", pointsize = 20)
  print(p1)
  dev.off()
}


## Plot D histogram (combining all simulation reps) ####################################################################
setwd("/fs/cbsubscb10/storage/rl683/TemporalScan/")
library(ggplot2)
library(tidyverse)
scenarios=c("NQTL10", "NQTL10_D0", "NQTL10_D1", "NQTL100", "NQTL10_Clustered", "NQTL100_Clustered")
for (j in 1:6) {
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
    df <- cbind(df, D=abs(df$FrequencyPlus-df$FrequencyMinus)) 
    if (i==1){
      gg_frame <- df
    }
    else {
      gg_frame<-rbind(gg_frame, df)
    }
  }
  gg_frame[,6]<-gg_frame[,6]/100
  
  EffectType <- rep("Neutral", dim(gg_frame)[1])
  EffectType[which(gg_frame$EffectSize>0)] <- "QTL"
  EffectType[which(gg_frame$EffectSize<0)] <- "QTL"
  gg_frame <-cbind(gg_frame, EffectType)
  gg_frame <- mutate(gg_frame, D_binned = cut(D, breaks=seq(0,1,0.1), include.lowest=T))
  
  p1 <- ggplot(data=gg_frame, aes(x=as.numeric(D_binned), colour = EffectType, fill=EffectType, group = EffectType)) +
    geom_bar(aes(y=..prop..), position='dodge') +
    scale_color_manual(values=c("Black", "Black"), label=c(" Neutral loci", " QTL")) +
    scale_fill_manual(values=c("Grey", "Black"), label=c(" Neutral loci", " QTL")) +
    xlab("D-value") +
    ylab("Proportion") +
    scale_x_continuous(breaks = c(seq(0, 10, by=1)+0.5),
                       labels = c(seq(0, 1, by=0.1))) +
    coord_cartesian(ylim=c(0,0.8)) +
    theme_bw() +
    theme(axis.text = element_text(size=30), axis.title = element_text(size=30)) +
    theme(text = element_text(size=30)) +
    theme(legend.position=c(0.75, 0.9)) +
    theme(legend.text=element_text(size=30)) +
    theme(legend.background = element_rect(fill = "transparent")) +
    theme(legend.title = element_blank()) +
    theme(panel.grid.minor.x = element_blank())
  
  
  png(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Figures/Misc/D_hist_", scenario, "_all.png"), width = 800, height = 750, units = "px", pointsize = 20)
  print(p1)
  dev.off()
}


## Plot D distribution in one single experiment #########################################################################################
library(ggplot2)
for (i in 1){
  setwd(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Simulations/NQTL100/SimRep", i, "/ExpRepPlus1"))
  df1 <- read.table("SampleFreqAllTimepoints.txt", sep = ",", header = T)
  setwd(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Simulations/NQTL100/SimRep", i, "/ExpRepMinus1"))
  df2 <- read.table("SampleFreqAllTimepoints.txt", sep = ",", header = T)
  
  df1 <- df1[,c(1, 3, 4, 16)]
  colnames(df1) <- c("PermanentID", "Postion", "EffectSize", "FrequencyPlus")
  df2 <- df2[,c(1, 3, 4, 16)]
  colnames(df2) <- c("PermanentID", "Position", "EffectSize", "FrequencyMinus")
  df <- merge(df1, df2, by="PermanentID", all=T)
  Position<-apply(df[,c(2,5)], 1, function(x){mean(x, na.rm = T)})
  EffectSize<-apply(df[,c(3,6)], 1, function(x){mean(x, na.rm = T)})
  df <- cbind(df[,c(1,4,7)], Position, EffectSize)
  df[is.na(df)]<-0
  df <- cbind(df, D=abs(df$FrequencyPlus-df$FrequencyMinus)) 
  EffectType <- rep("Neutral", dim(df)[1])
  EffectType[which(df$EffectSize>0)] <- "Plus"
  EffectType[which(df$EffectSize<0)] <- "Minus"
  df <-cbind(df, EffectType)
  
  # p1 <- ggplot() +
  #   geom_point(data=df[which(df[,5]==0),], aes(x=Position/1000000, y=D/100), colour = "Black", size=2, alpha=0.1) +
  #   geom_point(data=df[which(df[,5]>0),], aes(x=Position/1000000, y=D/100), colour = "Red", size=5) +
  #   geom_point(data=df[which(df[,5]<0),], aes(x=Position/1000000, y=D/100), colour = "Blue", size=5) +
  #   theme(legend.position="none") +
  #   xlab("Position on Chromosome (Mbp)") +
  #   ylab("D-value") +
  #   theme(axis.title = element_text(size = 30)) +
  #   theme(axis.text=element_text(size=30)) +
  #   theme(text = element_text(size=30))
  #   
  p2 <- ggplot() +
    geom_histogram(data=df, aes(x=D/100, y=..density.., colour = EffectType, fill=EffectType), position = "dodge") +
    scale_color_manual(values=c("Blue", "Grey", "Red")) +
    scale_fill_manual(values=c("Blue", "Grey", "Red")) +
    theme(legend.position="none")
  
  # p3 <- ggplot() +
  #   geom_point(data=df[which(df[,5]==0),], aes(x=FrequencyPlus, y=FrequencyMinus), colour = "Black", alpha=0.1) +
  #   geom_point(data=df[which(df[,5]>0),], aes(x=FrequencyPlus, y=FrequencyMinus), colour = "Red", size=3) +
  #   geom_point(data=df[which(df[,5]<0),], aes(x=FrequencyPlus, y=FrequencyMinus), colour = "Blue", size=3) +
  #   theme(legend.position="none")
  
  
  # setwd("/fs/cbsubscb10/storage/rl683/TemporalScan/Figures/Misc/")
  # png(paste0("D_NQTL100", i, ".png"), width = 800, height = 750, units = "px", pointsize = 20)
  # print(p1)
  # dev.off()
  
  setwd("/fs/cbsubscb10/storage/rl683/TemporalScan/Figures/Misc/")
  png(paste0("D_hist_NQTL100_", i, ".png"), width = 800, height = 750, units = "px", pointsize = 20)
  print(p2)
  dev.off()
  
  # setwd("/fs/cbsubscb10/storage/rl683/TemporalScan/Figures/Misc/")
  # png(paste0("EndFrequency_NQTL100", i, ".png"), width = 800, height = 750, units = "px", pointsize = 20)
  # print(p3)
  # dev.off()
  
}



## Plot LD decay ####################################################################

position <-seq(0,3e7, length.out =3001)
p<-4*1000*1e-8*position
LD<-(10+p)/(p^2+13*p+22)
plot(LD~position)
library(ggplot2)
df<- data.frame(position, LD)
ggplot(data=df, aes(x=position, y=LD))
ggplot(data=df, aes(x=position, y=LD)) +
  geom_line() +
  scale_x_continuous(limits = c(0,31000000), breaks = seq(0,3e7,length.out =16)) 
ggplot(data=df, aes(x=position, y=LD)) +
  geom_line() +
  scale_x_continuous(limits = c(0,2000000), breaks = seq(0,2000000,length.out =11)) 



