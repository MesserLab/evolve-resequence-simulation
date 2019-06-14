## This script is used to plot miscellaneous figures for the simulation. It now includes plots of mean trait trajectory, trait/genetic variance trajectory, allele frequency trajectory, and distribution of D-value on the chromosome. More plot will be added in the future. 
## This script is meant to run with interactive R sessions.

rm(list=ls(all=TRUE))
library(ggplot2)
library(reshape)

## A. Plot trait charateristics across all replicates in the standard model ##########################################################
## Prepare data
library(tidyverse)
library(ggplot2)
scenarios <-  c("NQTL10", "NQTL100")
for (scenario in scenarios[2]) {
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
    theme(legend.background = element_rect(fill = "transparent", color="black")) +
    guides(colour = guide_legend(override.aes = list(size=2))) +
    theme(legend.title = element_blank())
  
  png(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Figures/Misc/TraitMean_", scenario, ".png"), width = 800, height = 750, units = "px", pointsize = 20)
  print(p1)
  dev.off()
  
}


## B. Plot phenotype distribution change in one single experiment #########################################################################################
scenarios <-  c("NQTL10", "NQTL100")
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

for (scenario in scenarios[2]) {
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
    
    p2 <- ggplot(data=full_trait, aes(x=GenerationNumber, y=PhenotypeValue, fill=DirectionofSelection)) +
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
      theme(legend.background = element_rect(fill = "transparent", color = "black")) +
      theme(legend.title = element_blank())
    
    png(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Figures/Misc/Trait_dist_", scenario, "_", j, ".png"), width = 800, height = 750, units = "px", pointsize = 20)
    print(p2)
    dev.off()
  }
}


## C. Plot allele frequency trajectory in one single experiment #########################################################################################
library(ggplot2)
library(reshape)

scenarios <- c("NQTL100", "NQTL100_Clustered")
legend_positions <- list(c(0.27, 0.9), c(0.7, 0.9))
for (k in 1) {
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
    
    p3 <- mutate(df, EffectSize=as.factor(EffectSize)) %>%
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
      theme(legend.background = element_rect(fill = "transparent", color = "black")) +
      theme(legend.title = element_blank())
    
    png(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Figures/Misc/FrequencyTrajectory_", scenario, "_", i, ".png"), width = 800, height = 750, units = "px", pointsize = 20)
    print(p3)
    dev.off()
  }
}

## D. Plot D histogram (combining all simulation reps) ####################################################################
setwd("/fs/cbsubscb10/storage/rl683/TemporalScan/")
library(ggplot2)
library(tidyverse)
scenarios=c("NQTL10", "NQTL10_D0", "NQTL10_D1", "NQTL100", "NQTL10_Clustered", "NQTL100_Clustered")
for (j in 4) {
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
  
  p4 <- ggplot(data=gg_frame, aes(x=as.numeric(D_binned), colour = EffectType, fill=EffectType, group = EffectType)) +
    geom_bar(aes(y=..prop..), position='dodge') +
    scale_color_manual(values=c("Black", "Black"), label=c(" Neutral loci", " QTLs")) +
    scale_fill_manual(values=c("Grey", "Black"), label=c(" Neutral loci", " QTLs")) +
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
    theme(legend.background = element_rect(fill = "transparent", color="black")) +
    theme(legend.title = element_blank()) +
    theme(panel.grid.minor.x = element_blank())
  
  
  png(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Figures/Misc/D_hist_", scenario, "_all.png"), width = 800, height = 750, units = "px", pointsize = 20)
  print(p4)
  dev.off()
}

## E. Plot distribution of D on chromosome in one single experiment #######################################################################################

library(ggplot2)
library(reshape)
library(tidyverse)

scenarios <- c("NQTL100", "NQTL100_Clustered")
legend_positions <- list(c(0.64, 0.91), c(0.1, 0.08))
for (k in 1) {
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
    
    p5 <- mutate(df, EffectSize=as.factor(EffectSize), Position=Position/10^6) %>%
      ggplot()  +
      theme_bw() +
      geom_point(aes(x=Position, y=abs(D/100), group=PermanentID, colour = EffectSize, alpha=EffectSize, size=EffectSize)) +
      scale_alpha_manual(values=c(1, 0.3, 1), labels = c(" -1 allele", " Neutral allele", " +1 allele")) +
      scale_size_manual(values=c(5, 1, 5), labels = c(" -1 allele", " Neutral allele", " +1 allele")) +
      scale_color_manual(values=c("Blue", "Black", "Red"), labels = c(" -1 allele", " Neutral allele", " +1 allele")) +
      guides(colour = guide_legend(override.aes = list(size=c(2,1,2), alpha=c(1, 0.5, 1)))) +
      scale_x_continuous(name="Position on chromosome (Mbp)") +
      scale_y_continuous(breaks=seq(0.00,1.00,0.2), labels=seq(0.00,1.00,0.2), name="D-value") +
      coord_cartesian(ylim = c(0,1)) +
      theme(axis.text = element_text(size=30), axis.title = element_text(size=30)) +
      theme(text = element_text(size=30)) +
      theme(legend.position=lengend_position) +
      theme(legend.text=element_text(size=30)) +
      theme(legend.background = element_rect(fill = "transparent", color = "black")) +
      theme(legend.title = element_blank())
    
    
    png(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Figures/Misc/D_vs_position_", scenario, "_", i, ".png"), width = 1600, height = 750, units = "px", pointsize = 20)
    print(p5)
    dev.off()
  }
}

## Assemble these into one figure #########################################
library(cowplot)
figure_1 <- plot_grid(p1, p2, p3, p4, labels=c("A", "B", "C", "D"), nrow = 2, label_size=30) %>%
  plot_grid(., p5, labels = c('', 'E'), ncol = 1, rel_heights = c(2, 1), label_size=30)
png(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Figures/FiguresForPaper/Figure_1.png"), width = 800*2, height = 750*3, units = "px", pointsize = 20)
print(figure_1)
dev.off()

figure_1_up <- plot_grid(p1, p2, p3, labels=c("A", "B", "C"), nrow = 1, label_size=30)
figure_1_down <- plot_grid(p4, p5, labels=c("D", "E"), nrow = 1, label_size=30, rel_widths =  c(1, 2))
figure_1 <- plot_grid(figure_1_up, figure_1_down, nrow = 2)
png(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Figures/FiguresForPaper/Figure_1_v2.png"), width = 800*3, height = 750*2, units = "px", pointsize = 20)
print(figure_1)
dev.off()