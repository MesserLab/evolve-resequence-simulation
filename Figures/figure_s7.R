## Plot D ~ Starting Frequency (combining all simulation reps) ####################################################################
library(ggplot2)

scenarios <- c("NQTL100", "NQTL100_D0", "NQTL100_D1", "NQTL10", "NQTL10_D0", "NQTL10_D1", "NQTL10_FreqHigh5", "NQTL10_FreqLow5", "NQTL100_FreqHigh5", "NQTL100_FreqLow5", "NQTL100_ESDistE", "NQTL10_ESDistE")

for (k in c(4,7)) {
  scenario <- scenarios[k]
  print(scenario)
  for (i in 1:100){
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
    df <- mutate(df, D=abs(2*asin(sqrt(FrequencyPlus/100))-2*asin(sqrt(FrequencyMinus/100)))/pi, StartingFrequency=StartingFrequency/100) 
    if (i==1){
      gg_frame <- df
    }
    else {
      gg_frame<-rbind(gg_frame, df)
    }
  }
  
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
    ylab("Transformed D") +
    theme(axis.title = element_text(size = 30)) +
    theme(axis.text=element_text(size=30)) +
    theme(title = element_text(size=20)) +
    theme_bw() +
    theme(axis.text = element_text(size=30), axis.title = element_text(size=30)) +
    theme(text = element_text(size=30)) +
    theme(legend.text=element_text(size=30)) +
    theme(legend.background = element_rect(fill = "transparent")) +
    theme(legend.title = element_blank())
  
  setwd("/fs/cbsubscb10/storage/rl683/TemporalScan/Figures/Misc/")
  png(paste0("DvsFreq_", scenario, ".png"), width = 800, height = 750, units = "px", pointsize = 20)
  print(p1)
  dev.off()
  
  assign(paste0("panel_", k), p1)
}

library(cowplot)
figure_s7 <- plot_grid(panel_4, panel_7, labels=c("A", "B"), nrow = 1, label_size=30)
png(paste0("~/evolve-resequence-simulation/Figures/figure_s7.png"), width = 800*2, height = 750, units = "px", pointsize = 20)
print(figure_s7)
dev.off()
