library(ggplot2)
library(reshape)
library(tidyverse)
library(cowplot)

scenarios <- c("NQTL10", "NQTL100")
legend_positions <- list(c(0.27, 0.9), c(0.27, 0.9))
for (k in 1:2) {
  scenario <- scenarios[k]
  lengend_position <- legend_positions[[k]]
  for (i in 1){
    for (j in 1:5){
      setwd(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Simulations/", scenario, "/SimRep", i, "/ExpRepPlus", j))
      df1 <- read.table("SampleFreqAllTimepoints.txt", sep = ",", header = T)
      setwd(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Simulations/", scenario, "/SimRep", i, "/ExpRepMinus", j))
      df2 <- read.table("SampleFreqAllTimepoints.txt", sep = ",", header = T)
      
      colnames(df1)[7:16]<-0:9
      df1<- df1[,c(1, 4, 7:16)]
      df1 <- melt(df1, id.vars=c("ID_perm", "es"), value.name="Frequency", variable.name="Timepoint")
      
      colnames(df2)[7:16]<--9:0
      df2[,7:16]<-df2[,16:7]
      df2<- df2[,c(1, 4, 7:16)]
      df2 <- melt(df2, id.vars=c("ID_perm", "es"), value.name="Frequency", variable.name="Timepoint")
      
      df <- rbind(df2, df1)
      df[,4]<-df[,4]/100
      colnames(df) <- c("PermanentID", "EffectSize", "Timepoint", "Frequency")
      
      df_gen4 <- filter(df, Timepoint %in% -4:4)
        
      p_gen4 <- mutate(df_gen4, EffectSize=as.factor(EffectSize)) %>%
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
      
      p_gen10 <- mutate(df, EffectSize=as.factor(EffectSize)) %>%
        ggplot()  +
        theme_bw() +
        geom_line(aes(x=Timepoint, y=Frequency, group=PermanentID, colour = EffectSize, alpha=EffectSize, size=EffectSize)) +
        scale_alpha_manual(values=c(1, 0.012, 1), labels = c(" -1 allele", " Neutral allele", " +1 allele")) +
        scale_size_manual(values=c(1, 0.3, 1), labels = c(" -1 allele", " Neutral allele", " +1 allele")) +
        scale_color_manual(values=c("Blue", "Black", "Red"), labels = c(" -1 allele", " Neutral allele", " +1 allele")) +
        guides(colour = guide_legend(override.aes = list(size=c(2,1,2), alpha=c(1, 0.5, 1)))) +
        annotate("label", label="\"High\" line", x = 6.1*2, y = 0.8, size=10, colour='black', alpha=0.7) +
        annotate("label", label="\"Low\" line", x = 3.9*2, y = 0.8, size=10, colour='black', alpha=0.7) +
        annotate("segment", x=5*2, xend=5*2, y=0*2, yend=1.5, size=2) +
        annotate("segment", x=5.3*2, xend=6.9*2, y=0.75, yend=0.75, size=2, arrow = arrow(length = unit(0.5, "cm")), color ="black") +
        annotate("segment",x=4.7*2, xend=3.1*2, y=0.75, yend=0.75, size=2, arrow = arrow(length = unit(0.5, "cm")), color ="black") +
        scale_x_discrete(breaks=c(-9:9), labels=c(seq(9,0,-1), seq(1,9,1)), name="Generation number") +
        scale_y_continuous(breaks=seq(0.00,1.00,0.2), labels=seq(0.00,1.00,0.2), name="Allele frequency") +
        coord_cartesian(ylim = c(0,1)) +
        theme(axis.text = element_text(size=30), axis.title = element_text(size=30)) +
        theme(text = element_text(size=30)) +
        theme(legend.position=lengend_position) +
        theme(legend.text=element_text(size=30)) +
        theme(legend.background = element_rect(fill = "transparent", color = "black")) +
        theme(legend.title = element_blank())
      assign(paste0("p_gen4_", k, "_", j), p_gen4)
      assign(paste0("p_gen10_", k, "_", j), p_gen10)
    }
  }
}

#figure_s10_1 <- plot_grid(p_gen4_1_1, p_gen4_1_2, p_gen4_1_3, p_gen4_1_4, p_gen4_1_5, labels=c("A", "", "", "", ""), ncol = 1, label_size=30)
figure_s10_2 <- plot_grid(p_gen10_1_1, p_gen10_1_2, p_gen10_1_3, p_gen10_1_4, p_gen10_1_5, labels=c("A", "", "", "", ""), ncol = 1, label_size=30)
#figure_s10_3 <- plot_grid(p_gen4_2_1, p_gen4_2_2, p_gen4_2_3, p_gen4_2_4, p_gen4_2_5, labels=c("C", "", "", "", ""), ncol = 1, label_size=30)
figure_s10_4 <- plot_grid(p_gen10_2_1, p_gen10_2_2, p_gen10_2_3, p_gen10_2_4, p_gen10_2_5, labels=c("B", "", "", "", ""), ncol = 1, label_size=30)

#figure_s10 <- plot_grid(figure_s10_1, figure_s10_2, figure_s10_3, figure_s10_4, ncol=4)
figure_s10 <- plot_grid(figure_s10_2, figure_s10_4, ncol=2)
png(paste0("~/evolve-resequence-simulation/Figures/figure_s10.png"), width = 800*4, height = 750*5, units = "px", pointsize = 20)
print(figure_s10)
dev.off()
