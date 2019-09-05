## number of generation vs. dominance #####

rm(list=ls(all=TRUE))

setwd("/fs/cbsubscb10/storage/rl683/TemporalScan/")
filenames <- c("NQTL10", "NQTL10_D1", "NQTL10_D0")
n_gens <- c(9, 4)
dominance_cos <- c(0.5, 1, 0)
for (k in 1:2){
  n_gen <- n_gens[k]
  for (i in 1:3){
    path <- paste0("Simulations/", filenames[i], "/")
    dominance_co <- dominance_cos[i]
    filename <- filenames[i]
    for (j in 1:100){
      roc_temp <- read_csv(paste0(path, "SimRep", j, "/ExpRepMinus1/ROC_AllTimepoints_NGen", n_gen+1, "_TransformedD.txt"))[,1:3]
      colnames(roc_temp) <- c("Proportion_of_genetic_variance_in_gen_1", "Proportion_of_QTL_detected_weighted", "False_positive_rate")
      length <- dim(roc_temp)[1]
      roc_temp <- bind_cols(roc_temp, n=1:length, filename=rep(filename,length), n_gen=rep(n_gen, length), dominance_co=rep(dominance_co, length))
      if (i==1&j==1&k==1){
        roc <- roc_temp
      } else {
        roc <- bind_rows(roc, roc_temp)
      }
    }
  }
}


roc_final <- group_by(roc, n, filename, n_gen) %>%
  summarise(Proportion_of_genetic_variance_in_gen_1 = mean(Proportion_of_genetic_variance_in_gen_1), Proportion_of_QTL_detected_weighted = mean(Proportion_of_QTL_detected_weighted), False_positive_rate=mean(False_positive_rate), dominance_co=unique(dominance_co))


p1 <- ggplot(roc_final, aes(False_positive_rate,Proportion_of_QTL_detected_weighted, color=factor(dominance_co, levels = c(0.5,1,0), labels = c("Additive", "Dominant", "Recessive")), linetype=factor(n_gen, levels = c(4,9), labels = c("4", "9")))) +
  geom_line(size = 2) + 
  scale_x_continuous(breaks = seq(0,0.1,0.02)) + 
  scale_y_continuous(limits = c(0,1)) + 
  coord_cartesian(xlim = c(0,0.1)) +
  xlab("False positive rate") +
  ylab("Proportion of QTLs detected") +
  theme_bw() +
  theme(axis.title = element_text(size = 30)) +
  theme(axis.text=element_text(size=30)) +
  theme(text = element_text(size=30)) +
  ggtitle(NULL) +
  theme(legend.position = c(0.415, 0.865)) +
  theme(legend.text=element_text(size=30)) +
  theme(legend.background = element_rect(fill = "transparent", color="black")) +
  guides(color=guide_legend(title="Trait architecture", order = 1), linetype=guide_legend(title="Number of \ngenerations", order = 2)) +
  theme(legend.box = "horizontal") +
  theme(axis.title = element_text(size = 30)) +
  theme(axis.text=element_text(size=30))

png("Figures/Misc/dominance_vs_generation.png", width = 800, height = 750, units = "px", pointsize = 20)
print(p1)
dev.off()

png("~/evolve-resequence-simulation/Figures/figure_9.png", width = 800, height = 750, units = "px", pointsize = 20)
print(p1)
dev.off()

