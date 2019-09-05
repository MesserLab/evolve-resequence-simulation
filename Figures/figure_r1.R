## figure added for revision: comparison among different experimental designs #####

rm(list=ls(all=TRUE))

setwd("/fs/cbsubscb10/storage/rl683/TemporalScan/")
filenames <- c("NQTL10", "NQTL100")
n_qtls <-c(10, 100)
n_gens <- c(9, 4)
for (k in 1:2){
  n_gen <- n_gens[k]
  for (i in 1:2){
    n_qtl <- n_qtls[i]
    filename <- filenames[i]
    path <- paste0("Simulations/", filename, "/")
    for (j in 1:100){
      roc_temp_1 <- read_csv(paste0(path, "SimRep", j, "/ExpRepMinus1/ROC_AllTimepoints_NGen", n_gen+1, "_TransformedD.txt"))[,1:3] %>%
        mutate(n_rep=1)
      roc_temp_2 <- read_csv(paste0(path, "SimRep", j, "/ExpRepMinus2/ROC_AllTimepoints_NGen", n_gen+1, "_TransformedD.txt"))[,1:3] %>%
        mutate(n_rep=2)
      roc_temp_5 <- read_csv(paste0(path, "SimRep", j, "/ExpRepMinus5/ROC_AllTimepoints_NGen", n_gen+1, "_TransformedD.txt"))[,1:3] %>%
        mutate(n_rep=5)
      roc_temp <- rbind(roc_temp_1, roc_temp_2, roc_temp_5)
      colnames(roc_temp) <- c("Proportion_of_genetic_variance_in_gen_1", "Proportion_of_QTL_detected_weighted", "False_positive_rate", "n_rep")
      length <- dim(roc_temp)[1]/3
      roc_temp <- bind_cols(roc_temp, n=rep(1:length,3), filename=rep(filename,length*3), n_gen=rep(n_gen, length*3), n_qtl=rep(n_qtl, length*3))
      if (i==1&j==1&k==1){
        roc <- roc_temp
      } else {
        roc <- bind_rows(roc, roc_temp)
      }
    }
  }
}


roc_final <- group_by(roc, n, filename, n_gen, n_rep, n_qtl) %>%
  summarise(Proportion_of_genetic_variance_in_gen_1 = mean(Proportion_of_genetic_variance_in_gen_1), Proportion_of_QTL_detected_weighted = mean(Proportion_of_QTL_detected_weighted), False_positive_rate=mean(False_positive_rate))

i=1
for (j in c(10, 100)){
  for (k in c(4, 9)){
    p <- filter(roc_final, n_qtl==j, n_gen==k) %>%
      ggplot(aes(False_positive_rate,Proportion_of_QTL_detected_weighted, color=factor(n_rep, levels = c(1,2,5), labels = c("1", "2", "5")))) +
      geom_line(size = 2, alpha = 0.8, linetype="dotted") + 
      geom_line(data=filter(roc_final, n_qtl==j, n_gen==k, False_positive_rate>0, Proportion_of_QTL_detected_weighted>0), size = 2, alpha = 1) + 
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
      theme(legend.position = c(0.2, 0.85)) +
      theme(legend.text=element_text(size=30)) +
      theme(legend.background = element_rect(fill = "transparent", color="black")) +
      guides(color=guide_legend(title="Number of \nreplications", order = 1)) +
      theme(legend.box = "horizontal") +
      theme(axis.title = element_text(size = 30)) +
      theme(axis.text=element_text(size=30))
    assign(paste0("p", i), p)
    i=i+1
  }
}

library(cowplot)
figure_r1 <- plot_grid(p1, p2, p3, p4, labels=c("A", "B", "C", "D"), nrow = 2, label_size=30)

png("~/evolve-resequence-simulation/Figures/figure_r1.png", width = 800*2, height = 750*2, units = "px", pointsize = 20)
print(figure_r1)
dev.off()

