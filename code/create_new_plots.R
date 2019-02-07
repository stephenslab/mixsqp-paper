library(ggplot2)
library(cowplot)
dat <- read.csv("../output/simcompare4.csv")
# dat$robj <- pmax(1e-4,dat$robj)
p <- ggplot(dat,aes(x = time,y = robj,color = method)) +
    geom_point() +
    geom_line() +
    scale_y_log10() +
    scale_color_manual(values = c("darkorange","darkblue","dodgerblue",
                                  "firebrick")) +
    labs(x = "runtime (seconds)",
         y = "distance from minimum")
   
