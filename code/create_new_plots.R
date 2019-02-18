library(ggplot2)
library(cowplot)
dat <- read.csv("../output/simcompare2.csv")
dat.matlab <- read.csv("../output/pg-n=1000-m=20.csv.gz",header = FALSE)
dat.matlab <- cbind(dat.matlab,"minConf_SPG")
names(dat.matlab) <- c("time","robj","method")
f <- min(dat$obj - dat$robj)
dat.matlab$robj <- dat.matlab$robj - f
dat <- rbind(dat[c("time","robj","method")],
             dat.matlab)
dat$robj <- pmax(1e-4,dat$robj)
p <- ggplot(dat,aes(x = time,y = robj,color = method)) +
    geom_point() +
    geom_line() +
    xlim(c(0,0.1)) +
    scale_y_log10() +
    scale_color_manual(values = c("darkorange","darkblue","dodgerblue",
                                  "firebrick","magenta")) +
    labs(x = "runtime (seconds)",
         y = "distance from minimum",
         title = "n=1000, m=20")

