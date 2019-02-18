# Create a plot comparing the performance of mix-SQP against
# first-order methods (EM & projected gradient).
library(ggplot2)
library(cowplot)

# Read and prepare the data for plotting.
dat1 <- read.csv("../output/mixsqp-exact-n=20000-m=20.csv.gz",header = FALSE)
dat2 <- read.csv("../output/mixsqp-approx-n=20000-m=20.csv.gz",header = FALSE)
dat3 <- read.csv("../output/em-n=20000-m=20.csv.gz",header = FALSE)
dat4 <- read.csv("../output/pg-n=20000-m=20.csv.gz",header = FALSE)
dat  <- rbind(cbind(data.frame(method    = "mix-SQP-exact",
                               objective = dat1[[1]],
                               runtime   = cumsum(dat1[[2]]))),
              cbind(data.frame(method    = "mix-SQP-approx",
                               objective = dat2[[1]],
                               runtime   = cumsum(dat2[[2]]))),
              cbind(data.frame(method    = "EM",
                               objective = dat3[[1]],
                               runtime   = cumsum(dat3[[2]]))),
              cbind(data.frame(method   = "PG",
                               objective = dat4[[1]],
                               runtime   = cumsum(dat4[[2]]))))
rm(dat1,dat2,dat3,dat4)

# Find the distance between the current solutions and the best solution.
f   <- min(dat$objective)
dat <- transform(dat,objective = objective - f + 1e-8)

# Create the plot comparing performance of the methods over time.
p <- ggplot(dat,aes(x = runtime,y = objective,color = method)) +
  geom_line() +
  geom_point() +
  xlim(c(0,10)) +
  scale_y_log10() +
  labs(x = "runtime (seconds)",
       y = "distance from minimum",
       title = "n=2000, m=20")
