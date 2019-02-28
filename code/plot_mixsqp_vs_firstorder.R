# Create a plot comparing the performance of mix-SQP against
# first-order methods (EM & projected gradient). Note that the
# vertical axis in the plot shows the distance to the minimum
# *log-likelihood* (the mix-SQP objective multiplied by n).
library(ggplot2)
library(cowplot)

# The number of rows in the data set.
n <- 20000

# Read and prepare the data for plotting.
dat1 <- read.csv("../output/mixsqp-exact-n=20000-m=800.csv.gz",header = FALSE)
dat2 <- read.csv("../output/mixsqp-approx-n=20000-m=800.csv.gz",header = FALSE)
dat3 <- read.csv("../output/em-n=20000-m=800.csv.gz",header = FALSE)
dat4 <- read.csv("../output/pg-n=20000-m=800.csv.gz",header = FALSE)
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
f   <- min(subset(dat,method == "mix-SQP-exact" |
                      method == "mix-SQP-approx")$objective)
dat <- transform(dat,objective = pmax(n*(objective - f),1e-6))

# Create the plot comparing performance of the methods over time.
# For the 20,000 x 20 data set, use xlim(c(0,0.5)).
# For the 20,000 x 800 data set, use xlim(c(0,46)).
p <- ggplot(dat,aes(x = runtime,y = objective,color = method)) +
  geom_line() +
  geom_point() +
  xlim(c(0,46)) +
  scale_y_continuous(breaks = 10^seq(-6,4,2),trans = "log10") +
  labs(x = "runtime (seconds)",
       y = "distance to minimum",
       title = "n=20,000, m=20")
