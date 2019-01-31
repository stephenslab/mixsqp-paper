# Script to generate the plots for the paper.

# SET UP ENVIRONMENT
# ------------------
library(ggplot2)
library(cowplot)

# Colors used in some of the plots below.
colors <- c("#E69F00","#56B4E9","#009E73","#F0E442",
            "#0072B2","#D55E00","#CC79A7")

# LOAD RESULTS
# ------------
load("../output/results_for_plots.RData")

# CREATE PLOTS
# ------------
# Prepare the results for the first plot.
pdat <- with(dat1,
             rbind(data.frame(formulation = "Dual(linear constraint)",
                              method      = "IP (JuMP/MOSEK)",
                              n = n,runtime = t1),
                   data.frame(formulation = "Primal(simplex constraint)",
                              method      = "IP (JuMP/MOSEK)",
                              n = n,runtime = t2),
                   data.frame(formulation = "Primal(nonnegative constraint)",
                              method      = "IP (JuMP/MOSEK)",
                              n = n,runtime = t3),
                   data.frame(formulation = "Dual(linear constraint)",
                              method      = "SQP (JuMP/MOSEK)",
                              n = n,runtime = t4),
                   data.frame(formulation = "Primal(simplex constraint)",
                              method      = "SQP (JuMP/MOSEK)",
                              n = n,runtime = t5),
                   data.frame(formulation = "Primal(nonnegative constraint)",
                              method      = "SQP (JuMP/MOSEK)",
                              n = n,runtime = t6),
                   data.frame(formulation = "Dual(linear constraint)",
                              method      = "IP (KWDual/Rmosek)",
                              n = n,runtime = t7)))

# Create a plot comparing the computation time for solving 
# different formulations of the maximum-likelihood estimation problem
# with MOSEK (in JuMP) and the SQP algorithm: (1) the dual problem,
# (2) the primal problem with simple constraints, and (3) the primal
# problem with non-negativity constraints. Also
p1 <- ggplot(data = pdat,aes(x = n,y = runtime,color = method,
                             shape = formulation)) +
  geom_line(size = 0.5) +
  geom_point(size = 1.5) +
  scale_x_continuous(trans = "log10",breaks = c(40,100,1e3,1e4)) +
  scale_y_continuous(trans = "log10",breaks = c(0.01,0.1,1,10,100)) + 
  scale_color_manual(values = c(colors[1:2],"darkblue")) +
  scale_shape_manual(values = c(2,5,19)) +
  labs(x     = "number of rows (n)",
       y     = "runtime (seconds)",
       title = "Complexity of solving different problem formulations") +
  theme_cowplot(font_size = 12) +
  theme(plot.title  = element_text(face = "plain",size = 12),
        axis.line   = element_blank(),
        legend.text = element_text(size = 10))

# Create a plot comparing the runtime of the SQP solver using the full
# matrix L against the same SQP solver using low-rank approximations
# based on RRQR and tSVD factorizations of L.
p3 <- ggplot(data = dat2_1) +
  geom_line(aes(x = n,y = t1,color = "No approx."),size = 1) +
  geom_line(aes(x = n,y = t2,color = "tSVD"),size = 1) +
  geom_line(aes(x = n,y = t3,color = "RRQR"),size = 1) +
  geom_point(aes(x = n,y = t1,color = "No approx."),size = 3,shape = 20) +
  geom_point(aes(x = n,y = t2,color = "tSVD"),size = 3,shape = 20) +
  geom_point(aes(x = n,y = t3,color = "RRQR"),size = 3,shape = 20) +
  scale_x_continuous(trans = "log10",breaks = c(2e3,1e4,1e5,1e6)) +
  scale_y_continuous(trans = "log10",breaks = c(0.01,0.1,1,10,100)) +
  scale_color_manual(values = colors[c(1:2,6)],name = "") +
  labs(x     = "number of data matrix rows (n)",
       y     = "computation time (seconds)",
       title = "SQP with different low-rank approximations") +
  theme_cowplot(font_size = 12) +
  theme(legend.position = c(0.1,0.9),
        plot.title      = element_text(face = "plain",size = 12),
        axis.line       = element_blank())

# Create a plot comparing the accuracy of QR and SVD reconstructions
# of the matrix.
p4 <- ggplot(data = dat2_2) +
  geom_point(aes(x = n,y = 10^svd,color = "SVD",shape = "SVD"),size = 2) +
  geom_point(aes(x = n,y = 10^qr,color = "QR",shape = "QR"),size = 2) +
  scale_x_continuous(trans = "log10",breaks = c(2e3,1e4,1e5,1e6)) +
  scale_color_manual(values = colors[c(2,6)],name = "") +
  scale_shape_manual(values = c(19,4)) +
  labs(x = "number of data matrix rows (n)",
       y = "norm of exact L - approx. L",
       title = "Error in low-rank approximation of L") +
  theme_cowplot(font_size = 12) +
  theme(legend.position = c(0.1,0.9),
        plot.title      = element_text(face = "plain",size = 12),
        axis.line       = element_blank()) +
  guides(color = FALSE,shape = FALSE)

# Prepare the results for the next two plots. In particular, I merge
# the mix-SQP and REBayes results, and change the order of the factor
# levels for a more logical ordering in the plots below.
dat6_1 <- transform(dat6_1,label = as.character(label))
dat6_2 <- transform(dat6_2,label = as.character(label))
dat6_1 <- transform(dat6_1,
                    label = factor(label,
                                   c("QR factorization",
                                     "model fitting (mix-SQP)",
                                     "posterior calculations",
                                     "likelihood computation")))
rows <- which(dat6_2$label == "model fitting (REBayes)")
dat6_2[rows,"label"] <- "model fitting (KWDual)"
pdat <- rbind(transform(dat6_1,
                        x     = x - 5e3,
                        label = as.character(label)),
              transform(dat6_2,
                        x     = as.numeric(as.character(x)) + 5e3,
                        label = as.character(label)))
pdat <- transform(pdat,
                  label = factor(label,c("QR factorization",
                                         "model fitting (KWDual)",
                                         "model fitting (mix-SQP)",
                                         "posterior calculations",
                                         "likelihood computation")))

# Create a plot showing the computation breakdown of adaptive
# shrinkage with the REBayes (MOSEK) and mix-SQP solvers used to
# implement the model fitting.
p5 <- ggplot(pdat,aes(x = x,y = y,fill = label)) +
  geom_col(position = "stack",width = 7e3) +
  scale_fill_manual(name = "",
                    values = c("lightskyblue","orange","orangered","aliceblue",
                               "lightsteelblue")) +
  scale_x_continuous(breaks = seq(5e4,25e4,5e4)) +
  labs(x = "number of data matrix rows (n)",
       y = "computation time (seconds)",
       title = "Breakdown of adaptive shrinkage computation") +
  theme_cowplot(font_size = 12) +
  theme(legend.position = c(.05,0.85),
        plot.title   = element_text(face = "plain",size = 12),
        axis.line    = element_blank(),
        axis.ticks.x = element_blank())

# This is a zoomed-in version of the previous plot, for the mix-SQP
# results only.
p6 <- ggplot(dat6_1,aes(x = x,y = y,fill = label)) +
  geom_col(position = "stack",width = 7e3) +
  scale_fill_manual(name = "",
                    values = c("lightskyblue","orangered","aliceblue","lightsteelblue")) +
  scale_x_continuous(breaks = seq(5e4,25e4,5e4)) +
  labs(x = "number of data matrix rows (n)",
       y = "computation time (seconds)",
       title = "Zoomed version: mix-SQP only") +
  theme_cowplot(font_size = 12) +
  theme(legend.position = c(.05,0.9),
        plot.title   = element_text(face = "plain",size = 12),
        axis.line    = element_blank(),
        axis.ticks.x = element_blank())

# Create a plot showing the runtimes for the mix-SQP and KWDual
# (MOSEK) methods on simulated data sets with different numbers of
# samples (n) and different numbers of mixture components (m).
pdat <- data.frame(n      = rep(2^dat5$n,8),
                   m      = factor(rep(c(100,200,400,800),each = 20)),
                   solver = rep(rep(c("KWDual","mix-SQP"),each = 10),4),
                   time   = do.call(c,dat5[-(1:3)]))
pdat2 <- data.frame(n      = rep(dat5_1$n,2),
                    m      = factor(rep(c(100,200,400,800),2)),
                    solver = rep(c("KWDual","mix-SQP"),each = 4),
                    time   = c(dat5_1$rebayes,dat5_1$mixsqp))
pdat <- rbind(pdat,pdat2)
p7 <- ggplot(data = pdat,aes(x = n,y = time,color = m,shape = solver)) +
  geom_line(size = 0.5) +
  geom_point(size = 2) +
  scale_x_continuous(trans = "log10",breaks = c(2e3,2e4,2e5,2e6)) +
  scale_y_continuous(trans = "log10",breaks = c(0.01,0.1,1,10,100,1e3)) +
  scale_color_manual(values = c("lightskyblue","cornflowerblue",
                                "mediumblue","darkblue"),
                     name = "m (num. cols)") +
  labs(x = "n (number of rows in L)",
       y = "runtime (seconds)",
       title = "Comparison of mix-SQP and KWDual performance") +
  theme_cowplot(font_size = 12) +
  theme(plot.title   = element_text(face = "plain",size = 12),
        axis.line    = element_blank())

# Create a plot comparing the number of columns in the data matrix (L)
# against the "effective" rank of L. NOTE: This plot is not included
# in the manuscript.
p8 <- ggplot(data = dat3_1) +
  geom_line(aes(x = m,y = s,color = "synthetic"), size = 1) +
  geom_line(aes(x = m,y = s2,color = "GIANT"), size = 1) +
  geom_point(aes(x = m,y = s,color = "synthetic"), size = 3,shape = 20) +
  geom_point(aes(x = m,y = s2,color = "GIANT"), size = 3,shape = 20) +
  scale_x_continuous(trans = "log10",limits = c(20,1000),
                     breaks = c(20,100,1000)) +
  scale_y_continuous(breaks = c(0,10,20,30),limits = c(0,31)) +
  scale_color_manual(values = colors,name = "data") +
  labs(x = "number of columns of L (m)",
       y = "numeric rank of L",
       title = "Effective numeric rank of L") +
  theme(plot.title   = element_text(face = "plain",size = 12),
        axis.line    = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.6,0.2))

# Create a plot comparing the runtime of the interior-point (MOSEK)
# and active-set methods for solving the quadratic subproblem, for
# different settings of m (the number of columns in matrix L).
p9 <- ggplot(data = dat4_1) +
  geom_line(aes(x = m,y = t1,color = "interior point (MOSEK)"),size = 1) +
  geom_line(aes(x = m,y = t2,color = "active-set"),size = 1) +
  geom_point(aes(x = m,y = t1,color = "interior point (MOSEK)"),
             shape = 20,size = 3) +
  geom_point(aes(x = m,y = t2,color = "active-set"),shape = 20,size = 3) +
  scale_x_continuous(trans = "log10",breaks = c(10,30,100,500)) + 
  scale_y_continuous(trans = "log10",
                     breaks = c(0.0002,0.001,0.005,0.025,0.125,0.625)) +
  scale_color_manual(values = colors[c(6,2)],name = "") +
  labs(x     = "number of columns in L (m)",
       y     = "runtime (seconds)",
       title = "Average time solving the quadratic subproblem in m") +
  theme(plot.title      = element_text(face = "plain",size = 12),
        axis.line       = element_blank(),
        legend.position = c(0.05,0.9))

# Create a plot comparing the runtime of the interior-point (MOSEK)
# and active-st methods for solving the quadratic subproblem, for
# different settings of n (the number of rows in L).
p10 <- ggplot(data = dat4_2) +
  geom_line(aes(x = n,y = t1,color = "interior point (MOSEK)"),size = 1) +
  geom_line(aes(x = n,y = t2,color = "active-set"),size = 1) +
  geom_point(aes(x = n,y = t1,color = "interior point (MOSEK)"),
             shape = 20,size = 3) +
  geom_point(aes(x = n,y = t2,color = "active-set"),shape = 20,size = 3) +
  scale_x_continuous(trans = "log10",breaks = c(1e3,1e4,4e5)) +
  scale_y_continuous(trans = "log10",breaks = c(0.0005,0.001,0.002,0.004),
                     limits = c(0.0004,0.004)) +
  scale_color_manual(values = colors[c(6,2)],name = "") +
  labs(x     = "number of rows in L (n)",
       y     = "runtime (seconds)",
       title = "Average time solving the quadratic subproblem in n") +
  theme(plot.title      = element_text(face = "plain",size = 12),
        axis.line       = element_blank(),
        legend.position = c(0.05,0.5))

# Create a plot showing the number of backtracking line searches, and
# the number of nonzeros in q and y on convergence of the active set
# method, to illustrate behaviour of the SQP method.
p11 <- ggplot(data = dat4_3[-14,]) +
  geom_line(aes(x = iter,y = y,color = "nonzeros in q"),size = 1) +
  geom_line(aes(x = iter,y = q,color = "nonzeros in y"),size = 1) +
  geom_line(aes(x = iter,y = ls,color = "line search iterations"),size = 1) +
  geom_point(aes(x = iter,y = y,color = "nonzeros in q"),
             shape = 20,size = 3) +
  geom_point(aes(x = iter,y = q,color = "nonzeros in y"),
             shape = 20,size = 3) +
  geom_point(aes(x = iter,y = ls,color = "line search iterations"),
             shape = 20,size = 3) +
  scale_x_continuous(breaks = c(1,5,10,13)) +
  scale_y_continuous(breaks = c(0,5,10)) +
  scale_color_manual(values = c("limegreen","darkblue","skyblue"),
                     name = "") +
  labs(x     = "SQP iteration",
       y     = "count",
       title = "Important numbers in each iteration") +
  theme(plot.title      = element_text(face = "plain",size = 12),
        axis.line       = element_blank(),
        legend.position = c(0.4,0.89))

# Create a plot comparing r/m (effective rank r as estimated by RRQR
# over number of columns) against m (the number of columns).
p13 <- ggplot(data = dat3_1) +
  geom_line(aes(x = m,y = s/m,color = "GIANT data"),size = 1) +
  geom_line(aes(x = m,y = s2/m,color = "Synthetic data"),size = 1) +
  geom_point(aes(x = m,y = s/m, color = "GIANT data"),
             shape = 20,size = 3) +
  geom_point(aes(x = m,y = s2/m,color = "Synthetic data"),
             shape = 20,size = 3) +
  scale_x_continuous(trans = "log2",breaks = c(25,50,100,200,400,800)) +
  scale_y_continuous(trans="log2",breaks=c(0.5,0.25,0.125,0.0625,0.03125,1)) +
  scale_color_manual(values = colors[c(6,2)],name = "") +
  labs(x     = "number of cols in L (m)",
       y     = "effective rank/true rank (rank/m)",
       title = "Ratio of effective rank to true rank") +
  theme(plot.title      = element_text(face = "plain",size = 12),
        axis.line       = element_blank(),
        legend.position = c(0.05,0.3))

# Create a plot showing the L1 norm of differences in solutions
# returned by the interior point and SQP solvers applied to GIANT data
# set with different settings of m (the number of columns of L).
p14 <- ggplot(data = dat3_2) +
  geom_line(aes(x = m,y = rel_err,color = "||x_IP - x_SQP||_1"),size = 1) +
  geom_point(aes(x = m,y = rel_err,color = "||x_IP - x_SQP||_1"),
             shape = 20,size = 3) +
  scale_x_continuous(trans = "log2",breaks = c(25,50,100,200,400,800)) +
  scale_y_continuous(breaks = c(-8,-7,-6,-5),limits = c(-8,-4)) +
  scale_color_manual(values = colors[3],name = "") +
  labs(x     = "number of cols in L (m)",
       y     = "log10 of l1 difference (log10(diff))",
       title = "L1 difference between solutions") +
  theme(plot.title      = element_text(face = "plain",size = 12),
        axis.line       = element_blank(),
        legend.position = c(0.1,0.9))

# Create a plot showing the difference in the objective values at the
# solutions returned by the interior point and SQP solvers applied to
# the GIANT data set, with different settings of m (the number of
# columns of L).
p15 <- ggplot(data = dat3_3) +
  geom_line(aes(x = m,y = err,color = "|f(x_REBayes) - f(x_SQP)|"),
            size = 1) +
  geom_point(aes(x = m,y = err,color = "|f(x_REBayes) - f(x_SQP)|"),
             shape = 20,size = 3) +
  scale_x_continuous(trans = "log2",breaks = c(25,50,100,200,400,800)) +
  scale_y_continuous(breaks = c(-16,-14,-12,-10,-8,-6),limits = c(-16,-6)) +
  scale_color_manual(values = colors[5],name = "") +
  labs(x     = "number of cols in L (m)",
       y     = "log10 of difference (log10(diff))",
       title = "Difference between two objective values") +
  theme(plot.title      = element_text(face = "plain",size = 12),
        axis.line       = element_blank(),
        legend.position = c(0.1,0.9))

# Create plots showing the evolution of the objective for mix-SQP
# vs. two algorithms, projected gradient and EM, that do not benefit
# from the second-order information.
p21 <- ggplot(data = dat9_1,aes(x = time,y = objective,color = method,
              linetype = method)) + 
  geom_line(size = 0.75) +
  scale_x_continuous(limits = c(0,2)) +
  scale_y_continuous(trans = "log10",breaks = c(0.001,0.01,0.1,1,10,100,1e3)) +
  scale_color_manual(values = c("dodgerblue","darkblue","darkorange",
                                "darkorange")) +
  scale_linetype_manual(values = c("solid","solid","solid","dashed")) +
  labs(title = "number of columns (m) = 100",
       x     = "runtime (seconds)",
       y     = "objective value") + 
  theme(plot.title      = element_text(face = "plain",size = 12),
        axis.line       = element_blank())
p23 <- ggplot(data = dat9_3,aes(x = time,y = objective,color = method,
              linetype = method)) + 
  geom_line(size = 0.75) +
  scale_x_continuous(limits = c(0,10)) +
  scale_y_continuous(trans = "log10",breaks = c(0.001,0.01,0.1,1,10,100,1e3)) +
  scale_color_manual(values = c("dodgerblue","darkblue","darkorange",
                                "darkorange")) +
  scale_linetype_manual(values = c("solid","solid","solid","dashed")) +
  labs(title = "number of columns (m) = 400",
       x     = "runtime (seconds)",
       y     = "objective value") + 
  theme(plot.title      = element_text(face = "plain",size = 12),
        axis.line       = element_blank())

# Create plots to show the effect of the RRQR rank on the accuracy of
# the solution.
p18 <- ggplot(data = dat8_1) +
  geom_line(aes(x = x, y = df, color = m), size = 1) +
  geom_point(aes(x = x, y = df, color = m, shape = m), size = 2) +
  labs(x     = "the rank of approximated L (r)",
       y     = "f_approx* - f_true*",
       title = "Sub-optimality in objective vs. rank") +
  scale_x_continuous(limits = c(4,20)) +
  scale_y_continuous(trans = "log10") +
  theme(plot.title      = element_text(face = "plain",size = 12),
        axis.line       = element_blank(),
        legend.position = c(-10,-10))
p20 <- ggplot(data = dat8_1[c(1,3,5:17,18,20,22:34,35,37,39:51),]) +
  geom_line(aes(x = x, y = dx1, color = m), size = 1) +
  geom_point(aes(x = x, y = dx1, color = m, shape = m), size = 2) +
  scale_x_continuous(limits = c(4,20)) +
  scale_y_continuous(trans = "log10") +
  labs(x     = "the rank of approximated L (r)",
       y     = "||x_approx* - x_true*||_1",
       title = "L1-norm difference between solutions") +
  theme(plot.title      = element_text(face = "plain",size = 12),
        axis.line       = element_blank(),
        legend.position = c(0.5,0.8))

# SAVE PLOTS AS PDFs
# ------------------
ggsave("../output/compare-formulations.pdf",p1,height = 4,width = 7)
ggsave("../output/sqp-lowrank-approx.pdf",plot_grid(p3,p4),
       height = 4,width = 9)
ggsave("../output/low-rank-approx-error.pdf",
       plot_grid(p13,p14,p15,nrow = 1),height = 4,width = 12)
ggsave("../output/compare-quadratic-subproblem-solvers.pdf",
       plot_grid(p9,p10,p11,nrow = 1),height = 4,width = 13)
ggsave("../output/mixsqp-vs-kwdual.pdf",p7,height = 4,width = 7)
ggsave("../output/F6.pdf",plot_grid(p5,p6),height = 4,width = 9)
ggsave("../output/low-rank-approx-varying-rank.pdf",
       plot_grid(p18,p20,nrow = 1),height = 4,width = 9)
ggsave("../output/F9.pdf",p7,height = 4,width = 10)
ggsave("../output/mixsqp-vs-first-order.pdf",
       plot_grid(p21,p23,nrow = 1),height = 4,width = 9)
       
