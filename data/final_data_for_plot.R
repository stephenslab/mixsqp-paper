## figure 1

p <- ggplot(data = dat1_1[2:10,]) +
  geom_line(aes(x = log2(n), y=log2(t1),color = "S: simplex"), size = 1.2) +
  geom_line(aes(x = log2(n), y=log2(t2),color = "D: dual"), size = 1.2) +
  geom_line(aes(x = log2(n), y=log2(t3),color = "B: box"), size = 1.2)
p1 <- p + xlab("log2(n)") + ylab("log2(time)")
p1 <- p1 + scale_color_discrete(name = "") + ggtitle("formulation-IP-N-F with m = 10")  +
  theme(legend.position = c(.1,.9),legend.background = element_rect(fill = "transparent"))

p <- ggplot(data = dat1_2[2:10,]) +
  geom_line(aes(x = log2(n), y=log2(t1),color = "S: simplex"), size = 1.2) +
  geom_line(aes(x = log2(n), y=log2(t2),color = "D: dual"), size = 1.2) +
  geom_line(aes(x = log2(n), y=log2(t3),color = "B: box"), size = 1.2)
p2 <- p + xlab("log2(n)") + ylab("log2(time)")
p2 <- p2 + scale_color_discrete(name = "") + ggtitle("formulation-SQP-O-F with m = 40")  +
  theme(legend.position = c(.1,.9),legend.background = element_rect(fill = "transparent"))

multiplot(p1, p2, cols = 2)

## figure 2

p <- ggplot(data = dat2_1) +
  geom_line(aes(x = log2(n), y=log2(t1),color = "F"), size = 1.2) +
  geom_line(aes(x = log2(n), y=log2(t2),color = "SVD"), size = 1.2) +
  geom_line(aes(x = log2(n), y=log2(t3),color = "QR"), size = 1.2)
p <- p + xlab("log2(n)") + ylab("log2(time)")
p1 <- p + scale_color_discrete(name = "") + ggtitle("B-SQP-A-lowrankapprox")  +
  theme(legend.position = c(.1,.9),legend.background = element_rect(fill = "transparent"))

p <- ggplot(data = dat2_2) +
  geom_line(aes(x = log2(n), y=svd,color = "SVD"), size = 1.2) +
  geom_line(aes(x = log2(n), y=qr,color = "QR"), size = 1.2)
p <- p + xlab("log2(n)") + ylab("log10(frobenius error)")
p2 <- p + scale_color_discrete(name = "") + ggtitle("accuracy of low-rank approximation") + ylim(c(-13,-12))  +
  theme(legend.position = c(.1,.9),legend.background = element_rect(fill = "transparent"))
multiplot(p1, p2, cols = 2)

## figure 3

p <- ggplot(data = dat3_1) +
  geom_line(aes(x = log2(m), y= log2(s/m),color = "synthetic"), size = 1.2) +
  geom_line(aes(x = log2(m2), y = log2(s2/m),color = "GIANT"), size = 1.2)
p <- p + xlab("log2(m)") + ylab("log2(r/m)")
p <- p + scale_color_discrete(name = "") + ggtitle("log2(rank/m)") +
  theme(legend.position = c(.8,.95),legend.background = element_rect(fill = "transparent"))
p2 <- ggplot(data = dat3_2) + geom_line(aes(x = log2(m), y = rel_err),color = "#00BA38", size = 1.2)
p3 <- ggplot(data = dat3_3) + geom_line(aes(x = log2(m), y = rel_err ),color = "#619CFF", size = 1.2)
p2 <- p2 + xlab("log2(m)") + ylab("log10(||x_IP - x_SQP||_1)") + ggtitle("difference in l1 norm") + ylim(-8,-4)
p3 <- p3 + xlab("log2(m)") + ylab("log10(|f_IP-f_SQP|/|f_IP|)") + ggtitle("difference in objective") + ylim(-16,-8)
multiplot(p, p2, p3, cols=3)

## figure 4

p <- ggplot(data = dat4_1) +
  geom_line(aes(x = log2(m), y=log2(t1),color = "IP"), size = 1.2) +
  geom_line(aes(x = log2(m), y=log2(t2),color = "Act"), size = 1.2)
p <- p + xlab("log2(m)") + ylab("log2(time)")
p1 <- p + scale_color_discrete(name = "") + ggtitle("comptime of solving QP in m")  +
  theme(legend.position = c(.1,.95),legend.background = element_rect(fill = "transparent"))

p <- ggplot(data = dat4_2) +
  geom_line(aes(x = log2(n), y=log2(t1),color = "IP"), size = 1.2) +
  geom_line(aes(x = log2(n), y=log2(t2),color = "Act"), size = 1.2)
p <- p + xlab("log2(n)") + ylab("log2(time)")
p2 <- p + scale_color_discrete(name = "") + ggtitle("comptime of solving QP in n")  +
  theme(legend.position = c(.9,.6),legend.background = element_rect(fill = "transparent"))

p <- ggplot() +
  geom_line(aes(x = 1:17, y=dat4_3$q_nnz,color = "q_nnzs"), size = 1.2) +
  geom_line(aes(x = 1:17, y=dat4_3$y_nnz,color = "y_nnzs"), size = 1.2) +
  geom_line(aes(x = 1:17, y=dat4_3$linesearch,color = "# of ls"), size = 1.2) 
p <- p + xlab("iteration") + ylab("number")
p3 <- p + scale_color_discrete(name = "") + ggtitle("parameters in each iteration")  +
  theme(legend.position = c(.9,.9),legend.background = element_rect(fill = "transparent"))
multiplot(p1, p2, p3, cols = 3)


## figure 5

p <- ggplot(data = dat5) +
  geom_line(aes(x = n, y=log2(IP100),color = "pink",linetype = "D-IP-N-F"), size = 1.2) +
  geom_line(aes(x = n, y=log2(SQP100),color = "pink",linetype = "B-SQP-A-QR"), size = 1.2) + 
  geom_line(aes(x = n, y=log2(IP200),color = "turquoise",linetype = "D-IP-N-F"), size = 1.2) +
  geom_line(aes(x = n, y=log2(SQP200),color = "turquoise",linetype = "B-SQP-A-QR"), size = 1.2) +
  geom_line(aes(x = n, y=log2(IP400),color = "blue",linetype = "D-IP-N-F"), size = 1.2) +
  geom_line(aes(x = n, y=log2(SQP400),color = "blue",linetype = "B-SQP-A-QR"), size = 1.2) +
  geom_line(aes(x = n, y=log2(IP800),color = "salmon",linetype = "D-IP-N-F"), size = 1.2) +
  geom_line(aes(x = n, y=log2(SQP800),color = "salmon",linetype = "B-SQP-A-QR"), size = 1.2)
p1 <- p + xlab("log2(n)") + ylab("log2(time)")
fig5 <- p1 + scale_color_discrete(name = "m", breaks = c("pink","blue","turquoise","salmon"),
                                  labels = c("100", "200","400","800")) + ggtitle("Computation time of mixSQP versus REBayes")+
  theme(legend.position = c(.1,.75),legend.background = element_rect(fill = "transparent"))
fig5


## figure 6

#### figure 6_1
p1 <- ggplot(dat6_1,aes(x = x,y = y,fill = label)) + ggtitle("mixSQP comptime of each component in ASH") +
  theme(legend.position = c(.23,.87),legend.background = element_rect(fill = "transparent")) + 
  geom_col(position = "stack") + scale_fill_manual("", values = c("#C77CFF","#00BFC4","#7CAE00","#F8766D")) +
  xlab("n") + ylab("time")


#### figure 6_2
p2 <- ggplot(dat6_2,aes(x = x,y = y,fill = label)) + ggtitle("REBayes comptime of each component in ASH") +
  theme(legend.position = c(.23,.87),legend.background = element_rect(fill = "transparent")) + 
  geom_col(position = "stack") + 
  scale_fill_manual("", values = c("posterior computation" = "#C77CFF",
                                   "convex optimization (REBayes)" = "#00BFC4",
                                   "likelihood computation" = "#F8766D")) + xlab("n") + ylab("time")
p1$labels$fill = ""
p2$labels$fill = ""

fig6 <- multiplot(p1, p2, cols=2)
fig6
