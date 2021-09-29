library("methods")
library("foreign")
library("qqman")
library("ggplot2")

#################  Read GWAS result file ########################
MyData<- read.table("rSTATacute_10PCA.assoc.linear", header = TRUE)

#### Manhattan plot using "qqman" package ####

### As default suggestiveline = -log10(1e-05) and  genomewideline = -log10(5e-08)####

jpeg(filename = "Manhattan plot.rSTATacute HNC-UMCG.jpg") ### change file name and main based on your cohort and endpoint 
manhattan(MyData, main = "rSTATacute HNC-UMCG",
          ylim = c(0, 9), cex = 0.6,
          cex.axis = 0.9, col = c("blue4", "orange3"))
dev.off()

#### Function for QQ plot using "ggplot2" package#####
gg_qqplot <- function(ps, ci = 0.95) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed), shape = 1, size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
    geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po)
}


#### Lambda calculation ####
inflation <- function(ps) {
  chisq <- qchisq(1 - ps, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  lambda
}
##### generating QQ plot with 95% CI and Lambda ### change file name and main based on your cohort and endpoint 
jpeg(filename = "QQ plot.rSTATacute HNC-UMCG.jpg")
gg_qqplot(MyData$P) +
  theme_bw(base_size = 24) +
  annotate(
    geom = "text",
    x = -Inf,
    y = Inf,
    hjust = -0.15,
    vjust = 1 + 0.15 * 3,
    label = sprintf("?? = %.2f", inflation(MyData$P)),
    size = 8
  ) +
  labs(title = "rSTATacute HNC-UMCG") +
  theme(
    axis.ticks = element_line(size = 0.5),
    panel.grid = element_blank(),
  )
dev.off()
