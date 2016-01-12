library(ggplot2)
library(RColorBrewer)
library(gridExtra)

x <- seq(0, 38.5, by = 3.5)
yh <- x^2
yl <- rev(yh)
yi <- -2.67 * x^2 + 106.67 * x + 80
yb <- 2.67 * x^2 - 106.67 * x + 1080

lc <- brewer.pal(4, name = "Set1")

d1 <- as.data.frame(cbind(x, yh, yl, yi, yb))

d1long <- reshape(d1, varying = c("yh", "yl", "yi", "yb"), v.names = "Expression", direction = "long", times = c("High", "Low", "Intermediate", "Bimodal"))

str(d1long)

ggplot(d1long, aes(x = x, y = Expression, group = time)) +
  geom_line() +
  xlab(expression(Temperature ~degree~C)) +
  facet_wrap(~ time) +
  theme(legend.position = "none") +
  scale_colour_brewer(palette="Set1")

ggsave(filename = "../results/expression_types_Fig1.png")
