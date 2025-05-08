
# rd2d: illustration file
# Authors: M. D. Cattaneo, R. Titiunik, R. R. Yu

rm(list=ls(all=TRUE))

library(MASS)
library(ggplot2)
library(rdrobust)
library(latex2exp)
library(tidyr)
library(dplyr)
library(haven)
library(xtable)
library(expm)
library(rd2d)

################################## Load Data ###################################

data <- read.csv("spp.csv")
data$X <- NULL
colnames(data) <- c("x.1", "x.2","y","d")
na.ok <- complete.cases(data$x.1) & complete.cases(data$x.2)
data <- data[na.ok,]

neval <- 40
eval <- matrix(nrow = neval, ncol = 2)
for (i in 1: ceiling(neval * 0.5)){
  eval[i,] <- c(0, 50 - (i-1) * 50 / ceiling(neval * 0.5))
}
for (i in (ceiling(neval * 0.5)+1): neval){
  eval[i,] <- c((i - ceiling(neval * 0.5) - 1) *50 / (ceiling(neval * 0.5)),0)
}
eval <- data.frame(eval)
colnames(eval) <- c("x.1", "x.2")

####################### Scatter Plot and Boundary ############################## 

bound <- 40
eval_labeled <- eval %>%
  mutate(
    index = row_number(),
    label = as.list(if_else(
      index %in% c(1, 5, 10,15, 21,  25, 30, 35, 40),
      paste0("$\\textbf{b}_{", index, "}$"),  # e.g. "$x_{10}$"
      NA_character_
    ))
  )

eval_labeled_bound <- eval_labeled[eval_labeled$index <= bound,]

# Remove label for index 31
eval_labeled_bound_no_21 <- eval_labeled_bound %>%
  filter(index != 21)
annotation_data <- eval_labeled_bound_no_21 %>%
  filter(!is.na(label))

# Extract the coordinates of point 31
point_21 <- eval_labeled_bound %>% filter(index == 21)
x_21 <- point_21$x.1
y_21 <- point_21$x.2

p1.1 <- ggplot() +
  geom_point(
    data = data  %>% sample_frac(0.3),
    aes(x = x.1, y = x.2, color = factor(d)),
    alpha = 0.5, size = 0.5
  ) +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = 60),
               linetype = "solid", color = "grey", alpha = 1, size = 3) +
  geom_segment(aes(x = 0, xend = 60, y = 0, yend = 0),
               linetype = "solid", color = "grey", alpha = 1, size = 3) +
  
  # Plot evaluation points excluding index = 21
  geom_point(
    data = eval_labeled_bound_no_21,
    aes(x = x.1, y = x.2),  # Map color to factor(d)
    alpha = 1,
    size = 0.5
  ) +
  
  # Use scale_color_manual to define your own colors for d=0 and d=1
  scale_color_manual(
    name = NULL,  # Legend title (optional)
    values = c("0" = "#619CFF",  # Color for d=0
               "1" = "#F8766D"), # Color for d=1
    labels = c("0" = "Control", "1" = "Treatment")
  )

annotation_data <- eval_labeled_bound_no_21 %>%
  filter(!is.na(label))

# Loop through each row and add an annotation using TeX() for LaTeX parsing
for(i in seq_len(nrow(annotation_data))) {
  if (i <= 4){
    hjust <- -0.15
    vjust <- 0.1
  } else {
    hjust <- 0.1
    vjust <- -0.15
  }
  p1.1 <- p1.1 + annotate("text",
                          x = annotation_data$x.1[i],
                          y = annotation_data$x.2[i],
                          label = TeX(as.character(annotation_data$label[i])),
                          hjust = hjust,
                          vjust = vjust,
                          size = 6,
                          color = "black",
                          fontface = "bold")}

# Arrow pointing to point 21
p1.1 <- p1.1 +  geom_segment(
  aes(x = x_21 - 8, y = y_21 - 8, xend = x_21, yend = y_21),
  arrow = arrow(length = unit(0.2, "cm")),
  color = "black",
  size = 1
) +
  
  # Label next to the arrow
  annotate(
    "text",
    x = x_21 - 10,
    y = y_21 - 12,
    label = TeX("$\\textbf{b}_{21}$"),
    color = "black",
    size = 6,
    fontface = "bold"
  ) +

  annotate(
    "text",
    x = eval$x.1[40] + 10,
    y = eval$x.2[40] - 5,
    label = "Boundary",
    color = "black",
    size = 6,
    fontface = "bold"
  ) +
  
  labs(x = "Saber 11 Score", y = "Sisben Score") +
  coord_cartesian(xlim = c(-80, 100), ylim = c(-40, 60)) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    plot.title   = element_text(size = 20, hjust = 0.5),
    text = element_text(family="Times New Roman", face="bold"),
    axis.text.x  = element_text(face = "bold", size = 12),
    axis.text.y  = element_text(face = "bold", size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.2, 0.8),
    legend.background = element_rect(fill = "white", color = "black", linetype = "solid")
  ) +
  xlab("Saber 11") +
  ylab("Sisben")

print(p1.1)

ggsave("Results/fig1a-scatter.png", p1.1, width = 6, height = 5)

################################ Scale Data ####################################

scale <- TRUE
if (scale){
  x.1.scale <- 1
  x.2.scale <- sd(data$x.1) / sd(data$x.2)
  data$x.1 <- data$x.1 * x.1.scale
  data$x.2 <- data$x.2 * x.2.scale
  eval$x.1 <- eval$x.1 * x.1.scale
  eval$x.2 <- eval$x.2 * x.2.scale
} else{
  x.1.scale <- 1
  x.2.scale <- 1
}

############################# SPP: Using Bivariate Method ######################

Y <- data$y
X <- cbind(data$x.1, data$x.2)
t <- data$d
b <- eval
result.rd2d <- rd2d(Y, X, t, b, stdvars = FALSE)
summary(result.rd2d, CBuniform = TRUE, subset = c(1,5,10,15,21,25,30,35,40))

tau.hat.biv <- result.rd2d$tau.hat
CI.lower.biv <- result.rd2d$cb$CI.l
CI.upper.biv <- result.rd2d$cb$CI.r
CB.lower.biv <- result.rd2d$cb$CB.l
CB.upper.biv <- result.rd2d$cb$CB.r

############################# SPP: Using Distance Method #######################

D <- proxy::dist(X, eval, method = "euclidean")  # Use "euclidean" for Euclidean distances
d_expanded <- matrix(rep(2 * d - 1, times = ncol(D)), nrow = nrow(D), ncol = ncol(D))
D <- D * d_expanded

# kink off
result.dist.kinkoff <- rd2d.dist(Y,D,kink = "off")
summary(result.dist.kinkoff,CBuniform = FALSE, subset = c(1,5,10,15,21,25,30,35,40))
tau.hat.dist.kinkoff <- result.dist.kinkoff$tau.hat

# kink on
result.dist.kinkon <- rd2d.dist(Y,D,kink = "on")
summary(result.dist.kinkon,CBuniform = FALSE, subset = c(1,5,10,15,21,25,30,35,40))
tau.hat.dist.kinkon <- result.dist.kinkon$tau.hat

############################# SPP: Point Estimation ############################

indx <- c(1:neval)
df <- data.frame(
  indx = rep(indx, 3),
  y = c(tau.hat.biv, tau.hat.dist.kinkon, tau.hat.dist.kinkoff),
  label = rep(c("Bivariate","Dist (kink)", "Dist (no kink)"), each = length(indx))
)

temp_plot <- ggplot() + theme_bw()
temp_plot <- temp_plot + geom_point(data = df, aes(x = indx, y = y, color = label, shape = label))
temp_plot <- temp_plot + theme_minimal() +
  theme(
    axis.title.x = element_text(size = 15,face = "bold"),  # X-axis label size
    axis.title.y = element_text(size = 15,face = "bold"),  # Y-axis label size
    plot.title = element_text(size = 20, hjust = 0.5),  # Title size and centering
    text=element_text(family="Times New Roman", face="bold"),
    axis.text.x = element_text(face = "bold", size = 15),
    axis.text.y = element_text(face = "bold", size = 12),
    legend.position = c(0.76, 1),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", colour = NA),  # White background, no border
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

temp_plot <- temp_plot + scale_color_manual(
  values = c("Bivariate" = "blue", "Dist (kink)" = "grey", "Dist (no kink)" = "red"), name = NULL
) + scale_shape_manual(
  values = c("Bivariate" = 16, "Dist (kink)" = 1, "Dist (no kink)" = 4),  name = NULL
) + guides(
  color = guide_legend(order = 1),  # Keep "Uniform CB" first
  shape = guide_legend(order = 1)  # Ensure everything is in one group
)


temp_plot <- temp_plot + coord_cartesian(xlim = c(1, 40))
temp_plot <- temp_plot + xlab("Cutoffs on the Boundary") + ylab("Treatment Effect")

# Adjust legend
legend_order <- c("Bivariate","Dist (kink)", "Dist (no kink)")

temp_plot <- temp_plot + geom_vline(xintercept = c(21), color = "lightgrey", size = 1, linetype = "dotted")

temp_plot <- temp_plot +
  annotate(
    "text",
    x = 19,                # Same x as the vertical line
    y = 0.05,               # Adjust this so it appears where you want
    label = "kink",
    color = "dimgrey",
    size = 4,              # Text size
    vjust = -0.5,          # Vertical justification (pulls text below this y-value)
    fontface = "bold"
  )

temp_plot <- temp_plot + scale_x_continuous(
  breaks = c(1,5,10,15,21,25,30,35,40),
  labels = c(TeX("$\\textbf{b}_{1}$"), TeX("$\\textbf{b}_{5}$"), TeX("$\\textbf{b}_{10}$"), 
             TeX("$\\textbf{b}_{15}$"),TeX("$\\textbf{b}_{21}$"), TeX("$\\textbf{b}_{25}$"), 
             TeX("$\\textbf{b}_{30}$"), TeX("$\\textbf{b}_{35}$"),TeX("$\\textbf{b}_{40}$"))
)

temp_plot <- temp_plot + coord_cartesian(xlim = c(1, 40),ylim = c(0.05,0.55)) # ylim depends on confidence bands


# Print the plot
print(temp_plot)

ggsave("Results/fig3a-point-estimation.png", temp_plot, width = 6, height = 5)

############################# SPP: Confidence Bands ############################

indx <- c(1:neval)

# Create a data frame for plotting
bound <- 40

df <- data.frame( indx = indx, y = tau.hat.biv, label = rep(c("Bivariate"), each = length(indx)))

# Build the plot
temp_plot <- ggplot() + theme_bw()
df <- df[df$indx <= bound,]

# Scatter plot for "Bivariate" and "Distance"
temp_plot <- temp_plot + geom_point(data = df,
             aes(x = indx, y = y, color = label, shape = label, fill = label, linetype = label))

df_ribbon <- data.frame(
  indx = indx[c(1:bound)],
  ymin = CB.lower.biv[c(1:bound)],
  ymax = CB.upper.biv[c(1:bound)],
  label = "Uniform CB" 
)

temp_plot <- temp_plot + geom_ribbon(data = df_ribbon, aes(x = indx, ymin = ymin, ymax = ymax,
             color = label, shape = label, fill = label, linetype = label), alpha = 0.1)

df_errorbar <- data.frame(
  indx = indx[c(1:bound)],
  ymin = CI.lower.biv[c(1:bound)],
  ymax = CI.upper.biv[c(1:bound)],
  label = "Pointwise CI" 
)

temp_plot <- temp_plot + geom_errorbar(data = df_errorbar, aes(x = indx, ymin = ymin, ymax = ymax,
             color = label, shape = label, fill = label, linetype = label))

temp_plot <- temp_plot + xlab("Cutoffs on the Boundary") + ylab("Treatment Effect")

legend_order <- c("Bivariate", "Pointwise CI", "Uniform CB")

temp_plot <- temp_plot + scale_color_manual(
  values = c("Bivariate" = "blue", "Pointwise CI" = "blue","Uniform CB" = "dodgerblue4"),
  name = NULL, breaks = legend_order
) + scale_shape_manual(
  values = c("Bivariate" = 16, "Pointwise CI" = 124,"Uniform CB" = 0),  # 16 = filled circle, 17 = triangle
  name = NULL, breaks = legend_order
) + scale_fill_manual(
  values = c("Bivariate" = NA, "Pointwise CI" = NA, "Uniform CB" = "dodgerblue4"),
  name = NULL, breaks = legend_order
) + scale_linetype_manual(
  values = c("Bivariate" = 0, "Pointwise CI" = 5,"Uniform CB" = 0),
  name = NULL, breaks = legend_order) +
  guides(
    fill = guide_legend(order = 1),  # Keep "Uniform CB" first
    color = guide_legend(order = 1),  # Ensure everything is in one group
    shape = guide_legend(order = 1),   # Align shapes with color
    linetype = guide_legend(order = 1)
  )

temp_plot <- temp_plot + geom_vline(xintercept = c(21), color = "lightgrey", size = 1, linetype = "dotted")

temp_plot <- temp_plot +
  annotate(
    "text",
    x = 19,                # Same x as the vertical line
    y = 0.05,               # Adjust this so it appears where you want
    label = "kink",
    color = "dimgrey",
    size = 4,              # Text size
    vjust = -0.5,          # Vertical justification (pulls text below this y-value)
    fontface = "bold"
  )

# Place legend inside and adjust text sizes
temp_plot <- temp_plot + theme_minimal() +
  theme(
    axis.title.x = element_text(size = 15,face = "bold"),  # X-axis label size
    axis.title.y = element_text(size = 15,face = "bold"),  # Y-axis label size
    plot.title = element_text(size = 20, hjust = 0.5),  # Title size and centering
    text=element_text(family="Times New Roman", face="bold"),
    axis.text.x = element_text(face = "bold",
                               size = 15),
    axis.text.y = element_text(face = "bold",
                               size = 12),
    legend.position = c(0.8, 1),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", colour = NA),  # White background, no border
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

temp_plot <- temp_plot + scale_x_continuous(
  breaks = c(1,5,10,15,21,25,30,35,40),
  labels = c(TeX("$\\textbf{b}_{1}$"), TeX("$\\textbf{b}_{5}$"), TeX("$\\textbf{b}_{10}$"), 
             TeX("$\\textbf{b}_{15}$"),TeX("$\\textbf{b}_{21}$"), TeX("$\\textbf{b}_{25}$"), 
             TeX("$\\textbf{b}_{30}$"), TeX("$\\textbf{b}_{35}$"), TeX("$\\textbf{b}_{40}$"))
)

temp_plot <- temp_plot + coord_cartesian(xlim = c(1, 40), ylim = c(0.05, 0.55))

# Print the plot
print(temp_plot)

ggsave("Results/fig3b-ci-and-cb.png", temp_plot, width = 6, height = 5)

############################ SPP: heatmap ######################################

eval.original <- eval
scale <- TRUE
if (scale){
  eval.original$x.1 <- eval$x.1 / x.1.scale
  eval.original$x.2 <- eval$x.2 / x.2.scale
}

# heat map for treatment effect

data.plot <- cbind(eval.original, tau.hat.biv)
colnames(data.plot) <- c("x.1", "x.2", "tau.hat")

# Function to interpolate points and color values between two consecutive points
interpolate_points <- function(df, n=ninter){
  do.call(rbind, lapply(1:(nrow(df)-1), function(i){
    xseq <- seq(df$x.1[i], df$x.1[i+1], length.out = n+2)[2:(n+1)]
    yseq <- seq(df$x.2[i], df$x.2[i+1], length.out = n+2)[2:(n+1)]
    colorseq <- seq(df$tau.hat[i], df$tau.hat[i+1], length.out = n+2)[2:(n+1)]
    data.frame(x.1 = xseq, x.2 = yseq, tau.hat = colorseq)
  }))
}

# Generating interpolated points
ninter <- 10
interpolated_data <- interpolate_points(data.plot)

# Plotting

heat_wd <- 6.5
heatcol_low <- "blue"
heatcol_high <- "red"
heat_lab <- TeX("$\\tau(\\textbf{b})$")
heat_title <- NULL
xlabel <- "Saber11"
ylabel <- "Sisben"

augmented_data <-rbind(data.plot, interpolated_data)
ord <- order(augmented_data[,1], augmented_data[,2])
augmented_data <- augmented_data[ord,]

if (is.null(heat_title)) heat_title <- "Heat Map"
plot_heat <- ggplot(data.plot, aes(x=x.1, y=x.2)) +
  geom_segment(data = augmented_data,
               aes(xend = lead(x.1, order_by=x.1),
                   yend = lead(x.2, order_by=x.2),
                   color = tau.hat),
               size = heat_wd, lineend = "round") +
  scale_color_gradient(
    low=heatcol_low, 
    high=heatcol_high) +
  labs(color = heat_lab) +
  xlab(xlabel) +
  ylab(ylabel)

plot_heat <- plot_heat + geom_text(data=data.plot, aes(x=x.1, y=x.2, label=sprintf("%02d", 1:nrow(data.plot))), color="black", size=3)

plot_heat <- plot_heat + theme_minimal() +
  theme(
    axis.title.x = element_text(size = 15,face = "bold"),  # X-axis label size
    axis.title.y = element_text(size = 15,face = "bold"),  # Y-axis label size
    plot.title = element_text(size = 20, hjust = 0.5),  # Title size and centering
    text=element_text(family="Times New Roman", face="bold"),
    axis.text.x = element_text(face = "bold",
                               size = 15),
    axis.text.y = element_text(face = "bold",
                               size = 12),
    legend.position = c(0.8, 1),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", colour = NA),  # White background, no border
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + coord_cartesian(xlim = c(-10, 55), ylim = c(-10, 55))

print(plot_heat)

ggsave("Results/heat-spp.png", plot_heat, width = 6, height = 5)

# heatmap for p-value

data.plot$p.value <- result.rd2d$pvalues
data.plot$p.sig <- cut(data.plot$p.value,
                       breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                       labels = c("p < 0.001", "0.001 ≤ p < 0.01", "0.01 ≤ p < 0.05", "0.05 ≤ p < 0.1", "p ≥ 0.1"))
sig_colors <- c("p < 0.001" = "#d73027",       # red
                "0.001 ≤ p < 0.01" = "#fc8d59", # orange
                "0.01 ≤ p < 0.05" = "#fee08b",  # yellow
                "0.05 ≤ p < 0.1" = "#d9ef8b",   # light green
                "p ≥ 0.1" = "#91cf60")          # green
library(ggplot2)

plot_heat_pvalue <- ggplot(data.plot, aes(x = x.1, y = x.2, fill = p.sig)) +
  geom_tile(color = "white",show.legend = TRUE) +
  scale_fill_manual(values = sig_colors, name = "P-value",drop = FALSE) +
  labs(x = "Saber11", y = "Sisben") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    plot.title = element_text(size = 20, hjust = 0.5),
    text = element_text(family = "Times New Roman", face = "bold"),
    axis.text.x = element_text(face = "bold", size = 15),
    axis.text.y = element_text(face = "bold", size = 12),
    legend.position = c(0.8, 1),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) 

plot_heat_pvalue <- plot_heat_pvalue +  coord_cartesian(xlim = c(-10, 55), ylim = c(-10, 55)) 
plot_heat_pvalue <- plot_heat_pvalue + geom_text(data=data.plot, aes(x=x.1, y=x.2, label=sprintf("%02d", 1:nrow(data.plot))), color="black", size=3)

print(plot_heat_pvalue)

ggsave("Results/heat-pvalue.png", plot_heat_pvalue, width = 6, height = 5)

############################# Placebo: Using Bivariate Method ##################

data_placebo <- read_dta("spp.dta")
covariate <- "icfes_educm1"
data_placebo <- data_placebo[c("running_saber11","running_sisben",covariate)] 
colnames(data_placebo) <- c("x.1", "x.2","y")
na.ok <- complete.cases(data_placebo)
data_placebo <- data_placebo[na.ok,]
data_placebo <- as.data.frame(data_placebo)
data_placebo$d <- as.integer(data_placebo$x.1 >= 0 & data_placebo$x.2 >= 0)

scale <- TRUE
if (scale){
  x.1.scale <- 1
  x.2.scale <- sd(data_placebo$x.1) / sd(data_placebo$x.2)
  data_placebo$x.1 <- data_placebo$x.1 * x.1.scale
  data_placebo$x.2 <- data_placebo$x.2 * x.2.scale
} else{
  x.1.scale <- 1
  x.2.scale <- 1
}

Y <- data_placebo$y
X <- cbind(data_placebo$x.1, data_placebo$x.2)
d <- data_placebo$d
b <- eval
placebo.rd2d <- rd2d(Y, X, t, b)
summary(placebo.rd2d)

tau.hat.biv <- placebo.rd2d$tau.hat
CI.lower.biv <- placebo.rd2d$cb$CI.l
CI.upper.biv <- placebo.rd2d$cb$CI.r
CB.lower.biv <- placebo.rd2d$cb$CB.l
CB.upper.biv <- placebo.rd2d$cb$CB.r

############################# Placebo: Using Distance Method ###################

D <- proxy::dist(X, eval, method = "euclidean")  # Use "euclidean" for Euclidean distances
d_expanded <- matrix(rep(2 * d - 1, times = ncol(D)), nrow = nrow(D), ncol = ncol(D))
D <- D * d_expanded

# kink off
result.dist.kinkoff <- rd2d.dist(Y,D,kink = "off")
summary(result.dist.kinkoff)
tau.hat.dist.kinkoff <- result.dist.kinkoff$tau.hat

# kink on
result.dist.kinkon <- rd2d.dist(Y,D,kink = "on")
summary(result.dist.kinkon)
tau.hat.dist.kinkon <- result.dist.kinkon$tau.hat

############################# Placebo: Point Estimation ############################

indx <- c(1:neval)
df <- data.frame(
  indx = rep(indx, 3),
  y = c(tau.hat.biv, tau.hat.dist.kinkon, tau.hat.dist.kinkoff),
  label = rep(c("Bivariate","Dist (kink)", "Dist (no kink)"), each = length(indx))
)

temp_plot <- ggplot() + theme_bw()
temp_plot <- temp_plot + geom_point(data = df, aes(x = indx, y = y, color = label, shape = label))
temp_plot <- temp_plot + theme_minimal() +
  theme(
    axis.title.x = element_text(size = 15,face = "bold"),  # X-axis label size
    axis.title.y = element_text(size = 15,face = "bold"),  # Y-axis label size
    plot.title = element_text(size = 20, hjust = 0.5),  # Title size and centering
    text=element_text(family="Times New Roman", face="bold"),
    axis.text.x = element_text(face = "bold", size = 15),
    axis.text.y = element_text(face = "bold", size = 12),
    legend.position = c(0.76, 1),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", colour = NA),  # White background, no border
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

temp_plot <- temp_plot + scale_color_manual(
  values = c("Bivariate" = "blue", "Dist (kink)" = "grey", "Dist (no kink)" = "red"), name = NULL
) + scale_shape_manual(
  values = c("Bivariate" = 16, "Dist (kink)" = 1, "Dist (no kink)" = 4),  name = NULL
) + guides(
  color = guide_legend(order = 1),  # Keep "Uniform CB" first
  shape = guide_legend(order = 1)  # Ensure everything is in one group
)


temp_plot <- temp_plot + coord_cartesian(xlim = c(1, 40))
temp_plot <- temp_plot + xlab("Cutoffs on the Boundary") + ylab("Treatment Effect")

# Adjust legend
legend_order <- c("Bivariate","Dist (kink)", "Dist (no kink)")

temp_plot <- temp_plot + geom_vline(xintercept = c(21), color = "lightgrey", size = 1, linetype = "dotted")

temp_plot <- temp_plot +
  annotate(
    "text",
    x = 19,                # Same x as the vertical line
    y = -0.16,               # Adjust this so it appears where you want
    label = "kink",
    color = "dimgrey",
    size = 4,              # Text size
    vjust = -0.5,          # Vertical justification (pulls text below this y-value)
    fontface = "bold"
  )

temp_plot <- temp_plot + scale_x_continuous(
  breaks = c(1,5,10,15,21,25,30,35,40),
  labels = c(TeX("$\\textbf{b}_{1}$"), TeX("$\\textbf{b}_{5}$"), TeX("$\\textbf{b}_{10}$"), 
             TeX("$\\textbf{b}_{15}$"),TeX("$\\textbf{b}_{21}$"), TeX("$\\textbf{b}_{25}$"), 
             TeX("$\\textbf{b}_{30}$"), TeX("$\\textbf{b}_{35}$"),TeX("$\\textbf{b}_{40}$"))
)

temp_plot <- temp_plot + coord_cartesian(xlim = c(1, 40), ylim = c(-0.16, 0.16)) # ylim depends on confidence bands


# Print the plot
print(temp_plot)

ggsave("Results/fig4a-point-estimation.png", temp_plot, width = 6, height = 5)

############################# Placebo: Confidence Bands ############################

indx <- c(1:neval)

# Create a data frame for plotting
bound <- 40

df <- data.frame( indx = indx, y = tau.hat.biv, label = rep(c("Bivariate"), each = length(indx)))

# Build the plot
temp_plot <- ggplot() + theme_bw()
df <- df[df$indx <= bound,]

# Scatter plot for "Bivariate" and "Distance"
temp_plot <- temp_plot + geom_point(data = df,
                                    aes(x = indx, y = y, color = label, shape = label, fill = label, linetype = label))

df_ribbon <- data.frame(
  indx = indx[c(1:bound)],
  ymin = CB.lower.biv[c(1:bound)],
  ymax = CB.upper.biv[c(1:bound)],
  label = "Uniform CB" 
)

temp_plot <- temp_plot + geom_ribbon(data = df_ribbon, aes(x = indx, ymin = ymin, ymax = ymax,
                                                           color = label, shape = label, fill = label, linetype = label), alpha = 0.1)

df_errorbar <- data.frame(
  indx = indx[c(1:bound)],
  ymin = CI.lower.biv[c(1:bound)],
  ymax = CI.upper.biv[c(1:bound)],
  label = "Pointwise CI" 
)

temp_plot <- temp_plot + geom_errorbar(data = df_errorbar, aes(x = indx, ymin = ymin, ymax = ymax,
                                                               color = label, shape = label, fill = label, linetype = label))

temp_plot <- temp_plot + xlab("Cutoffs on the Boundary") + ylab("Treatment Effect")

legend_order <- c("Bivariate", "Pointwise CI", "Uniform CB")

temp_plot <- temp_plot + scale_color_manual(
  values = c("Bivariate" = "blue", "Pointwise CI" = "blue","Uniform CB" = "dodgerblue4"),
  name = NULL, breaks = legend_order
) + scale_shape_manual(
  values = c("Bivariate" = 16, "Pointwise CI" = 124,"Uniform CB" = 0),  # 16 = filled circle, 17 = triangle
  name = NULL, breaks = legend_order
) + scale_fill_manual(
  values = c("Bivariate" = NA, "Pointwise CI" = NA, "Uniform CB" = "dodgerblue4"),
  name = NULL, breaks = legend_order
) + scale_linetype_manual(
  values = c("Bivariate" = 0, "Pointwise CI" = 5,"Uniform CB" = 0),
  name = NULL, breaks = legend_order) +
  guides(
    fill = guide_legend(order = 1),  # Keep "Uniform CB" first
    color = guide_legend(order = 1),  # Ensure everything is in one group
    shape = guide_legend(order = 1),   # Align shapes with color
    linetype = guide_legend(order = 1)
  )

temp_plot <- temp_plot + geom_vline(xintercept = c(21), color = "lightgrey", size = 1, linetype = "dotted")

temp_plot <- temp_plot +
  annotate(
    "text",
    x = 19,                # Same x as the vertical line
    y = -0.16,               # Adjust this so it appears where you want
    label = "kink",
    color = "dimgrey",
    size = 4,              # Text size
    vjust = -0.5,          # Vertical justification (pulls text below this y-value)
    fontface = "bold"
  )

# Place legend inside and adjust text sizes
temp_plot <- temp_plot + theme_minimal() +
  theme(
    axis.title.x = element_text(size = 15,face = "bold"),  # X-axis label size
    axis.title.y = element_text(size = 15,face = "bold"),  # Y-axis label size
    plot.title = element_text(size = 20, hjust = 0.5),  # Title size and centering
    text=element_text(family="Times New Roman", face="bold"),
    axis.text.x = element_text(face = "bold",
                               size = 15),
    axis.text.y = element_text(face = "bold",
                               size = 12),
    legend.position = c(0.8, 1),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", colour = NA),  # White background, no border
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

temp_plot <- temp_plot + scale_x_continuous(
  breaks = c(1,5,10,15,21,25,30,35,40),
  labels = c(TeX("$\\textbf{b}_{1}$"), TeX("$\\textbf{b}_{5}$"), TeX("$\\textbf{b}_{10}$"), 
             TeX("$\\textbf{b}_{15}$"),TeX("$\\textbf{b}_{21}$"), TeX("$\\textbf{b}_{25}$"), 
             TeX("$\\textbf{b}_{30}$"), TeX("$\\textbf{b}_{35}$"), TeX("$\\textbf{b}_{40}$"))
)

temp_plot <- temp_plot + coord_cartesian(xlim = c(1, 40), ylim = c(-0.16, 0.16))

# Print the plot
print(temp_plot)

ggsave("Results/fig4b-ci-and-cb.png", temp_plot, width = 6, height = 5)

######################### Placebo: heatmap #####################################

eval.original <- eval 
scale <- TRUE
if (scale){
  eval.original$x.1 <- eval$x.1 / x.1.scale
  eval.original$x.2 <- eval$x.2 / x.2.scale
} 

# heat map for treatment effect

data.plot <- cbind(eval.original, tau.hat.biv)
colnames(data.plot) <- c("x.1", "x.2", "tau.hat")

# Function to interpolate points and color values between two consecutive points
interpolate_points <- function(df, n=ninter){
  do.call(rbind, lapply(1:(nrow(df)-1), function(i){
    xseq <- seq(df$x.1[i], df$x.1[i+1], length.out = n+2)[2:(n+1)]
    yseq <- seq(df$x.2[i], df$x.2[i+1], length.out = n+2)[2:(n+1)]
    colorseq <- seq(df$tau.hat[i], df$tau.hat[i+1], length.out = n+2)[2:(n+1)]
    data.frame(x.1 = xseq, x.2 = yseq, tau.hat = colorseq)
  }))
}

# Generating interpolated points
ninter <- 10
interpolated_data <- interpolate_points(data.plot)

# Plotting

heat_wd <- 6.5
heatcol_low <- "blue"
heatcol_high <- "red"
heat_lab <- TeX("$\\tau(\\textbf{b})$")
heat_title <- NULL
xlabel <- "Saber11"
ylabel <- "Sisben"

augmented_data <-rbind(data.plot, interpolated_data)
ord <- order(augmented_data[,1], augmented_data[,2])
augmented_data <- augmented_data[ord,]

if (is.null(heat_title)) heat_title <- "Heat Map"
plot_heat <- ggplot(data.plot, aes(x=x.1, y=x.2)) +
  geom_segment(data = augmented_data,
               aes(xend = lead(x.1, order_by=x.1),
                   yend = lead(x.2, order_by=x.2),
                   color = tau.hat),
               size = heat_wd, lineend = "round") +
  scale_color_gradient(low=heatcol_low, high=heatcol_high) +
  labs(color = heat_lab) +
  xlab(xlabel) +
  ylab(ylabel)

plot_heat <- plot_heat + geom_text(data=data.plot, aes(x=x.1, y=x.2, label=sprintf("%02d", 1:nrow(data.plot))), color="black", size=3)

plot_heat <- plot_heat + theme_minimal() +
  theme(
    axis.title.x = element_text(size = 15,face = "bold"),  # X-axis label size
    axis.title.y = element_text(size = 15,face = "bold"),  # Y-axis label size
    plot.title = element_text(size = 20, hjust = 0.5),  # Title size and centering
    text=element_text(family="Times New Roman", face="bold"),
    axis.text.x = element_text(face = "bold",
                               size = 15),
    axis.text.y = element_text(face = "bold",
                               size = 12),
    legend.position = c(0.8, 1),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", colour = NA),  # White background, no border
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + coord_cartesian(xlim = c(-10, 55), ylim = c(-10, 55))

print(plot_heat)

ggsave("Results/heat-spp-placebo.png", plot_heat, width = 6, height = 5)

# heatmap for p-value

data.plot$p.value <- placebo.rd2d$pvalues
data.plot$p.sig <- cut(data.plot$p.value,
                       breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                       labels = c("p < 0.001", "0.001 ≤ p < 0.01", "0.01 ≤ p < 0.05", "0.05 ≤ p < 0.1", "p ≥ 0.1"))
sig_colors <- c("p < 0.001" = "#d73027",       # red
                "0.001 ≤ p < 0.01" = "#fc8d59", # orange
                "0.01 ≤ p < 0.05" = "#fee08b",  # yellow
                "0.05 ≤ p < 0.1" = "#d9ef8b",   # light green
                "p ≥ 0.1" = "#91cf60")          # green
library(ggplot2)

plot_heat_pvalue <- ggplot(data.plot, aes(x = x.1, y = x.2, fill = p.sig)) +
  geom_tile(color = "white", show.legend = TRUE) +
  scale_fill_manual(values = sig_colors, name = "P-value", drop = FALSE) +
  labs(x = "Saber11", y = "Sisben") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    plot.title = element_text(size = 20, hjust = 0.5),
    text = element_text(family = "Times New Roman", face = "bold"),
    axis.text.x = element_text(face = "bold", size = 15),
    axis.text.y = element_text(face = "bold", size = 12),
    legend.position = c(0.8, 1),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

plot_heat_pvalue <- plot_heat_pvalue +  coord_cartesian(xlim = c(-10, 55), ylim = c(-10, 55)) 
plot_heat_pvalue <- plot_heat_pvalue + geom_text(data=data.plot, aes(x=x.1, y=x.2, label=sprintf("%02d", 1:nrow(data.plot))), color="black", size=3)

print(plot_heat_pvalue)

ggsave("Results/heat-pvalue-placebo.png", plot_heat_pvalue, width = 6, height = 5)

########################### conditional mean plot ##############################

library(ggplot2)
library(latex2exp)
library(extrafont)
library(grid) # For custom drawing

# Load fonts into R session
loadfonts(device = "pdf")  # Use "win" for Windows, "pdf" for PDF devices

# Define the piecewise function
f <- function(x) {
  ifelse(x <= 3/4, 
         2/pi * x, 
         (x + 3/4)/(pi - acos(0.75 / x))
  )
}

point_x <- c(0.75)

point_y <- f(point_x)
point_data <- data.frame(x = point_x, y = point_y)
point_data$label <- "0.75"

# Create a sequence of x values from 0 to 1
x_values <- seq(0, 1, length.out = 1000)

# Calculate y values using the piecewise function
y_values <- sapply(x_values, f)

# Create a data frame for plotting
data <- data.frame(x = x_values, y = y_values)


# Plot the function using ggplot2
# Plot the function using ggplot2
plot <- ggplot(data, aes(x = x, y = y)) +
  geom_line(color = "dimgrey", size = 1) +
  geom_point(data = point_data, 
             aes(x = x, y = y),
             color = "blue", 
             size = 3) +
  geom_vline(xintercept = 0.75, linetype = "dashed", color = "lightgrey") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(family = "Times", hjust = 0.5, size = 20, face = "bold"),
    axis.text.x = element_text(family = "Times", size = 15, face = "bold"),
    axis.text.y = element_text(family = "Times", size = 12, face = "bold"),
    plot.margin = margin(20, 20, 20, 20),
    panel.grid.major = element_line(color = "grey80", size = 0.5),
    panel.grid.minor = element_line(color = "grey90", size = 0.25),
    axis.line.x = element_line(color = "black", size = 0.5, 
                               arrow = grid::arrow(type = "closed", length = unit(0.15, "inches"))),
    axis.line.y = element_line(color = "black", size = 0.5, 
                               arrow = grid::arrow(type = "closed", length = unit(0.15, "inches")))
  ) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 0.75)
  ) +
  labs(x = TeX("$r$"), y = TeX("$\\theta_{1,(3/4,0)}(r)$"))
print(plot)
# Save the plot as a PDF
ggsave(sprintf("Results/conditional_mean_in_r_sa.png",bound), plot, width = 6, height = 5)


