# rd2d: post analysis of Monte Carlo
# Authors: M. D. Cattaneo, R. Titiunik, R. R. Yu
# Last update: 
# rm(list=ls(all=TRUE))
# library(binsreg); library(ggplot2)

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

x <- read.csv("CIT_2023_CUP_multiscore-nongeo.csv")
x$X <- NULL
colnames(x) <- c("x.1", "x.2","y","d")
na.ok <- complete.cases(x$x.1) & complete.cases(x$x.2)
x <- x[na.ok,]

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

# Polynomial Fits to the Data

# Define the 2nd order polynomial
# formula <- y ~ x.1 + x.2 + I(x.1^2) + I(x.1 * x.2) + I(x.2^2)
formula <- y ~ x.1 + x.2

# Fit models separately for d = 0 and d = 1
model_d0 <- lm(formula, data = x[x$d == 0, ])
model_d1 <- lm(formula, data = x[x$d == 1, ])

# Calculate true treatment effects
true_effect <- data.frame(eval) # create copy of eval
true_effect$mu.0 <- predict(model_d0, newdata = true_effect)
true_effect$mu.1 <- predict(model_d1, newdata = true_effect)
true_effect$tau <- (true_effect$mu.1 - true_effect$mu.0) * 2

m <- 1000

# subset <- c(1,5,10,15,21,25,30,35,40)
# 
# neval_subset <- length(subset)

subset <- c(1:neval)
neval_subset <- neval

result <- list("rd2d" = matrix(NA, nrow = neval_subset, ncol = 12),
               "dist_kinkon" = matrix(NA, nrow = neval_subset, ncol = 12),
               "dist_kinkoff" = matrix(NA, nrow = neval_subset, ncol = 12))

for (method in c("rd2d", "dist_kinkon", "dist_kinkoff")){
  
  h <- rep(0,neval_subset)
  tau.hat <- rep(0,neval_subset)
  ci.lower <- rep(0,neval_subset)
  ci.upper <- rep(0,neval_subset)
  cb.lower <- rep(0,neval_subset)
  cb.upper <- rep(0,neval_subset)
  mse <- rep(0,neval_subset)
  ec.ptw <- rep(0,neval_subset)
  il.ptw <- rep(0,neval_subset)
  ec.unif <- 0
  il.unif <- rep(0,neval_subset)
  
  tau.hat.full <- matrix(0, m, neval_subset)
  
  shift <- 1
  
  for (i in c(1:m)){
    dat.method <- as.matrix(read.csv(sprintf("rd2d_monte/%s_sim%d_m1000_sd01_n20000_linear.csv", method, i))) 
    h <- h + dat.method[subset,1]
    tau.hat <- tau.hat + dat.method[subset,1 + shift]
    tau.hat.full[i,] <- dat.method[subset,1 + shift]
    ci.lower <- ci.lower + dat.method[subset,2 + shift]
    ci.upper <- ci.upper + dat.method[subset,3 + shift]
    cb.lower <- cb.lower + dat.method[subset,4 + shift]
    cb.upper <- cb.upper + dat.method[subset,5 + shift]
    
    mse <- mse + (dat.method[subset,1 + shift] - true_effect$tau[subset])^2
    ec.ptw <- ec.ptw + as.integer((dat.method[subset,2 + shift] <= true_effect$tau[subset]) & (true_effect$tau[subset] <= dat.method[subset,3 + shift]))
    il.ptw <- il.ptw + dat.method[subset,3 + shift] - dat.method[subset,2 + shift]
    ec.unif <- ec.unif + as.integer(sum((dat.method[subset,4 + shift] <= true_effect$tau[subset]) & (true_effect$tau[subset] <= dat.method[subset,5 + shift])) == neval_subset)
    il.unif <- il.unif + dat.method[subset,5 + shift] - dat.method[subset,4 + shift]
  }
  
  h <- h / m
  tau.hat <- tau.hat / m
  ci.lower <- ci.lower / m
  ci.upper <- ci.upper / m
  cb.lower <- cb.lower / m
  cb.upper <- cb.upper / m
  bias <- tau.hat - true_effect$tau[subset]
  mse <- mse / m
  rmse <- sqrt(mse)
  se <- apply(tau.hat.full, 2, sd)
  ec.ptw <- ec.ptw / m
  il.ptw <- il.ptw / m
  ec.unif <- rep(ec.unif/m, neval_subset)
  il.unif <- il.unif / m
  
  # tau.full.df <- as.data.frame(tau.hat.q)
  # se <- apply(tau.full.df, 2, sd)
  
  result[[method]]  <- cbind(h, tau.hat, ci.lower, ci.upper, cb.lower, cb.upper,
                             bias, se, rmse, il.ptw, ec.ptw, il.unif, ec.unif)
}

subset <- c(1, 5, 10, 15, 21, 25, 30, 35, 40)

# Extract the desired subset and columns
data_subset <- result$dist_kinkoff[subset, -c(2:6)]
ec.unif <- data_subset[,]

num_cols <- ncol(data_subset)
digits_vector <- c(0, rep(3, num_cols))

# Create the xtable object with scientific notation (2 significant digits)
xt <- xtable(data_subset, digits = digits_vector)

# Print the LaTeX table with math-style exponents
print(xt, type = "latex")

indx <- c(1:neval)

# Only plot for bound = 40 and config = 4
bound <- 40
config <- 4

# Extract estimates
tau.hat.biv <- result$rd2d[, 2]
tau.hat.dist.kinkoff <- result$dist_kinkoff[, 2]
tau.hat.dist.kinkon <- result$dist_kinkon[, 2]

# Prepare data: swap order of kinkon and kinkoff
df <- data.frame(
  indx = rep(indx, 4),
  y = c(tau.hat.biv, tau.hat.dist.kinkon, tau.hat.dist.kinkoff, true_effect$tau),
  label = rep(c("Bivariate", "Dist (kink)", "Dist (no kink)", "Population"), each = length(indx))
)

df <- df[df$indx <= bound, ]

# Initialize plot
temp_plot <- ggplot() + theme_bw()

# Add population line
temp_plot <- temp_plot + geom_line(
  data = df[df$label == "Population", ],
  aes(x = indx, y = y, color = label, shape = label, linetype = label),
  size = 0.5, show.legend = TRUE
)

# Add Dist (no kink) points
temp_plot <- temp_plot + geom_point(
  data = df[df$label == "Dist (no kink)", ],
  aes(x = indx, y = y, color = label, shape = label, linetype = label)
)

# Add Dist (kink) points
temp_plot <- temp_plot + geom_point(
  data = df[df$label == "Dist (kink)", ],
  aes(x = indx, y = y, color = label, shape = label, linetype = label)
)

# Add Bivariate points
temp_plot <- temp_plot + geom_point(
  data = df[df$label == "Bivariate", ],
  aes(x = indx, y = y, color = label, shape = label, linetype = label)
)

# Axis labels
temp_plot <- temp_plot + xlab("Cutoffs on the Boundary") + ylab("Treatment Effect")

# Legend customization (updated order)
legend_order <- c("Population", "Dist (no kink)", "Dist (kink)", "Bivariate")

temp_plot <- temp_plot +
  scale_color_manual(
    values = c("Bivariate" = "blue", "Dist (no kink)" = "red", "Dist (kink)" = "grey", "Population" = "black"),
    name = NULL, breaks = legend_order
  ) +
  scale_shape_manual(
    values = c("Bivariate" = 16, "Dist (no kink)" = 4, "Dist (kink)" = 1, "Population" = NA),
    name = NULL, breaks = legend_order
  ) +
  scale_linetype_manual(
    values = c("Bivariate" = 0, "Dist (no kink)" = 0, "Dist (kink)" = 0, "Population" = 1),
    name = NULL, breaks = legend_order
  ) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 1),
    linetype = guide_legend(order = 1)
  )

# Vertical line and annotation
temp_plot <- temp_plot +
  geom_vline(xintercept = 21, color = "lightgrey", size = 1, linetype = "dotted") +
  annotate(
    "text",
    x = 19, y = 0.69,
    label = "kink",
    color = "dimgrey",
    size = 4,
    vjust = -0.5,
    fontface = "bold"
  )

# Theme
temp_plot <- temp_plot + theme_minimal() +
  theme(
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    plot.title = element_text(size = 20, hjust = 0.5),
    text = element_text(family = "Times New Roman", face = "bold"),
    axis.text.x = element_text(face = "bold", size = 15),
    axis.text.y = element_text(face = "bold", size = 12),
    legend.position = c(0.8, 1),
    legend.justification = c(0, 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# x-axis labels
temp_plot <- temp_plot + scale_x_continuous(
  breaks = c(1, 5, 10, 15, 21, 25, 30, 35, 40),
  labels = c(
    TeX("$\\textbf{b}_{1}$"), TeX("$\\textbf{b}_{5}$"), TeX("$\\textbf{b}_{10}$"),
    TeX("$\\textbf{b}_{15}$"), TeX("$\\textbf{b}_{21}$"), TeX("$\\textbf{b}_{25}$"),
    TeX("$\\textbf{b}_{30}$"), TeX("$\\textbf{b}_{35}$"), TeX("$\\textbf{b}_{40}$")
  )
)

# y-axis no labels
temp_plot <- temp_plot + scale_y_continuous(breaks = NULL)

# Adjust coordinate limits
temp_plot <- temp_plot + coord_cartesian(xlim = c(1, 40), ylim = c(0.69,0.83))

# Print and save the plot
print(temp_plot)

ggsave("Results/bias_optimal_bw_simulation.png", temp_plot, width = 6, height = 5)

############################# Main paper: Bias #################################

x <- read.csv("CIT_2023_CUP_multiscore-nongeo.csv")
x$X <- NULL
colnames(x) <- c("x.1", "x.2","y","d")
na.ok <- complete.cases(x$x.1) & complete.cases(x$x.2)
x <- x[na.ok,]

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

# Polynomial Fits to the Data

# Define the 2nd order polynomial
formula <- y ~ x.1 + x.2 + I(x.1^2) + I(x.1 * x.2) + I(x.2^2)
# formula <- y ~ x.1 + x.2

# Fit models separately for d = 0 and d = 1
model_d0 <- lm(formula, data = x[x$d == 0, ])
model_d1 <- lm(formula, data = x[x$d == 1, ])

# Calculate true treatment effects
true_effect <- data.frame(eval) # create copy of eval
true_effect$mu.0 <- predict(model_d0, newdata = true_effect)
true_effect$mu.1 <- predict(model_d1, newdata = true_effect)
true_effect$tau <- (true_effect$mu.1 - true_effect$mu.0) * 2

m <- 1000

# subset <- c(1,5,10,15,21,25,30,35,40)
# 
# neval_subset <- length(subset)

subset <- c(1:neval)
neval_subset <- neval

result <- list("rd2d" = matrix(NA, nrow = neval_subset, ncol = 12),
               "dist_kinkon" = matrix(NA, nrow = neval_subset, ncol = 12),
               "dist_kinkoff" = matrix(NA, nrow = neval_subset, ncol = 12))

for (method in c("rd2d", "dist_kinkon", "dist_kinkoff")){
  
  h <- rep(0,neval_subset)
  tau.hat <- rep(0,neval_subset)
  ci.lower <- rep(0,neval_subset)
  ci.upper <- rep(0,neval_subset)
  cb.lower <- rep(0,neval_subset)
  cb.upper <- rep(0,neval_subset)
  mse <- rep(0,neval_subset)
  ec.ptw <- rep(0,neval_subset)
  il.ptw <- rep(0,neval_subset)
  ec.unif <- 0
  il.unif <- rep(0,neval_subset)
  
  tau.hat.full <- matrix(0, m, neval_subset)
  
  shift <- 1
  
  for (i in c(1:m)){
    dat.method <- as.matrix(read.csv(sprintf("rd2d_monte/%s_sim%d_m1000_sd01_n20000_hc3.csv", method, i))) 
    h <- h + dat.method[subset,1]
    tau.hat <- tau.hat + dat.method[subset,1 + shift]
    tau.hat.full[i,] <- dat.method[subset,1 + shift]
    ci.lower <- ci.lower + dat.method[subset,2 + shift]
    ci.upper <- ci.upper + dat.method[subset,3 + shift]
    cb.lower <- cb.lower + dat.method[subset,4 + shift]
    cb.upper <- cb.upper + dat.method[subset,5 + shift]
    
    mse <- mse + (dat.method[subset,1 + shift] - true_effect$tau[subset])^2
    ec.ptw <- ec.ptw + as.integer((dat.method[subset,2 + shift] <= true_effect$tau[subset]) & (true_effect$tau[subset] <= dat.method[subset,3 + shift]))
    il.ptw <- il.ptw + dat.method[subset,3 + shift] - dat.method[subset,2 + shift]
    ec.unif <- ec.unif + as.integer(sum((dat.method[subset,4 + shift] <= true_effect$tau[subset]) & (true_effect$tau[subset] <= dat.method[subset,5 + shift])) == neval_subset)
    il.unif <- il.unif + dat.method[subset,5 + shift] - dat.method[subset,4 + shift]
  }
  
  h <- h / m
  tau.hat <- tau.hat / m
  ci.lower <- ci.lower / m
  ci.upper <- ci.upper / m
  cb.lower <- cb.lower / m
  cb.upper <- cb.upper / m
  bias <- tau.hat - true_effect$tau[subset]
  mse <- mse / m
  rmse <- sqrt(mse)
  se <- apply(tau.hat.full, 2, sd)
  ec.ptw <- ec.ptw / m
  il.ptw <- il.ptw / m
  ec.unif <- rep(ec.unif/m, neval_subset)
  il.unif <- il.unif / m
  il.unif <- mean(abs(il.unif))
  il.unif <- rep(il.unif, neval_subset)
  
  # tau.full.df <- as.data.frame(tau.hat.q)
  # se <- apply(tau.full.df, 2, sd)
  
  result[[method]]  <- cbind(h, tau.hat, ci.lower, ci.upper, cb.lower, cb.upper,
                             bias, se, rmse, ec.ptw, il.ptw, ec.unif, il.unif)
}

subset <- c(1, 5, 10, 15, 21, 25, 30, 35, 40)

# Extract the desired subset and columns
data_subset <- result$rd2d[subset, -c(2:6)]
ec.unif <- data_subset[1, ncol(data_subset) - 1]
ec.unif <- as.numeric(ec.unif)
il.unif <- data_subset[1, ncol(data_subset)]
il.unif <- as.numeric(il.unif)

data_subset[, c(ncol(data_subset) - 1, ncol(data_subset))] <- NA
new_row <- rep(NA, ncol(data_subset))
new_row[length(new_row) - 1] <- ec.unif
new_row[length(new_row)] <- il.unif

data_subset <- rbind(data_subset, new_row)
rownames(data_subset) <- NULL

nrow_data <- nrow(data_subset)
left_col <- rep("", nrow_data)
# left_col[1] <- sprintf("\\multirow{%d}{*}{\\shortstack{\\code{rd2d.dist}\\\\\\code{kink = \"off\"}}}", nrow_data)
left_col[1] <- sprintf("\\multirow{%d}{*}{\\code{rd2d}}", nrow_data)
index_full <- c(1:neval)
index_col <- c(as.character(index_full[subset]),"") 
data_to_print <- data.frame(left_col, index_col, data_subset)
digits_vec <- c(0, 0, 0, rep(3, ncol(data_subset)))
tab <- xtable(data_to_print, align = c("l", "l", "l", rep("r", ncol(data_subset))),digits = digits_vec)
print(tab, include.rownames = FALSE, sanitize.text.function = identity)


indx <- c(1:neval)

# Only plot for bound = 40 and config = 4
bound <- 40
config <- 4

# Extract estimates
tau.hat.biv <- result$rd2d[, 2]
tau.hat.dist.kinkoff <- result$dist_kinkoff[, 2]
tau.hat.dist.kinkon <- result$dist_kinkon[, 2]

# Prepare data: swap order of kinkon and kinkoff
df <- data.frame(
  indx = rep(indx, 4),
  y = c(tau.hat.biv, tau.hat.dist.kinkon, tau.hat.dist.kinkoff, true_effect$tau),
  label = rep(c("Bivariate", "Dist (kink)", "Distance", "Population"), each = length(indx))
)

df <- df[df$indx <= bound, ]

# Initialize plot
temp_plot <- ggplot() + theme_bw()

# Add population line
temp_plot <- temp_plot + geom_line(
  data = df[df$label == "Population", ],
  aes(x = indx, y = y, color = label, shape = label, linetype = label),
  size = 0.5, show.legend = TRUE
)

# Add Dist (no kink) points
temp_plot <- temp_plot + geom_point(
  data = df[df$label == "Distance", ],
  aes(x = indx, y = y, color = label, shape = label, linetype = label)
)

# # Add Dist (kink) points
# temp_plot <- temp_plot + geom_point(
#   data = df[df$label == "Dist (kink)", ],
#   aes(x = indx, y = y, color = label, shape = label, linetype = label)
# )

# Add Bivariate points
temp_plot <- temp_plot + geom_point(
  data = df[df$label == "Bivariate", ],
  aes(x = indx, y = y, color = label, shape = label, linetype = label)
)

# Axis labels
temp_plot <- temp_plot + xlab("Cutoffs on the Boundary") + ylab("Treatment Effect")

# Legend customization (updated order)
legend_order <- c("Population", "Bivariate", "Distance", "Dist (kink)")

temp_plot <- temp_plot +
  scale_color_manual(
    values = c("Bivariate" = "blue", "Distance" = "red", "Dist (kink)" = "grey", "Population" = "black"),
    name = NULL, breaks = legend_order
  ) +
  scale_shape_manual(
    values = c("Bivariate" = 16, "Distance" = 4, "Dist (kink)" = 1, "Population" = NA),
    name = NULL, breaks = legend_order
  ) +
  scale_linetype_manual(
    values = c("Bivariate" = 0, "Distance" = 0, "Dist (kink)" = 0, "Population" = 1),
    name = NULL, breaks = legend_order
  ) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 1),
    linetype = guide_legend(order = 1)
  )

# Vertical line and annotation
temp_plot <- temp_plot +
  geom_vline(xintercept = 21, color = "lightgrey", size = 1, linetype = "dotted") +
  annotate(
    "text",
    x = 19, y = 0.69,
    label = "kink",
    color = "dimgrey",
    size = 4,
    vjust = -0.5,
    fontface = "bold"
  )

# Theme
temp_plot <- temp_plot + theme_minimal() +
  theme(
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    plot.title = element_text(size = 20, hjust = 0.5),
    text = element_text(family = "Times New Roman", face = "bold"),
    axis.text.x = element_text(face = "bold", size = 15),
    axis.text.y = element_text(face = "bold", size = 12),
    legend.position = c(0.8, 1),
    legend.justification = c(0, 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# x-axis labels
temp_plot <- temp_plot + scale_x_continuous(
  breaks = c(1, 5, 10, 15, 21, 25, 30, 35, 40),
  labels = c(
    TeX("$\\textbf{b}_{1}$"), TeX("$\\textbf{b}_{5}$"), TeX("$\\textbf{b}_{10}$"),
    TeX("$\\textbf{b}_{15}$"), TeX("$\\textbf{b}_{21}$"), TeX("$\\textbf{b}_{25}$"),
    TeX("$\\textbf{b}_{30}$"), TeX("$\\textbf{b}_{35}$"), TeX("$\\textbf{b}_{40}$")
  )
)

# y-axis no labels
temp_plot <- temp_plot + scale_y_continuous(breaks = NULL)

# Adjust coordinate limits
temp_plot <- temp_plot + coord_cartesian(xlim = c(1, 40), ylim = c(0.69,0.83))

# Print and save the plot
print(temp_plot)

ggsave("Results/bias_optimal_bw_simulation_main_paper.png", temp_plot, width = 6, height = 5)
