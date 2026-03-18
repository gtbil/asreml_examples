# set file location -------------------------------------------------------
here::i_am("code/01_RepeatedMeasuresByPlot.R")

# set constants -----------------------------------------------------------
n.reps <- 10
n.genos <- 30
n.years <- 5
n.plot <- n.reps * n.genos * n.years

var.rep <- 1
var.year <- 1
var.geno <- 2
var.resid2 <- 0.1

# generate the random components ------------------------------------------
set.seed(123)
var.resids <- truncnorm::rtruncnorm(n = n.years, 
                                    a = 0.5, b = 1.5, mean = 0.5, sd = 0.2)
D <- diag(sqrt(var.resids))
R <- clusterGeneration::rcorrmatrix(n.years)

# simulate a strictly positive covariance matrix
# L <- matrix(abs(rnorm(n.years^2)), n.years, n.years)
# R_raw <- tcrossprod(L)
# R <- cov2cor(R_raw)

vcov.resid <- D %*% R %*% D

eigen(R)$values
eigen(vcov.resid)$values

# # ensure pos definitite
# vcov.resid <- Matrix::nearPD(vcov.resid)$mat |>
#   as.matrix()

# simulate data -----------------------------------------------------------
sim.data <- expand.grid(year = 1:n.years, rep = 1:n.reps, geno = 1:n.genos) |>
  data.frame() |>
  dplyr::arrange(year, rep, geno)

sim.data$plot <- with(sim.data, paste0("rep", 
                                       sprintf("%02d", rep), "_geno", 
                                       sprintf("%03d", geno)) |>
                        as.factor() |>
                        as.integer())

# generate random effects
effect.rep <- rnorm(n.reps, mean = 0, sd = sqrt(var.rep))
effect.geno <- rnorm(n.genos, mean = 0, sd = sqrt(var.geno))
effect.year <- rnorm(n.years, mean = 0, sd = sqrt(var.year))
effect.plot <- MASS::mvrnorm(n = n.reps * n.genos, 
                             mu = rep(0, n.years), 
                             Sigma = vcov.resid) |>
  c()

sim.effects <- data.frame(rep = effect.rep[sim.data$rep],
                          geno = effect.geno[sim.data$geno],
                          year = effect.year[sim.data$year],
                          plot = effect.plot[sim.data$plot])

y <- rowSums(sim.effects) + rnorm(n.plot, mean = 0, sd = sqrt(var.resid2))

data <- cbind(lapply(sim.data, factor) |> data.frame(), y = y)

# fit some simple models --------------------------------------------------
model1 <- lme4::lmer(y ~ (1|rep) + (1|geno) + year, data = data)

plot(lme4::ranef(model1)$geno[,1], effect.geno) +
  abline(a = 0, b = 1, col = "red")

# fit an asreml model that correctly matches this simulation
mod <- asreml::asreml(
  fixed = y ~ 1,
  random = ~ geno + year + rep,
  residual = ~ plot:corgh(year),
  data = data |>
    dplyr::arrange(plot, year)
)

summary(mod)$varcomp

# geno is not separable here... because geno and year are confounded with plot