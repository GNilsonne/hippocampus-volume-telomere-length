# Script for meta-analysis of the relationship between leukocyte telomere length and hippocampus volume
# Gustav Nilsonne 2015-07-07

# REQUIRE PACKAGES

require(compute.es)
require(metafor)
require(xtable)
require(psych)
require(pwr)

# DETERMINE EFFECT SIZES FOR EACH PAPER

# Grodstein 2008
# We have n and p
n <- 26
df <- n - 1
p <- 0.038
t <- qt((1-p/2), df)
r <- sqrt(t^2/(t^2 + df))
res(r = r, n = n)
data <- data.frame(author = "Grodstein", year = 2008, ri = r, ni = n, covariates = "age, educational attainment")

# Wikgren 2012
# We have r and n for both groups
data <- rbind(data, data.frame(author = "Wikgren_e3e3", year = 2012, ri = -0.519, ni = 29, covariates = "age, body size"))
data <- rbind(data, data.frame(author = "Wikgren_e4", year= 2012, ri = 0.134, ni = 28, covariates = "age, body size"))

# King 2014
# We have r2 and n
data <- rbind(data, data.frame(author = "King", year = 2014, ri = sqrt(0.0091), ni = 1960, covariates = "age, sex, race/ethnicity"))

# Wolkowitz 2015
# We have r and n for both groups
data <- rbind(data, data.frame(author = "Wolkowitz_MDD", year = 2015, ri = 0.05, ni = 19, covariates = "age, sex"))
data <- rbind(data, data.frame(author = "Wolkowitz_control", year = 2015, ri = 0.35, ni = 17, covariates = "age, sex"))

# Jacobs 2015
# We have r2 and n
data <- rbind(data, data.frame(author = "Jacobs", year = 2015, ri = sqrt(0.16), ni = 28, covariates = "age, education, BMI"))

# WRITE TABLE TO LATEX
print(xtable(data))

# PERFORM META-ANALYSIS

# Perform analysis
res <- rma(ri = ri, ni = ni, measure = "COR", slab=paste(author, year, sep=" "), data = data)
res

# Make forest plot
forest(res)

# Make trim-and-fill analysis and show funnel plot
taf <- trimfill(res)
taf
funnel(taf)

# Make q-q- plot
qqnorm(res)

# Perform analysis again, but excluding Wiksten 2012
res2 <- rma(ri = ri, ni = ni, measure = "COR", slab=paste(author, year, sep=" "), data = data[data$year != 2012, ])
res2

# Make forest plot
forest(res2)

# Make trim-and-fill analysis and show funnel plot
taf2 <- trimfill(res2)
taf2
funnel(taf2)

# Make q-q- plot
qqnorm(res2)

# PERFORM META-REGRESSION
# Age and sex are the two identified candidate predictors

data$age_mean <- c(79.2, 62.1, 61.1, 50, 37.8, 34.9, 58.0) #Grodstein did not give age for imaging subsample, age for whole healthy subgroup substituted. King gave median age.
data$age_sd <- c(2.2, 8.5, 8.1, NA, 12.0, 9.6, 4.7)
data$age_median <- c(NA)
data$age_range <- c(NA)
data$sex_f_n <- c(26, 18, 18, 1153, 10, 17, 28)
data$sex_f_frac <- c(1, 0.62, 0.64, 0.59, 0.63, 0.59, 1)

res3 <- rma(ri = ri, ni = ni, measure = "COR", slab=paste(author, year, sep=" "), data = data, mods = age_mean)
res3

res4 <- rma(ri = ri, ni = ni, measure = "COR", slab=paste(author, year, sep=" "), data = data, mods = sex_f_frac)
res4

# ESTIMATE OBSERVED POWER
pwr.r.test(n = 30, r = 0.12)
pwr.r.test(n = 30, r = 0.23)

# MAKE FIGURES FOR PUBLICATION

pdf("figure1.pdf")
forest(res, alim = c(-0.8, 0.8), digits = c(2, 1), xlab = "r", ylim=c(-2.5,10), cex = 1.2, steps = 3)
addpoly(res2, row=-2, mlab="RE Model, Wikgren 2012 excluded", cex = 1.2)
dev.off()

pdf("figure2.pdf")
par(mfrow=c(2,2))
funnel(taf, main = "Full sample", at = c(-0.6, -0.3, 0, 0.3, 0.6), steps = 3, digits = c(1, 3))
funnel(taf2, main = "Excluding Wikgren 2012", at = c(-0.6, -0.3, 0, 0.3, 0.6), steps = 3, digits = c(1, 3))
qqnorm(res, main = "Full sample", xaxt = "n", yaxt = "n")
axis(1, at = c(-2, -1, 0, 1))
axis(2, at = c(-2, -1, 0, 1))
qqnorm(res2, main = "Excluding Wikgren 2012", xaxt = "n", yaxt = "n")
axis(1, at = c(-2, -1, 0, 1))
axis(2, at = c(-2, -1, 0, 1))
dev.off()
