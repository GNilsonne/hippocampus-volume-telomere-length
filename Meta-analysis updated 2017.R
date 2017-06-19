# Script for meta-analysis of the relationship between leukocyte telomere length and hippocampus volume
# Gustav Nilsonne 2017-06-17

# This script updates a previous meta-analysis by adding more data

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
data <- data.frame(study = "Grodstein08", group.in.study = "1.in.Grodstein08", author = "Grodstein", year = 2008, ri = r, ni = n, covariates = "age, educational attainment")

# Wikgren 2012
# We have r and n for both groups
data <- rbind(data, data.frame(study = "Wikgren12", group.in.study = "1.in.Wikgren12", author = "Wikgren_e3e3", year = 2012, ri = -0.519, ni = 29, covariates = "age, body size"))
data <- rbind(data, data.frame(study = "Wikgren12", group.in.study = "2.in.Wikgren12", author = "Wikgren_e4", year= 2012, ri = 0.134, ni = 28, covariates = "age, body size"))

# King 2014
# We have r2 and n
data <- rbind(data, data.frame(study = "King14", group.in.study = "1.in.King14", author = "King", year = 2014, ri = sqrt(0.0091), ni = 1960, covariates = "age, sex, race/ethnicity"))

# Wolkowitz 2015
# We have r and n for both groups
data <- rbind(data, data.frame(study = "Wolkowitz15", group.in.study = "1.in.Wolkowitz15", author = "Wolkowitz_MDD", year = 2015, ri = 0.05, ni = 19, covariates = "age, sex"))
data <- rbind(data, data.frame(study = "Wolkowitz15", group.in.study = "2.in.Wolkowitz15", author = "Wolkowitz_control", year = 2015, ri = 0.35, ni = 17, covariates = "age, sex"))

# Jacobs 2015
# We have r2 and n
data <- rbind(data, data.frame(study = "Jacobs15", group.in.study = "1.in.Jacobs15", author = "Jacobs", year = 2015, ri = sqrt(0.16), ni = 28, covariates = "age, education, BMI"))

# Sleepy Brain 2017
# Effects calculated in separate script
data <- rbind(data, data.frame(study = "SleepyBrain17", group.in.study = "1.in.SleepyBrain17", author = "SleepyBrain17_younger", year = 2017, ri = 0.0887, ni = 44, covariates = "age, ICV"))
data <- rbind(data, data.frame(study = "SleepyBrain17", group.in.study = "2.in.SleepyBrain17", author = "SleepyBrain17_older", year = 2017, ri = 0.0752, ni = 33, covariates = "age, ICV"))

# WRITE TABLE TO LATEX
print(xtable(data))

# PERFORM META-ANALYSIS

# Perform analyses
res <- rma(ri = ri, ni = ni, measure = "COR", slab=paste(author, year, sep=" "), data = data)
res
res_old <- rma(ri = ri, ni = ni, measure = "COR", slab=paste(author, year, sep=" "), data = data[data$study != "SleepyBrain17", ])
res_old

# Make forest plots
forest(res)
forest(res_old)

# Make trim-and-fill analysis and show funnel plot
taf <- trimfill(res)
taf
funnel(taf)
taf_old <- trimfill(res_old)
taf_old
funnel(taf_old)

# Make q-q- plot
qqnorm(res)

# PERFORM MULTILEVEL ANALYSIS
# Since the dataset contains effects nested within studies, multilevel modelling is more appropriate in theory. 
# Because the rma.mv function yields an object for which trim-and-fill analyses cannot be performed, we add this just to verify that this alternate specification gives similar results to the main analysis above.
data <- data.frame(data, escalc(ri = data$ri, ni = data$ni, measure = "COR"))
res.ml <- rma.mv(yi, vi, random = list(~ 1 | study, ~ 1 | group.in.study), slab=paste(author, year, sep=" "), data = data)
summary.rma(res.ml)
forest(res.ml)
profile(res.ml, sigma2=1)
profile(res.ml, sigma2=2)

# ESTIMATE POWER
pwr.r.test(n = 33, r = 0.12)
pwr.r.test(n = 44, r = 0.12)

# MAKE FIGURES FOR PUBLICATION

pdf("figure_newforest.pdf") # Note: First polygon should be removed before publication using Inkscape or similar
forest(res.ml, alim = c(-0.8, 0.8), digits = c(2, 1), xlab = "r", ylim=c(-6,12), cex = 1.2, steps = 3)
addpoly(res_old, row=-2, mlab="RE without SB17", cex = 1.2)
addpoly(taf_old, row=-3, mlab="RE without SB17, TAF adjusted", cex = 1.2)
addpoly(res.ml, row=-4, mlab="RE", cex = 1.2)
addpoly(taf, row=-5, mlab="RE, TAF adjusted", cex = 1.2)
dev.off()