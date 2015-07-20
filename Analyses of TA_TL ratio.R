# Script for meta-analysis of the relationship between leukocyte telomere length and hippocampus volume
# Gustav Nilsonne 2015-07-07

# REQUIRE PACKAGES

require(RCurl) # To read data from GitHub
require(psych)

# VERIFY DATA EXTRACTION

# Read data
JacobsL <- read.csv(text = getURL("https://raw.githubusercontent.com/GNilsonne/hippocampus-volume-telomere-length/master/JacobsL.csv", ssl.verifypeer = F), header = T)
JacobsR <- read.csv(text = getURL("https://raw.githubusercontent.com/GNilsonne/hippocampus-volume-telomere-length/master/JacobsR.csv", ssl.verifypeer = F), header = T)
Wolkowitz_TA_Con <- read.csv(text = getURL("https://raw.githubusercontent.com/GNilsonne/hippocampus-volume-telomere-length/master/Wolkowitz_TA_Con.csv", ssl.verifypeer = F), header = T)
Wolkowitz_TA_MDD <- read.csv(text = getURL("https://raw.githubusercontent.com/GNilsonne/hippocampus-volume-telomere-length/master/Wolkowitz_TA_MDD.csv", ssl.verifypeer = F), header = T)
Wolkowitz_TL_Con <- read.csv(text = getURL("https://raw.githubusercontent.com/GNilsonne/hippocampus-volume-telomere-length/master/Wolkowitz_TL_Con.csv", ssl.verifypeer = F), header = T)
Wolkowitz_TL_MDD <- read.csv(text = getURL("https://raw.githubusercontent.com/GNilsonne/hippocampus-volume-telomere-length/master/Wolkowitz_TL_MDD.csv", ssl.verifypeer = F), header = T)
Wolkowitz_merge_Con <- read.csv(text = getURL("https://raw.githubusercontent.com/GNilsonne/hippocampus-volume-telomere-length/master/Wolkowitz_merge_Con.csv", ssl.verifypeer = F), header = T)
Wolkowitz_merge_MDD <- read.csv(text = getURL("https://raw.githubusercontent.com/GNilsonne/hippocampus-volume-telomere-length/master/Wolkowitz_merge_MDD.csv", ssl.verifypeer = F), header = T)

# Make plots to inspect extracted data vs original plots
# Statistical results for legends are from the published papers and from the sections below
plot(TL ~ hippvol1, data = JacobsL, frame.plot = F, pch = 16, xlim = c(2.1, 3.9), ylim = c(4500, 7000),xaxs="i", yaxs="i", xaxt = "n", xlab = "Hippocampus volume", ylab = "Telomere length")
grid(nx = NA, ny = NULL, lty = 1)
abline(lm(TL ~ hippvol1, data = JacobsL))
points(TL ~ hippvol1, data = JacobsL, pch = 16)
legend("bottomright", legend = c(expression('r'[adj]*' = 0.40'), expression('r'[unadj]*' = 0.48'), expression('p'[adj]*' = 0.02'), expression('p'[unadj]*' = 0.01')))
axis(1, at = c(2.4, 2.7, 3.0, 3.3, 3.6))

plot(ratio ~ hippvol1, data = JacobsL, frame.plot = F, pch = 16, xlim = c(2.2, 3.7), ylim = c(0, 25),xaxs="i", yaxs="i", xaxt = "n", xlab = "Hippocampus volume", ylab = "Activity/length ratio")
grid(nx = NA, ny = NULL, lty = 1)
abline(lm(ratio ~ hippvol1, data = JacobsL))
points(ratio ~ hippvol1, data = JacobsL, pch = 16)
legend("topright", legend = c(expression('r'[adj]*' = -0.51'), expression('r'[unadj]*' = -0.37'), expression('p'[adj]*' = 0.003'), expression('p'[unadj]*' = 0.05')))
axis(1, at = c(2.3, 2.6, 2.9, 3.2, 3.5))

plot(TL ~ hippvol1, data = JacobsR, frame.plot = F, pch = 16, xlim = c(2.1, 3.9), ylim = c(4500, 7000),xaxs="i", yaxs="i", xaxt = "n", xlab = "Hippocampus volume", ylab = "Telomere length")
grid(nx = NA, ny = NULL, lty = 1)
abline(lm(TL ~ hippvol1, data = JacobsR))
points(TL ~ hippvol1, data = JacobsR, pch = 16)
legend("bottomright", legend = c(expression('r'[adj]*' = 0.40'), expression('r'[unadj]*' = 0.56'), expression('p'[adj]*' = 0.002'), expression('p'[unadj]*' = 0.002')))
axis(1, at = c(2.4, 2.7, 3.0, 3.3, 3.6))

plot(ratio ~ hippvol1, data = JacobsR, frame.plot = F, pch = 16, xlim = c(2.2, 3.7), ylim = c(0, 25),xaxs="i", yaxs="i", xaxt = "n", xlab = "Hippocampus volume", ylab = "Activity/length ratio")
grid(nx = NA, ny = NULL, lty = 1)
abline(lm(ratio ~ hippvol1, data = JacobsR))
points(ratio ~ hippvol1, data = JacobsR, pch = 16)
legend("topright", legend = c(expression('r'[adj]*' = -0.42'), expression('r'[unadj]*' = -0.35'), expression('p'[adj]*' = 0.001'), expression('p'[unadj]*' = 0.07')))
axis(1, at = c(2.3, 2.6, 2.9, 3.2, 3.5))

plot(TL ~ telomerase, data = JacobsL, frame.plot = F, xlim = c(5000, 120000), ylim = c(4000, 7000), xlab = "Telomerase activity", ylab = "Telomere length", main = "Comparison of estimates from left and right sides")
points(TL ~ telomerase, data = JacobsR, col = "red")
legend("topright", col = c("black", "red"), pch = 1, legend = c("left", "right"))

plot(hippvol ~ TA_z, data = Wolkowitz_TA_Con, pch = 22, xlim = c(-2.5, 2), ylim = c(3, 6.4), yaxs="i", xlab = "Telomerase activity", ylab = "Hippocampus volume", xaxt = "n", yaxt = "n")
axis(1, at = c(-2.5, -1.5, -0.5, 0.5, 1.5))
axis(2, at = c(3, 4, 5, 6))
abline(lm(hippvol ~ TA_z, data = Wolkowitz_TA_Con), lty = 2, col = "blue")
points(hippvol ~ TA_z, data = Wolkowitz_TA_MDD, pch = 16)
abline(lm(hippvol ~ TA_z, data = Wolkowitz_TA_MDD))

plot(hippvol ~ TL_z, data = Wolkowitz_TL_Con, pch = 22, xlim = c(-3.3, 3.3), ylim = c(3, 6.4), xaxs="i", yaxs="i", xlab = "Telomere length", ylab = "Hippocampus volume")
abline(lm(hippvol ~ TL_z, data = Wolkowitz_TL_Con), lty = 2)
points(hippvol ~ TL_z, data = Wolkowitz_TL_MDD, pch = 16)
abline(lm(hippvol ~ TL_z, data = Wolkowitz_TL_MDD))

# Verify matching of data points by hippocampus volume
plot(hippvol2 ~ hippvol1, data = JacobsL, frame.plot = F)
cor.test(JacobsL$hippvol2, JacobsL$hippvol1)

plot(hippvol2 ~ hippvol1, data = JacobsR, frame.plot = F)
cor.test(JacobsR$hippvol2, JacobsL$hippvol1)

plot(hippvol ~ hippvol.1, data = Wolkowitz_merge_Con, frame.plot = F)
cor.test(Wolkowitz_merge_Con$hippvol, Wolkowitz_merge_Con$hippvol.1)

plot(hippvol ~ hippvol.1, data = Wolkowitz_merge_MDD, frame.plot = F)
cor.test(Wolkowitz_merge_MDD$hippvol, Wolkowitz_merge_MDD$hippvol.1)

# ANALYSE TELOMERASE/TELOMERE LENGTH RATIO

# Attempt to reconstruct results from plot in paper by Jacobs
plot(telomerase ~ hippvol1, data = JacobsL, frame.plot = F)
cor.test(JacobsL$telomerase, JacobsL$hippvol1)
cor.test(JacobsL$ratio, JacobsL$hippvol1)
cor.test(JacobsL$TL, JacobsL$hippvol1)

plot(telomerase ~ hippvol1, data = JacobsR, frame.plot = F)
cor.test(JacobsR$telomerase, JacobsR$hippvol1)
cor.test(JacobsR$ratio, JacobsR$hippvol1)
cor.test(JacobsR$TL, JacobsR$hippvol1)

# Find correlations between telomerase activity, telomere length, and their ratio
cor.test(JacobsL$telomerase, JacobsL$TL)
cor.test(JacobsL$telomerase, JacobsL$ratio)
cor.test(JacobsL$TL, JacobsL$ratio)

cor.test(JacobsR$telomerase, JacobsR$TL)
cor.test(JacobsR$telomerase, JacobsR$ratio)
cor.test(JacobsR$TL, JacobsR$ratio)

cor.test(Wolkowitz_merge_MDD$hippvol, Wolkowitz_merge_MDD$TA)
cor.test(Wolkowitz_merge_MDD$hippvol, Wolkowitz_merge_MDD$ratio)
cor.test(Wolkowitz_merge_MDD$TA, Wolkowitz_merge_MDD$ratio)
cor.test(Wolkowitz_merge_MDD$hippvol, Wolkowitz_merge_MDD$TL)
cor.test(Wolkowitz_merge_MDD$TL, Wolkowitz_merge_MDD$ratio)

cor.test(Wolkowitz_merge_Con$hippvol, Wolkowitz_merge_Con$TA)
cor.test(Wolkowitz_merge_Con$hippvol, Wolkowitz_merge_Con$ratio)
cor.test(Wolkowitz_merge_Con$TA, Wolkowitz_merge_Con$ratio)
cor.test(Wolkowitz_merge_Con$hippvol, Wolkowitz_merge_Con$TL)
cor.test(Wolkowitz_merge_Con$TL, Wolkowitz_merge_Con$ratio)

# Compare r:s between predictors
r.test(n = 28, r12 = -0.33, r13 = -0.37, r23 = 0.99) # Jacobs, left hippocampus. r12 is telomerase-volume, r13 is ratio-volume, and r23 is telomerase-ratio
r.test(n = 28, r12 = 0.48, r13 = 0.37, r23 = 0.32) # Jacobs, left hippocampus. r12 is TL-volume, r13 is ratio-volume, and r23 is TL-ratio
r.test(n = 28, r12 = -0.28, r13 = -0.35, r23 = 0.99) # Jacobs, right hippocampus. r12 is telomerase-volume, r13 is ratio-volume, and r23 is telomerase-ratio
r.test(n = 28, r12 = 0.56, r13 = 0.35, r23 = 0.27) # Jacobs, right hippocampus. r12 is TL-volume, r13 is ratio-volume, and r23 is TL-ratio

r.test(n = 19, r12 = 0.58, r13 = 0.41, r23 = 0.29) # Wolkowitz MDD. r12 is telomerase-volume, r13 is ratio-volume, and r23 is telomerase-ratio
r.test(n = 19, r12 = 0.07, r13 = 0.41, r23 = 0.16) # Wolkowitz MDD. r12 is TL-volume, r13 is ratio-volume, and r23 is TL-ratio
r.test(n = 16, r12 = 0.18, r13 = 0, r23 = 0.39) # Wolkowitz controls. r12 is telomerase-volume, r13 is ratio-volume, and r23 is telomerase-ratio
r.test(n = 16, r12 = 0.09, r13 = 0, r23 = 0.10) # Wolkowitz controls. r12 is TL-volume, r13 is ratio-volume, and r23 is TL-ratio

# Investigate relation between telomere length and telomerase
plot(TL_z ~ TA_z, data = Wolkowitz_merge_MDD, frame.plot = F)
cor.test(Wolkowitz_merge_MDD$TL, Wolkowitz_merge_MDD$TA)

plot(TL_z ~ TA_z, data = Wolkowitz_merge_Con, frame.plot = F)
cor.test(Wolkowitz_merge_Con$TL, Wolkowitz_merge_Con$TA)
