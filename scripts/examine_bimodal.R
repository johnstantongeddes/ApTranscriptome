# Script to 

# Transcript with maximum standard deviation of expression
tmax <- A22_bim_sd[which(A22_bim_sd$exp_sd == max(A22_bim_sd$exp_sd)), "Transcript"]
tmin <- A22_bim_sd[which(A22_bim_sd$exp_sd == min(A22_bim_sd$exp_sd)), "Transcript"]
tran <- A22_bim_sd[sample(nrow(A22_bim_sd), 1), "Transcript"]

# extract data for these transcripts
setkey(A22.TPM.bim.dt, Transcript)
dfmax <- A22.TPM.bim.dt[tmax]
dfmin <- A22.TPM.bim.dt[tmin]
dfran <- A22.TPM.bim.dt[tran]

# fit linear models to each
lmax <- lm(TPM ~ val + I(val^2), data = dfmax)
lmin <- lm(TPM ~ val + I(val^2), data = dfmin)
lran <- lm(TPM ~ val + I(val^2), data = dfran)

# coefficients
coef(lmax)
coef(lmin)
coef(lran)

# new data to predict
temp.int <- seq(from = 0, to = 38.5, by = 0.5)
temp.df <- data.frame(val = temp.int)

# predict
pout.max <- predict(lmax, newdata=temp.df)
pout.min <- predict(lmin, newdata=temp.df)
pout.ran <- predict(lran, newdata=temp.df)

max(c(pout.max, pout.min, pout.ran))

# plot
plot(temp.int, pout.min, col = "blue", ylim = c(0,0.5), type = "l", ylim = c(0, ))
lines(temp.int, pout.max, col = "red", lty = 2)
lines(temp.int, pout.ran, col = "green", lty = 3)

# Plot shows that the transcript with *minimum* standard deviation is almost a flat line, while the transcript with the *maximum* standard deviation of expression is highly responsive.