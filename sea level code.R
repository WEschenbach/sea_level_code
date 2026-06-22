
# obtaining autocorrelation adjusted quadratic fits
fit_car1 <- gls(
	y ~ time_frac+I(time_frac^2/2),
	data        = dat,
	correlation = corCAR1(form = ~ time_frac)
)

# obtaining autocorrelation adjusted linear fits
fit_car12 <- gls(
	y ~ time_frac,
	data        = dat,
	correlation = corCAR1(form = ~ time_frac)
)

# comparing AIC/BIC for quadratic and piecewise linear fits to full dataset

fit_quadratic_vs_piecewise <- function(t, y) {
	# "t" is time, "y" is avereage sea level
	# Ensure required package
	if (!requireNamespace("strucchange", quietly = TRUE)) {
		stop("Package 'strucchange' is required but not installed.")
	}
	
	# Put in data frame
	dat <- data.frame(t = as.numeric(t), y = as.numeric(y))
	
	# ----- 1. Quadratic model -----
	# y = a0 + a1 * t + (a2/2) * t^2  )
	# This is equivalent (for AIC/BIC) to y ~ t + I(t^2)

	fit_car1 <- gls(
		y ~ t + I(t^2),
		data        = dat,
		correlation = corCAR1(form = ~ t)
	)
	thefit=summary(fit_car1)
	
	quad_aic <- thefit$AIC
	quad_bic <- thefit$AIC
	
	# ----- 2. Piecewise linear model with one breakpoint -----
	# Use strucchange to estimate a single structural break in the linear trend
	# Model: y ~ 1 + t, with break in the regression parameters
	library(strucchange)
	
	# Compute breakpoints (allowing at most 1 break)
	bp <- strucchange::breakpoints(y ~ t, data = dat, breaks = 1, h = 0.15)
	
	
	
	# Select the number of breaks
	# Here we force a single break, but we still let breakpoints() estimate its location.
	bp_index <- bp$breakpoints[1]  # location (index) of the single break
	
	if (is.na(bp_index)) {
		warning("No breakpoint detected; fitting only a single linear model.")
		piece_fit <- lm(y ~ t, data = dat)
	} else {
		# Construct a piecewise linear term:
		# y = beta0 + beta1 * t + beta2 * (t - t_break)
		t_break <- dat$t[bp_index]
		dat$t_break <- pmax(0, dat$t - t_break)
		
		# Fit piecewise linear model
		piece_fit <- lm(y ~ t + t_break, data = dat)
		piece_fit <- gls(
			y ~ t + t_break,
			data        = dat,
			correlation = corCAR1(form = ~ t)
		)
		thefit2=summary(piece_fit)
		# thefit$tTable
		thefit2$AIC
		thefit2$BIC
	}
	
	quad_k   <- length(coef(quad_fit))
	piece_k  <- length(coef(piece_fit))
	
	
	piece_aic <- thefit2$AIC
	piece_bic <- thefit2$BIC
	
	# ----- Return results -----
	list(
		quadratic = list(
			fit = quad_fit,
			AIC = quad_aic,
			BIC = quad_bic
		),
		piecewise = list(
			fit = piece_fit,
			break_index = if (exists("bp_index")) bp_index else NA_integer_,
			break_time  = if (exists("t_break")) t_break else NA_real_,
			AIC = piece_aic,
			BIC = piece_bic
		)
	)
}

