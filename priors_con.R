#Function for more conservative priors

# Priors for Bernoulli models (logit link: coefficients are in log-odds)
prior_bern <- c(
  prior(normal(0, 2), class = "b"),                  # Allows large effects on odds (OR ~ 0.02 to 55)
  prior(normal(-1.5, 1.5), class = "Intercept")      # Baseline infection prob ~18% (range ~3% to 82%)
)

# Priors for Bernoulli model with continuous duration predictor
priors_bern_dur <- c(
  prior(normal(0, 2), class = "b"),                  # As above: broad prior on other categorical predictors
  prior(normal(0, 0.05), class = "b", coef = "DurMin"), # Tight prior: ~10% change in odds per minute
  prior(normal(-1.5, 1.5), class = "Intercept") )     # # Baseline infection prob in reference group ~18% (but allows wide range from near 0% to ~80%)

# Priors for Gamma models (log link: coefficients are multiplicative effects on EPG)
prior_gamma <- c(
  prior(normal(0, 3), class = "b"),                  # Allows 0.05× to 20× effects on mean EPG
  prior(normal(2.5, 1), class = "Intercept"),        # Baseline mean EPG ~12 (range ~3.5 to 33)
  prior(exponential(1), class = "shape")             # Shape prior favours overdispersion, weakly regularising
)

# Priors for Gamma model with continuous duration predictor
prior_dur_model <- c(
  prior(normal(0, 0.05), class = "b", coef = "DurMin"), # Tight prior: ~5% change in mean EPG per minute
  prior(normal(0, 2), class = "b"),                     # Broader prior for other categorical effects (0.14× to 7.4×)
  prior(normal(2.5, 1), class = "Intercept"),           # As above: baseline mean EPG ~12 (on natural scale)
  prior(exponential(1), class = "shape")                # As above: flexible dispersion control
)

# Model family prior selector
get_custom_prior_con <- function(family) {
  switch(
    family,
    "bernoulli" = prior_bern,
    "bern_dur"  = priors_bern_dur,
    "gamma"     = prior_gamma,
    "gam_dur"   = prior_dur_model,
    stop(paste("Unsupported family:", family))
  )
}

