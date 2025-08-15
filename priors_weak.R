#Function for weak priors


log_mean_epg <- log(mean(filt$mean_sm[filt$mean_sm > 0], na.rm = TRUE))
log_sd_epg   <- log(sd(filt$mean_sm[filt$mean_sm > 0], na.rm = TRUE))



# Function to return priors based on model family
get_custom_prior <- function(family) {

  switch(
    family,
    
    # Bernoulli model: Infection status
    "bernoulli" = c(
      prior(normal(0, 10), class = "b"),                   # Weakly informative, broad prior for effects
      prior(normal(-1.5, 2.5), class = "Intercept") # Prior belief: baseline infection ~18%
    ),
    
    # Bernoulli model with continuous duration predictor
    "bern_dur" = c(
      prior(normal(0, 10), class = "b"),                    # Broad prior for all effects
      prior(normal(-1.5, 2.5), class = "Intercept")
    ),
    
    # Gamma model: Egg burden
    "gamma" = c(
      prior(normal(0, 10), class = "b"),                    # Broad prior for categorical effects
      prior(normal(2.4, 3), class = "Intercept"),  # Informed by data scale
      prior(exponential(1), class = "shape")                # Favour overdispersion
    ),
    
    # Gamma model with continuous duration predictor
    "gam_dur" = c(
      prior(normal(0, 10), class = "b"),                    # Broad prior for all effects
      prior(normal(2.4, 3), class = "Intercept"),
      prior(exponential(1), class = "shape")
    ),
    
    stop(paste("Unsupported family:", family))
  )
}

