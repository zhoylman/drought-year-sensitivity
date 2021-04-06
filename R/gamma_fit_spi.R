gamma_fit_spi = function(x, export_opts = 'SPI') {
  #load the package needed for these computations
  library(lmomco)
  #first try gamma
  tryCatch(
    {
      x = as.numeric(x)
      #if precip is 0, replace it with 0.01mm Really Dry
      if(any(x == 0, na.rm = T)){
        index = which(x == 0)
        x[index] = 0.01
      }
      #Unbiased Sample Probability-Weighted Moments (following Beguer ́ıa et al 2014)
      pwm = pwm.ub(x)
      #Probability-Weighted Moments to L-moments
      lmoments_x = pwm2lmom(pwm)
      #fit gamma
      fit.gam = pargam(lmoments_x)
      #compute probabilistic cdf 
      fit.cdf = cdfgam(x, fit.gam)
      #compute spi
      spi = qnorm(fit.cdf, mean = 0, sd = 1)
      if(export_opts == 'CDF'){
        return(fit.cdf[length(fit.cdf)]) 
      }
      if(export_opts == 'params'){
        return(fit.gam) 
      }
      if(export_opts == 'SPI'){
        return(spi[length(spi)]) 
      }
    },
    #else return NA
    error=function(cond) {
      return(NA)
    })
}
