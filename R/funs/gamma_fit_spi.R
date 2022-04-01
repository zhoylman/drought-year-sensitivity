#function to fit a gamma distrobuion to a vector of data
#export options (export_opts) allows the user to return 
#SPI valules if export_opts = 'SPI', CDF values if export_opts = 'CDF'
#or the gamma distrobution paramters if export_opts = 'params'.
#the function also allows the user to return either the latest
#CDF or SPI values when return_latest = T. when return_latest = F
#the entire SPI or CDF vector is returned. Default is to return latest. 

gamma_fit_spi = function(x, export_opts = 'SPI', return_latest = T) {
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
      if(return_latest == T){
        if(export_opts == 'CDF'){
          return(fit.cdf[length(fit.cdf)]) 
        }
        if(export_opts == 'params'){
          return(fit.gam) 
        }
        if(export_opts == 'SPI'){
          return(spi[length(spi)]) 
        }
      }
      if(return_latest == F){
        if(export_opts == 'CDF'){
          return(fit.cdf) 
        }
        if(export_opts == 'params'){
          return(fit.gam) 
        }
        if(export_opts == 'SPI'){
          return(spi) 
        }
      }
      
    },
    #else return NA
    error=function(cond) {
      return(NA)
    })
}

gamma_fit_spi_mod = function(x, export_opts = 'SPI', return_latest = T) {
  #load the package needed for these computations
  library(lmomco)
  #first try gamma
  tryCatch(
    {
      x = as.numeric(x)
      # find zeros for replacement
      if(any(x == 0, na.rm = T)){
        index = which(x == 0)
        non_index =  which(x != 0)
        x_cont = x[-index]
        #Unbiased Sample Probability-Weighted Moments (following Beguer ́ıa et al 2014)
        pwm = pwm.ub(x_cont)
        #Probability-Weighted Moments to L-moments
        lmoments_x = pwm2lmom(pwm)
        #fit gamma
        fit.gam = pargam(lmoments_x)
        #compute probabilistic cdf 
        fit.cdf = cdfgam(x_cont, fit.gam)
        #reconstruct values based on mixed distrobution assumption (Wu et. al., 2007)
        q = length(index)/length(x)
        fit.cdf.mod = x
        #replace values based on index defined above
        fit.cdf.mod[index] = q
        fit.cdf.mod[non_index] = q+((1-q)*fit.cdf)
        #compute spi
        spi = qnorm(fit.cdf.mod, mean = 0, sd = 1)
        # #compute t
        # t = rep(NA, length(fit.cdf.mod))
        # t = ifelse(fit.cdf.mod > 0.5, 
        #            sqrt(log(1/((1 - fit.cdf.mod)^2))),
        #            sqrt(log(1/(fit.cdf.mod ^2))))
        # #compute spi
        # c0=2.515517
        # c1=0.802853
        # c2=0.010328
        # d1=1.432788
        # d2=0.189269
        # d3=0.001308
        # 
        # spi = rep(NA, length(fit.cdf.mod))
        # spi = ifelse(fit.cdf.mod > 0.5,
        #              ((t -(c0 + c1*t + c2*t^2)/(1 + d1*t + d2*t^2 + d3*t^3))),
        #              -((t -(c0 + c1*t + c2*t^2)/(1 + d1*t + d2*t^2 + d3*t^3))))
        # 
      } else {
        x = as.numeric(x)
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
      }
      if(return_latest == T){
        if(export_opts == 'SPI'){
          return(spi[length(spi)]) 
        }
      }
      if(return_latest == F){
        if(export_opts == 'SPI'){
          return(spi) 
        }
      }
      
    },
    #else return NA
    error=function(cond) {
      return(NA)
    })
}
