################################################################################
# "Working with dynamic models for agriculture"
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2013-03-25
################################ FUNCTIONS #####################################
#' @title Calcule multiple goodness-of-fit criteria
#' @description This function is depreciated and will be remove from the package in future versions. Please use goodness.of.fit
#' @param Ypred prediction values from the model
#' @param Yobs observed values
#' @param draw.plot draw evaluation plot
#' @return data.frame with the different evaluation criteria
#' @importFrom stats var cov
#' @importFrom graphics abline barplot par
#' @export
#' @examples
#' # observed and simulated values
#' obs<-c(78,110,92,75,110,108,113,155,150)
#' sim<-c(126,126,126,105,105,105,147,147,147)
#' evaluation.criteria(sim,obs,draw.plot=TRUE)
evaluation.criteria=function(Ypred,Yobs,draw.plot=FALSE){
  .Deprecated(new = "goodness.of.fit")

  # Call new implementation
  res_gof <- goodness.of.fit(Yobs = Yobs, Ypred = Ypred, draw.plot = draw.plot)

  # we keep only Ypred and Yobs where both are not NA
  select=!is.na(Ypred)&!is.na(Yobs)
  Ypred_clean=Ypred[select]
  Yobs_clean=Yobs[select]
  Nobs<-length(Yobs_clean)

  var.Yobs<- var(Yobs_clean)*(Nobs-1)/Nobs
  std.Yobs<-var.Yobs^0.5

  var.Ypred<- var(Ypred_clean)*(Nobs-1)/Nobs
  std.Ypred<-var.Ypred^0.5

  SSE <- res_gof$MSE * res_gof$N

  covBYobs<-cov(Ypred_clean, Yobs_clean)*(Nobs-1)/Nobs
  r<-covBYobs/(var.Yobs^0.5*var.Ypred^0.5)

  return(data.frame(
    Nobs=res_gof$N,
    mean.Yobs=res_gof$mean.Yobs,
    mean.Ypred=res_gof$mean.Ypred,
    std.Yobs=std.Yobs,
    std.Ypred=std.Ypred,
    SSE=SSE,
    MSE=res_gof$MSE,
    RMSE=res_gof$RMSE,
    r=r,
    bias.Squared=res_gof$bias.squared,
    SDSD=res_gof$SDSD,
    LCS=res_gof$LCS,
    EF=res_gof$EF
  ))
}

################################################################################
# End of file
