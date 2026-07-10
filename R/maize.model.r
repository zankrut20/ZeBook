################################################################################
# "Working with dynamic models for agriculture"
# R script for practical work
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2010-08-09
# Model described in the book, Appendix. Models used as illustrative examples: description and R code
################################ FUNCTIONS #####################################
#' @title The basic Maize model.
#' @description \strong{Model description.}
#' This model is a dynamic model of crop growth for Maize cultivated in potential conditions.
#' The crop growth is represented by three state variables, leaf area per unit ground area (leaf area index, LAI), total biomass (B) and cumulative thermal time since plant emergence (TT). It is based on key concepts included in most crop models, at least for the "potential production" part. In fact, this model does not take into account any effects of soil water, nutrients, pests, or diseases,... 
#' @details The tree state variables are dynamic variables depending on days after emergence: TT(day), B(day), and LAI(day). The model has a time step dt of one day.\cr
#' The model is defined by a few equations, with a total of seven parameters for the described process.
#' \cr (1) \eqn{TT(day+1) = TT(day)+dTT(day)}{TT(day+1) = TT(day)+dTT(day)}
#' \cr (2) \eqn{B(day+1) = B(day)+dB(day)}{B(day+1) = B(day)+dB(day)}
#' \cr (3) \eqn{LAI(day+1) = LAI(day)+dLAI(day)}{LAI(day+1) = LAI(day)+dLAI(day)}
#' \cr (4) \eqn{dTT(day) = \max(\frac{TMIN(day)+TMAX(day)}{2}-Tbase;0)}{dTT(day) = max((TMIN(day)+TMAX(day))/2-Tbase ; 0)}
#' \cr (5) \eqn{dB(day) = RUE*(1-e^{-K*LAI(day)*I(day)}),\ if\ TT(day)\le TTM}{dB(day) = RUE*(1-e^{-K*LAI(day)*I(day)}), if TT(day)<= TTM} 
#' \cr \eqn{dB(day) = 0,\ if\  TT(day)>TTM}{dB(day) = 0, if TT(day)>TTM}
#' \cr (6) \eqn{dLAI(day) = alpha*dTT(day)*LAI(day)*\max(LAImax-LAI(day);0),\ if \ TT(day)\le TTL }{alpha*dTT(day)*LAI(day)*max(LAImax-LAI(day);0), if TT(day)<= TTL }
#' \cr \eqn{dLAI(day) = 0,\ if\  TT(day)>TTL }{dLAI(day) = 0 if TT(day)>TTL}
#' @param Tbase parameter the baseline temperature for growth (degreeCelsius)
#' @param TTM parameter temperature sum for crop maturity (degreeC.day)
#' @param TTL parameter temperature sum at the end of leaf area increase (degreeC.day)
#' @param K parameter extinction coefficient (relation between leaf area index and intercepted radiation) (-)
#' @param RUE parameter radiation use efficiency (?)
#' @param alpha parameter the relative rate of leaf area index increase for small values of leaf area index (?)
#' @param LAImax parameter maximum leaf area index (-)
#' @param weather weather data.frame for one single year
#' @param sdate sowing date
#' @param ldate last date
#' @return data.frame with daily TT, LAI,B
#' @seealso \code{\link{maize.model2}}, \code{\link{maize.define.param}}, \code{\link{maize.simule}}, \code{\link{maize.multisy}},
#' \code{\link{maize.simule240}},\code{\link{maize.simule_multisy240}}
#' @export
#' @examples 
#' weather = maize.weather(working.year=2010, working.site=30,weather_all=weather_EuropeEU)
#' maize.model(Tbase=7, RUE=1.85, K=0.7, alpha=0.00243, LAImax=7, TTM=1200, TTL=700,
#'   weather, sdate=100, ldate=250)
# ── Internal simulation engine ─────────────────────────────────────────────────
#
# All four exported maize model variants share the identical simulation loop;
# they differ only in:
#   (a) how RUE is computed (constant vs temperature-dependent via maize.RUEtemp)
#   (b) whether CumInt (cumulative intercepted radiation) is tracked
#   (c) whether BE (ear biomass) is tracked
#
# This engine centralises the loop once.  Exported functions are thin wrappers
# that set feature flags and pick the columns they need from the result.
#
# @param Tbase,RUE_or_max,K,alpha,LAImax,TTM,TTL  model parameters
# @param weather  data.frame with columns Tmin, Tmax, I
# @param sdate,ldate  integer day-of-year indices
# @param temp_rue    logical — use maize.RUEtemp() instead of constant RUE
# @param track_cumint logical — accumulate CumInt state variable
# @param track_ear   logical — accumulate BE (ear biomass) state variable
# @return named list with vectors TT, LAI, B (and optionally CumInt, BE),
#   each of length (ldate - sdate + 1) representing days sdate..ldate
# @keywords internal
.maize_engine <- function(Tbase, RUE_or_max, K, alpha, LAImax, TTM, TTL,
                          weather, sdate, ldate,
                          temp_rue    = FALSE,
                          track_cumint = FALSE,
                          track_ear   = FALSE)
{
  # Pre-allocate full-length state vectors (indexed 1..ldate)
  TT  <- rep(NA_real_, ldate)
  B   <- rep(NA_real_, ldate)
  LAI <- rep(NA_real_, ldate)

  # Optional state variables — only allocate when needed
  if (track_cumint) CumInt <- rep(NA_real_, ldate)
  if (track_ear)    BE     <- rep(NA_real_, ldate)

  # Initial conditions at sowing date
  TT[sdate]  <- 0.0
  B[sdate]   <- 1.0
  LAI[sdate] <- 0.01
  if (track_cumint) CumInt[sdate] <- 0.0
  if (track_ear)    BE[sdate]     <- 0.0

  # Extract weather columns once to avoid repeated $ dispatch inside the loop
  w_Tmin <- weather$Tmin
  w_Tmax <- weather$Tmax
  w_I    <- weather$I

  # Simulation loop
  for (day in sdate:(ldate - 1)) {
    # Thermal time increment
    dTT <- max((w_Tmin[day] + w_Tmax[day]) / 2 - Tbase, 0)

    # Radiation-use efficiency: constant or temperature-dependent
    if (temp_rue) {
      tday <- (w_Tmin[day] + w_Tmax[day]) / 2
      RUE  <- maize.RUEtemp(tday, RUE_or_max, 6.2, 16.5, 33, 44)
    } else {
      RUE  <- RUE_or_max
    }

    # Absorbed PAR term (shared by dB and dCumInt)
    fI <- 1 - exp(-K * LAI[day])

    # Biomass increment: zero after maturity
    dB <- if (TT[day] <= TTM) RUE * fI * w_I[day] else 0

    # LAI increment: zero after TTL
    dLAI <- if (TT[day] <= TTL)
      alpha * dTT * LAI[day] * max(LAImax - LAI[day], 0)
    else 0

    # Update core state variables
    TT[day  + 1] <- TT[day]  + dTT
    B[day   + 1] <- B[day]   + dB
    LAI[day + 1] <- LAI[day] + dLAI

    # Optional state updates
    if (track_cumint)
      CumInt[day + 1] <- CumInt[day] + w_I[day] * fI

    if (track_ear) {
      dBE <- if (TT[day] > TTL) dB else 0
      BE[day + 1] <- BE[day] + dBE
    }
  }
  # End simulation loop

  # Return only the simulation window sdate..ldate
  idx <- sdate:ldate
  out <- list(TT = TT[idx], LAI = LAI[idx], B = B[idx])
  if (track_cumint) out$CumInt <- CumInt[idx]
  if (track_ear)    out$BE     <- BE[idx]
  out
}

################################################################################
#' @title The basic Maize model.
#' @description \strong{Model description.}
#' This model is a dynamic model of crop growth for Maize cultivated in
#' potential conditions. The crop growth is represented by three state
#' variables, leaf area per unit ground area (leaf area index, LAI), total
#' biomass (B) and cumulative thermal time since plant emergence (TT).
#' @details The tree state variables are dynamic variables depending on days
#' after emergence: TT(day), B(day), and LAI(day). The model has a time
#' step dt of one day.
#' \cr (1) \eqn{TT(day+1) = TT(day)+dTT(day)}
#' \cr (2) \eqn{B(day+1) = B(day)+dB(day)}
#' \cr (3) \eqn{LAI(day+1) = LAI(day)+dLAI(day)}
#' \cr (4) \eqn{dTT(day) = \max(\frac{TMIN+TMAX}{2}-Tbase;0)}
#' \cr (5) \eqn{dB(day) = RUE*(1-e^{-K*LAI})*I,\ if\ TT\le TTM;\  0\ otherwise}
#' \cr (6) \eqn{dLAI(day) = alpha*dTT*LAI*\max(LAImax-LAI;0),\ if\ TT\le TTL;\ 0\ otherwise}
#' @param Tbase baseline temperature for growth (degreeC)
#' @param RUE radiation use efficiency (g.MJ-1)
#' @param K extinction coefficient (-)
#' @param alpha relative rate of LAI increase for small LAI ((degreeC.day)-1)
#' @param LAImax maximum leaf area index (m2 leaf/m2 soil)
#' @param TTM temperature sum for crop maturity (degreeC.day)
#' @param TTL temperature sum at end of leaf area increase (degreeC.day)
#' @param weather weather data.frame for one single year
#' @param sdate sowing date
#' @param ldate last date
#' @return data.frame with daily TT, LAI, B
#' @seealso \code{\link{maize.model2}}, \code{\link{maize.define.param}},
#'   \code{\link{maize.simule}}, \code{\link{maize.multisy}},
#'   \code{\link{maize.simule240}}, \code{\link{maize.simule_multisy240}}
#' @export
#' @examples
#' weather = maize.weather(working.year=2010, working.site=30,
#'   weather_all=weather_EuropeEU)
#' maize.model(Tbase=7, RUE=1.85, K=0.7, alpha=0.00243, LAImax=7,
#'   TTM=1200, TTL=700, weather, sdate=100, ldate=250)
maize.model <- function(Tbase, RUE, K, alpha, LAImax, TTM, TTL,
                        weather, sdate, ldate)
{
  e <- .maize_engine(Tbase, RUE, K, alpha, LAImax, TTM, TTL,
                     weather, sdate, ldate,
                     temp_rue = FALSE, track_cumint = FALSE, track_ear = FALSE)
  data.frame(day = sdate:ldate, TT = e$TT, LAI = e$LAI, B = e$B)
}

################################################################################
#' @title The basic Maize model for use with maize.simule
#' @description Wrapper pour maize.model
#' @param param : a vector of parameters
#' @param weather : weather data.frame for one single year
#' @param sdate : sowing date
#' @param ldate : last date
#' @return data.frame with daily TT, LAI,B
#' @export
#' @examples 
#' weather = maize.weather(working.year=2010, working.site=30,weather_all=weather_EuropeEU)
#' maize.model2(maize.define.param()["nominal",], weather, sdate=100, ldate=250)
maize.model2<-function(param, weather,sdate,ldate)
{
  # 7 Parameter values of the model, read from the param vector
  Tbase <- param["Tbase"]
  RUE <- param["RUE"]
  K <- param["K"]
  alpha <- param["alpha"]
  LAImax <- param["LAImax"]
  TTM <- param["TTM"]
  TTL <- param["TTL"]
  # use maize.model function to run the model
  return(maize.model(Tbase, RUE, K, alpha, LAImax, TTM, TTL, weather, sdate, ldate))
}
################################################################################
#' @title Define values of the parameters for the Maize model
#' @description Define parameters values
#' @return matrix with parameter values (nominal, binf, bsup)
#' @export
maize.define.param <- function()
{
# nominal, binf, bsup
# Tbase  : the baseline temperature for growth (degreeC)
Tbase <- c(7, 6, 8)
# RUE : radiation use efficiency (g.MJ-1)
RUE <- c(1.85,1.5,2.5)
# K : extinction coefficient (-)
K <- c(0.7,0.6,0.8)
#alpha : the relative rate of leaf area index increase for small values of leaf area index ((degreeC.day)-1)
alpha <- c(0.00243,0.002,0.003)
#LAImax : maximum leaf area index (m2 leaf/m2 soil)
LAImax <- c(7.0,6.0,8.0)
#TTM :  temperature sum for crop maturity (degreeC.day)
TTM <- c(1200,1100,1400)
#TTL : temperature sum at the end of leaf area increase (degreeC.day)
TTL <- c(700,600,850)
param <- data.frame(Tbase, RUE, K, alpha, LAImax, TTM, TTL)
return(.make_param_matrix(param))
}
################################################################################
#' @title Wrapper function to run Maize model for multiple sets of parameter values
#' @description wrapper for maize.model2
#' @param X : matrix of n row vectors of 7 parameters
#' @param weather : weather data.frame for one single year
#' @param sdate : sowing date
#' @param ldate : last date
#' @param all : if you want a matrix combining X and output (default = FALSE)
#' @return matrix with maximum biomass for each parameter vector
#' @export
maize.simule <- function(X, weather, sdate, ldate, all=FALSE){
  # output: maximum biomass only
  .apply_simule(X,
    fn          = function(v) max(maize.model2(v[1:7], weather, sdate, ldate)$B, na.rm=TRUE),
    output_name = "B",
    all         = all)
}
#' @title Wrapping function to run maize model on several site-years
#' @description Wrapping function to run maize model on several site-years
#' @param param : a vector of parameters
#' @param list_site_year : vector of site-year
#' @param sdate : sowing date
#' @param ldate : last date
#' @param weather_all : weather data.frame for corresponding site-years
#' @return a data.frame with simulation for all site-years, with the first column sy indicating the site-years
#' @export
maize.multisy <- function(param, list_site_year, sdate, ldate, weather_all=NA){
  n    <- length(list_site_year)
  # Pre-allocate a list of the correct length; assign results by index.
  # Avoids O(n^2) incremental rbind copying that grows the data.frame each step.
  sims <- vector("list", n)
  for (i in seq_len(n)) {
    sy       <- list_site_year[i]
    sy_parts <- .parse_siteyear(sy)
    weather  <- maize.weather(working.year = sy_parts["year"],
                              working.site = sy_parts["site"],
                              weather_all  = weather_all)
    result   <- maize.model2(param, weather, sdate, ldate)
    sims[[i]] <- cbind(sy, result)
  }
  do.call(rbind, sims)
}
################################################################################
#' @title Wrapper function to run Maize model multiple times for multiple sets of parameter values and give Biomass at day240
#' @description Wrapper function for multiple simulation with Maize model
#' @param X : matrix of n row vectors of 7 parameters
#' @param weather : weather data.frame for one single year
#' @param sdate : sowing date
#' @param ldate : last date
#' @param all : if you want a matrix combining X and output (default = FALSE)
#' @return a matrix of biomass at day=240 for all combinations of parameters of X
#' @export
#' @examples sy="18-2006"
#' weather = maize.weather(working.year=strsplit(sy,"-")[[1]][2],
#'   working.site=strsplit(sy,"-")[[1]][1],weather_all=weather_EuropeEU)
#' maize.simule240(maize.define.param(),weather, sdate=100, ldate=250, all=FALSE)
maize.simule240 <- function(X, weather, sdate, ldate, all=FALSE){
  # output: biomass at day 240
  .apply_simule(X,
    fn          = function(v) maize.model2(v[1:7], weather, sdate, ldate)[240-sdate+1, "B"],
    output_name = "B",
    all         = all)
}
################################################################################
#' @title Wrapper function to run Maize model for multiple sets of input variables (site-year) and give Biomass at day240.
#' @description Wrapper function to run Maize model for multiple sets of input variables (site-year) and give Biomass at day240.
#' @param param a vector of parameters
#' @param liste_sy vector of site-year
#' @param sdate sowing date
#' @param ldate last date
#' @param weather_all weather data table used
#' @return mean biomass at day=240
#' @export
#' @examples maize.multisy240(maize.define.param()["nominal",],c("18-2006","64-2004") 
#' , sdate=100, ldate=250, weather_all=weather_EuropeEU)
maize.multisy240 <- function(param, liste_sy, sdate, ldate, weather_all=NA){
  Y <- sapply(liste_sy, function(sy) {
    sy_parts <- .parse_siteyear(sy)
    maize.model2(param,
                 maize.weather(working.year=sy_parts["year"],
                               working.site=sy_parts["site"],
                               weather_all=weather_all),
                 sdate, ldate)[240-sdate+1, "B"]
  })
  return(mean(as.matrix(Y)))
}
################################################################################
#' @title Wrapper function to run Maize model for multiple sets of parameter values (virtual design) and multiple sets of input variables (site-year) and give Biomass at day240
#' @description Wrapper function to run Maize model for multiple sets of input variables (site-year) and give Biomass at day240.
#' @param X matrix of n row vectors of 7 parameters
#' @param liste_sy vector of site-year
#' @param sdate sowing date
#' @param ldate last date
#' @param weather_all Weather data base
#' @param all if you want a matrix combining X and output (default = FALSE)
#' @return a matrix of mean biomass at day=240 for all combinations of parameters of X
#' @export
#' @examples maize.simule_multisy240(maize.define.param(),c("18-2006","64-2004"),
#'   sdate=100, ldate=250,weather_all=weather_EuropeEU,all=FALSE)
maize.simule_multisy240 <- function(X, liste_sy, sdate, ldate, weather_all=NA, all=FALSE){
  .apply_simule(X,
    fn          = function(v) maize.multisy240(v[1:7], liste_sy, sdate, ldate, weather_all),
    output_name = "B",
    all         = all)
}
################################################################################
#' @title The Maize model with additional state variable CumInt
#' @description Variant of the maize model
#' @param Tbase parameter the baseline temperature for growth (degreeCelsius)
#' @param TTM parameter temperature sum for crop maturity (degreeC.day)
#' @param TTL parameter temperature sum at the end of leaf area increase (degreeC.day)
#' @param K parameter extinction coefficient (relation between leaf area index and intercepted radiation) (-)
#' @param RUE parameter radiation use efficiency (?)
#' @param alpha parameter the relative rate of leaf area index increase for small values of leaf area index (?)
#' @param LAImax parameter maximum leaf area index (-)
#' @param  weather weather data.frame for one single year
#' @param sdate sowing date
#' @param ldate last date
#' @return data.frame with daily TT, LAI,B
#' @export
################################################################################
#' @title The Maize model with additional state variable CumInt
#' @description Variant of the maize model that also tracks cumulative
#'   intercepted radiation (CumInt). Internally delegates to the shared
#'   `.maize_engine()` with `track_cumint = TRUE`.
#' @param Tbase baseline temperature for growth (degreeC)
#' @param RUE radiation use efficiency (g.MJ-1)
#' @param K extinction coefficient (-)
#' @param alpha relative rate of LAI increase for small LAI ((degreeC.day)-1)
#' @param LAImax maximum leaf area index (m2 leaf/m2 soil)
#' @param TTM temperature sum for crop maturity (degreeC.day)
#' @param TTL temperature sum at end of leaf area increase (degreeC.day)
#' @param weather weather data.frame for one single year
#' @param sdate sowing date
#' @param ldate last date
#' @return data.frame with daily TT, LAI, B, CumInt
#' @export
maize_cir.model <- function(Tbase, RUE, K, alpha, LAImax, TTM, TTL,
                            weather, sdate, ldate)
{
  e <- .maize_engine(Tbase, RUE, K, alpha, LAImax, TTM, TTL,
                     weather, sdate, ldate,
                     temp_rue = FALSE, track_cumint = TRUE, track_ear = FALSE)
  data.frame(day = sdate:ldate, TT = e$TT, LAI = e$LAI,
             B = e$B, CumInt = e$CumInt)
}


###############################################################################
#' @title Calculate effect of temperature on RUE for Maize
#' @description Function to compute effect of temperature on RUE 
#' @param T : temperature
#' @param RUE_max : maximum value for RUE
#' @param T0 : temperature parameter
#' @param T1 : temperature parameter
#' @param T2 : temperature parameter
#' @param T3 : temperature parameter
#' @return RUE value
#' @export
maize.RUEtemp <- function(T, RUE_max,T0,T1,T2,T3)
	{
	RUE = ((T>=T0)*(T<T1))* RUE_max*(T-T0)/(T1-T0) + ((T>=T1)*(T<T2))*RUE_max + ((T>=T2)*(T<T3))*RUE_max*(T3-T)/(T3-T2)
	}
###############################################################################
#' @title The Maize model with temperature dependent RUE and CumInt
#' @description Variant of the maize.model
#' @param Tbase parameter the baseline temperature for growth (degreeCelsius)
#' @param TTM parameter temperature sum for crop maturity (degreeC.day)
#' @param TTL parameter temperature sum at the end of leaf area increase (degreeC.day)
#' @param K parameter extinction coefficient (relation between leaf area index and intercepted radiation) (-)
#' @param RUE_max parameter maximum radiation use efficiency (?)
#' @param alpha parameter the relative rate of leaf area index increase for small values of leaf area index (?)
#' @param LAImax parameter maximum leaf area index (-)
#' @param weather weather data.frame for one single year
#' @param sdate sowing date
#' @param ldate last date
#' @return data.frame with daily TT, LAI,B
#' @export
################################################################################
#' @title The Maize model with temperature dependent RUE and CumInt
#' @description Variant of the maize model where RUE is a function of
#'   temperature (via \code{\link{maize.RUEtemp}}) and CumInt is tracked.
#'   Internally delegates to the shared `.maize_engine()` with
#'   `temp_rue = TRUE, track_cumint = TRUE`.
#' @param Tbase baseline temperature for growth (degreeC)
#' @param RUE_max maximum radiation use efficiency (g.MJ-1)
#' @param K extinction coefficient (-)
#' @param alpha relative rate of LAI increase for small LAI ((degreeC.day)-1)
#' @param LAImax maximum leaf area index (m2 leaf/m2 soil)
#' @param TTM temperature sum for crop maturity (degreeC.day)
#' @param TTL temperature sum at end of leaf area increase (degreeC.day)
#' @param weather weather data.frame for one single year
#' @param sdate sowing date
#' @param ldate last date
#' @return data.frame with daily TT, LAI, B, CumInt
#' @export
maize_cir_rue.model <- function(Tbase, RUE_max, K, alpha, LAImax, TTM, TTL,
                                weather, sdate, ldate)
{
  e <- .maize_engine(Tbase, RUE_max, K, alpha, LAImax, TTM, TTL,
                     weather, sdate, ldate,
                     temp_rue = TRUE, track_cumint = TRUE, track_ear = FALSE)
  data.frame(day = sdate:ldate, TT = e$TT, LAI = e$LAI,
             B = e$B, CumInt = e$CumInt)
}


###############################################################################
#' @title The Maize model with temperature dependent RUE, CumInt and ear growth
#' @description Variant of the maize.model
#' @param Tbase parameter the baseline temperature for growth (degreeCelsius)
#' @param TTM parameter temperature sum for crop maturity (degreeC.day)
#' @param TTL parameter temperature sum at the end of leaf area increase (degreeC.day)
#' @param K parameter extinction coefficient (relation between leaf area index and intercepted radiation) (-)
#' @param RUE_max parameter maximum radiation use efficiency (?)
#' @param alpha parameter the relative rate of leaf area index increase for small values of leaf area index (?)
#' @param LAImax parameter maximum leaf area index (-)
#' @param  weather weather data.frame for one single year
#' @param sdate sowing date
#' @param ldate last date
#' @return data.frame with daily TT, LAI,B
#' @export
################################################################################
#' @title The Maize model with temperature dependent RUE, CumInt and ear growth
#' @description Variant of the maize model that additionally tracks BE, the
#'   biomass allocated to the ear after TTL.  Internally delegates to the
#'   shared `.maize_engine()` with
#'   `temp_rue = TRUE, track_cumint = TRUE, track_ear = TRUE`.
#' @param Tbase baseline temperature for growth (degreeC)
#' @param RUE_max maximum radiation use efficiency (g.MJ-1)
#' @param K extinction coefficient (-)
#' @param alpha relative rate of LAI increase for small LAI ((degreeC.day)-1)
#' @param LAImax maximum leaf area index (m2 leaf/m2 soil)
#' @param TTM temperature sum for crop maturity (degreeC.day)
#' @param TTL temperature sum at end of leaf area increase (degreeC.day)
#' @param weather weather data.frame for one single year
#' @param sdate sowing date
#' @param ldate last date
#' @return data.frame with daily TT, LAI, B, CumInt, BE
#' @export
maize_cir_rue_ear.model <- function(Tbase, RUE_max, K, alpha, LAImax, TTM, TTL,
                                    weather, sdate, ldate)
{
  e <- .maize_engine(Tbase, RUE_max, K, alpha, LAImax, TTM, TTL,
                     weather, sdate, ldate,
                     temp_rue = TRUE, track_cumint = TRUE, track_ear = TRUE)
  data.frame(day = sdate:ldate, TT = e$TT, LAI = e$LAI,
             B = e$B, CumInt = e$CumInt, BE = e$BE)
}


###############################################################################
#' @title Read weather data for the Maize model
#' @description Function to read weather data and format them for maize.model
#' @param working.year year for the subset of weather data (default=NA : all the year)
#' @param working.site site for the subset of weather data (default=NA : all the site)
#' @param weather_all weather data base (default=weather_FranceWest)
#' @return data.frame with daily weather data for one or several site(s) and for one or several year(s)
#' @export
# Reading Weather data function
maize.weather <- function(working.year=NA, working.site=NA, weather_all=NA)
    {
    # WEYR => year  WEDAY => day  SRAD => I  TMAX => Tmax  TMIN => Tmin
    # select only useful columns and rename to model-expected names
    weather <- weather_all[, c("idsite","GPSlatitude","GPSlongitude",
                               "WEYR","WEDAY","SRAD","TMAX","TMIN")]
    names(weather)[names(weather) == "WEDAY"] <- "day"
    names(weather)[names(weather) == "WEYR"]  <- "year"
    names(weather)[names(weather) == "SRAD"]  <- "I"
    names(weather)[names(weather) == "TMAX"]  <- "Tmax"
    names(weather)[names(weather) == "TMIN"]  <- "Tmin"
    # delegate site/year subsetting to the shared internal helper
    weather <- .filter_site_year(weather, working.year, working.site)
    return(weather)
    }

################################################################################
# End of file
