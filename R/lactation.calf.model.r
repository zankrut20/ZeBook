################################ FUNCTIONS #####################################
# Contribution of Juliette Adrian, Master2 internship, january-jully 2013
## ── Internal simulation engine ─────────────────────────────────────────────────
#
# Consolidates the timestep logic for both lactation models.
#
# @keywords internal
.lactation_engine <- function(cu, kdiv, kdl, kdh, km, ksl, kr, ks, ksm, mh, mm, p, mum,
                              duration, dt, is_machine = FALSE,
                              rc = NULL, rma = NULL, t1 = NULL, t2 = NULL, t3 = NULL,
                              t4 = NULL, t5 = NULL, t6 = NULL, CSi = 520, Mi = 0.0) {
  n_steps <- as.integer(duration / dt)
  n_alloc <- n_steps + 2

  H    <- rep(NA_real_, n_alloc)
  CS   <- rep(NA_real_, n_alloc)
  M    <- rep(NA_real_, n_alloc)
  Mmoy <- rep(NA_real_, n_alloc)
  RM   <- rep(NA_real_, n_alloc)

  H[1]    <- 1.0
  CS[1]   <- CSi
  M[1]    <- Mi
  Mmoy[1] <- 0.0

  if (is_machine) {
    RM[1] <- 0.0
  }

  i <- 1
  for (t in seq(0, duration, by = dt)) {
    if (is_machine) {
      t_mod <- t %% 1
      if ((t_mod > t1 && t_mod < t2) || (t_mod > t3 && t_mod < t4) || (t_mod > t5 && t_mod < t6)) {
        mach <- rma
      } else {
        mach <- 0
      }
    } else {
      mach <- rc
    }

    dH    <- - kdh * H[i] * dt
    dCS   <- (mum * (H[i] / (kdiv + H[i])) * cu - (ks + ksm * ((Mmoy[i] / mh)^p / (1 + (Mmoy[i] / mh)^p))) * CS[i]) * dt
    dM    <- (km * CS[i] * ((mm - M[i]) / (mm - M[i] + ksl)) - (M[i] / (kdl + M[i])) * mach) * dt
    dMmoy <- kr * (M[i] - Mmoy[i]) * dt

    H[i + 1]    <- H[i] + dH
    CS[i + 1]   <- CS[i] + dCS
    M[i + 1]    <- M[i] + dM
    Mmoy[i + 1] <- Mmoy[i] + dMmoy

    if (is_machine) {
      RM[i + 1] <- (M[i + 1] / (kdl + M[i + 1])) * mach
    } else {
      RM[i] <- (M[i] / (kdl + M[i])) * mach
    }

    i <- i + 1
  }

  day  <- seq(dt, duration, by = dt)
  week <- day %/% 7

  results1 <- data.frame(
    M    = M[1:n_steps],
    Mmoy = Mmoy[1:n_steps],
    CS   = CS[1:n_steps],
    RM   = RM[1:n_steps],
    day  = day,
    week = week
  )

  result <- by(results1[, c("week", "M", "Mmoy", "CS", "RM")], results1$week, function(x) apply(x, 2, mean))
  matrix(unlist(result), ncol = 5, byrow = TRUE, dimnames = list(NULL, c("week", "M", "Mmoy", "CS", "RM")))
}

#' @title The Lactation model
#' @description \strong{Model description.}
#' This model is a model of lactating mammary glands of cattle described by Heather et al. (1983). This model was then inspired more complex models based on these principles.
#' This model simulates the dynamics of the production of cow's milk.
#' the system is represented by 6 state variables: change in hormone levels (H), the production and loss of milk secreting cells (CS), and removing the secretion of milk (M), the average quantity of milk contained in the animal (Mmean), the amount of milk removed (RM) and yield (Y).
#' The model has a time step dt = 0.1 for regular consumption of milk by a calf.
#' The model is defined by a few equations, with a total of fourteen parameters for the described process.
#' @param cu : number of undifferentiated cells
#' @param kdiv : cell division rate, Michaelis-Menten constant
#' @param kdl : constant degradation of milk
#' @param kdh : rate of decomposition of the hormone
#' @param km :  constant secretion of milk
#' @param ksl : milk secretion rate, Michaelis-Menten constant
#' @param kr : average milk constant
#' @param ks : rate of degradation of the basal cells
#' @param ksm : constant rate of degradation of milk secreting cells
#' @param mh : parameter
#' @param mm : storage Capacity milk the animal
#' @param p : parameter
#' @param mum : setting the maximum rate of cell division
#' @param rc : parameter of milk m (t) function
#' @param duration : duration of simulation
#' @param dt : time step
#' @return data.frame with CS, M, Mmoy, RM, day, week
#' @examples lactation.calf.model2(lactation.define.param()["nominal",],300,0.1)
#' @export
lactation.calf.model <- function(cu, kdiv, kdl, kdh, km, ksl, kr, ks, ksm, mh, mm, p, mum, rc, duration, dt) {
  .lactation_engine(cu, kdiv, kdl, kdh, km, ksl, kr, ks, ksm, mh, mm, p, mum,
                    duration, dt, is_machine = FALSE, rc = rc, CSi = 520, Mi = 0.0)
}
################################################################################
#' @title The Lactation model for use with lactation.calf.simule
#' @description see lactation.calf.model for model description.
#' @param param : a vector of parameters
#' @param duration : duration of simulation
#' @param dt : time step
#' @return data.frame with CS, M, Mmoy, RM, day, week
#' @examples sim=lactation.calf.model2(lactation.define.param()["nominal",],6+2*7, 0.1)
#' @export
lactation.calf.model2 <- function(param,duration,dt){
  if ("rc" %in% names(param)) {
    return(lactation.calf.model(param["cu"],param["kdiv"],param["kdl"],param["kdh"],param["km"],param["ksl"],param["kr"],param["ks"],param["ksm"],param["mh"],param["mm"],param["p"],param["mum"],param["rc"],duration,dt))
  } else {
    .lactation_engine(param["cu"],param["kdiv"],param["kdl"],param["kdh"],param["km"],param["ksl"],param["kr"],param["ks"],param["ksm"],param["mh"],param["mm"],param["p"],param["mum"],duration,dt,is_machine=TRUE,rma=param["rma"],t1=param["t1"],t2=param["t2"],t3=param["t3"],t4=param["t4"],t5=param["t5"],t6=param["t6"])
  }
}
################################################################################
#' @title Wrapper function to run the Lactation model for multiple sets of parameter values
#' @description Wrapper function to run the Lactation model for multiple sets of parameter values
#' @param X : parameter matrix
#' @param duration : duration of simulation
#' @param dt : time step
#' @return data.frame with : number of paramter vector (line number from X), week, CS, M, Mmoy, RM, day, week
#' @export
lactation.calf.simule = function(X, duration, dt){
# output : all
#sim <- apply(X,1,function(v) lactation.calf.model2(v,duration, dt))
#sim=do.call(rbind, sim)
sim <- lapply(1:dim(X)[1], function(id) cbind(id,lactation.calf.model2(X[id,],duration, dt)))
sim=do.call(rbind, sim)
return(sim)
}
################################################################################
#' @title Define values of the parameters for the Lactation model
#' @description values from Heather et al. (1983) for different scenarios
#' @param type : for which model version ? "calf" or "machine"
#' @return matrix with parameter values (nominal, binf, bsup)
#' @examples lactation.define.param()
#' @export
lactation.define.param <- function(type="calf")
{
# nominal, binf, bsup
#cu : number of undifferentiated cells (Unit ?)
cu=c(1000, NA, NA)
#kdiv : cell division rate, Michaelis-Menten constant (Unit ?)
kdiv=c(0.2, NA, NA)
#kdl : constant degradation of milk (Unit ?)
kdl=c(4.43, NA, NA)
#kdh : rate of decomposition of the hormone (Unit ?)
kdh=c(0.01, NA, NA)
#km : constant secretion of milk ksl : milk secretion rate, Michaelis-Menten constant (Unit ?)
km=c(0.005, NA, NA)
#
ksl=c(3.0, NA, NA)
#kr : average milk constant (Unit ?)
kr=c(0.048, NA, NA)
#ks : rate of degradation of the basal cells (Unit ?)
ks=c(0.1, NA, NA)
#ksm : constant rate of degradation of milk secreting cells (Unit ?)
ksm=c(0.2, NA, NA)
#mh : parameter mm : storage Capacity milk the animal (Unit ?)
mh=c(27, NA, NA)
#
mm=c(30, NA, NA)
#p : parameter mum : setting the maximum rate of cell division (Unit ?)
p=c(10, NA, NA)
#
mum=c(1, NA, NA)
  if (type == "calf") {
    rc=c(40, NA, NA)
    param <- data.frame(cu,kdiv,kdl,kdh,km,ksl,kr,ks,ksm,mh,mm,p,mum,rc)
  } else if (type == "machine") {
    rma=c(80, NA, NA)
    t1=c(0.25, NA, NA)
    t2=c(0.30, NA, NA)
    t3=c(0.75, NA, NA)
    t4=c(0.80, NA, NA)
    t5=c(0, NA, NA)
    t6=c(0, NA, NA)
    param <- data.frame(cu,kdiv,kdl,kdh,km,ksl,kr,ks,ksm,mh,mm,p,mum,rma,t1,t2,t3,t4,t5,t6)
  }
  return(.make_param_matrix(param))

}
# end of file
