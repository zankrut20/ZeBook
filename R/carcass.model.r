
# ── Internal simulation engine ─────────────────────────────────────────────────
#
# Consolidates the timestep logic for both carcass models.
#
# @keywords internal
.carcass_engine <- function(protcmax, protncmax, alphac, alphanc, gammac, gammanc,
                            lip0, lipc1, lipnc1, beta, delta, b0c, b1c, b0nc, b1nc, c0, c1,
                            duration, init_ProtC, init_LipC, init_ProtNC, init_LipNC,
                            cem = NULL, k = NULL, amW = NULL, energie = NULL) {
  # Pre-allocate state vectors
  ProtC  <- rep(NA_real_, duration)
  LipC   <- rep(NA_real_, duration)
  ProtNC <- rep(NA_real_, duration)
  LipNC  <- rep(NA_real_, duration)
  PV     <- rep(NA_real_, duration)

  # Initial conditions
  ProtC[1]  <- init_ProtC
  LipC[1]   <- init_LipC
  ProtNC[1] <- init_ProtNC
  LipNC[1]  <- init_LipNC

  has_energie <- !is.null(energie)

  # Simulation loop
  for (t in seq_len(duration)) {
    # Body-weight sub-model
    MDC     <- b0c  * ProtC[t]^b1c
    PoidsC  <- LipC[t]  + MDC
    MDNC    <- b0nc * ProtNC[t]^b1nc
    PoidsNC <- LipNC[t] + MDNC

    PVV   <- PoidsC + PoidsNC
    PV[t] <- c0 * PVV^c1

    # Energy calculations
    if (has_energie) {
      EMI <- energie$y[t]
      CPM <- amW * PV[t]^0.75
    } else {
      EMI <- cem * (0.0157 * (PV[t]^0.9) + 3.3161)
      CPM <- k * PV[t]^0.75
    }

    # Cache energy-availability ratio
    emi_ratio <- EMI / (CPM + EMI)

    # Cache log-growth terms for proteins
    logC  <- log(protcmax  / ProtC[t])
    logNC <- log(protncmax / ProtNC[t])

    # Carcass lipids
    LipCmax  <- (lip0 + lipc1  * (ProtC[t]  / protcmax))  * PoidsC
    logLipC  <- log(LipCmax  / LipC[t])

    # Non-carcass lipids
    LipNCmax <- (lip0 + lipnc1 * (ProtNC[t] / protncmax)) * PoidsNC
    logLipNC <- log(LipNCmax / LipNC[t])

    # Rates of change
    dPC  <- alphac  * ProtC[t]  * logC  * emi_ratio  - gammac  * ProtC[t]  * logC
    dLC  <- beta    * LipC[t]   * logLipC  * emi_ratio  - delta   * LipC[t]   * logLipC
    dPNC <- alphanc * ProtNC[t] * logNC * emi_ratio  - gammanc * ProtNC[t] * logNC
    dLNC <- beta    * LipNC[t]  * logLipNC * emi_ratio  - delta   * LipNC[t]  * logLipNC

    # Update state variables
    if (t < duration) {
      ProtC[t+1]  <- ProtC[t]  + dPC
      LipC[t+1]   <- LipC[t]   + dLC
      ProtNC[t+1] <- ProtNC[t] + dPNC
      LipNC[t+1]  <- LipNC[t]  + dLNC
    }
  }

  data.frame(
    time   = seq_len(duration),
    ProtC  = ProtC,
    LipC   = LipC,
    ProtNC = ProtNC,
    LipNC  = LipNC,
    PV     = PV
  )
}

################################ FUNCTIONS #####################################
#' @title The Carcass (growth of beef cattle) model
#' @description \strong{Model description.}
#'
#'
#'
#' The model is defined by 20 equations, with a total of 19 parameters for the described process.
#' @param protcmax amounts of protein in the carcass of the adult animal (kg)
#' @param protncmax amounts of protein in the 5th district of the adult animal (kg)
#' @param alphac maximum protein synthesis rate in the frame (excluding basal metabolism) (j-1)
#' @param alphanc maximum rate of protein synthesis in the 5th district (except basal metabolism) (j-1)
#' @param gammac Maximum rate of protein degradation in the frame (excluding basal metabolism) (j-1)
#' @param gammanc maximum rate of protein degradation in the 5th district (except basal metabolism) (j-1)
#' @param lip0 maximum lipid concentration to the theoretical physiological age (percent)
#' @param lipc1 increase coefficient of the maximum lipid concentration with the physiological age of the carcase (percent)
#' @param lipnc1 increase coefficient of the highest lipid concentration with physiological age area in the 5th (percent)
#' @param beta lipid synthesis rate (j-1)
#' @param delta lipid degradation rate (d-1)
#' @param k Parameter coefficient between the half-saturation of the Michaelis-Menten equation of the metabolic weight (MJ.kg^0.75)
#' @param b0c coefficient of the allometric equation linking mass and lipid-protein carcass
#' @param b1c exponent allometric equation linking mass and defatted protein carcass
#' @param b0nc coefficient of the allometric equation linking mass and lipid-protein 5th district
#' @param b1nc exponent allometric equation linking mass and lipid-protein 5th district
#' @param c0 coefficient of the allometric equation between live weight and live weight empty
#' @param c1 exponent allometric equation linking body weight and live weight empty
#' @param cem parameter for body weight
#' @param duration duration of simulation
#' @return data.frame with ProtC,LipC,ProtNC,LipNC,PV
#' @export
carcass.model <- function(protcmax, protncmax, alphac, alphanc, gammac, gammanc,
                          lip0, lipc1, lipnc1, beta, delta, k,
                          b0c, b1c, b0nc, b1nc, c0, c1, cem, duration)
{
  .carcass_engine(protcmax, protncmax, alphac, alphanc, gammac, gammanc,
                  lip0, lipc1, lipnc1, beta, delta, b0c, b1c, b0nc, b1nc, c0, c1,
                  duration,
                  init_ProtC = 30, init_LipC = 15, init_ProtNC = 15, init_LipNC = 8,
                  cem = cem, k = k)
}
# End of file
