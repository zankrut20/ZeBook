################################################################################
# "Working with dynamic models for agriculture"
# R script for pratical work
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2010-08-09
# Model described in the book, Appendix. Models used as illustrative examples: description and R code
################################ FUNCTIONS #####################################
#' @title Generate a random plan as a data frame. Columns are parameters. Values have uniform distribution
#' @description according to minimal and maximal values defined in a model.factors matrix 
#' @param model.factors : matrix defining minimal (binf) and maximal values (bsup) for a set of p parameters
#' @param N : size of sample
#' @return parameter matrix of dim = (N, p)
#' @importFrom stats runif
#' @export
param.runif = function(model.factors, N){
  p  <- ncol(model.factors)
  nm <- colnames(model.factors)
  # Pre-allocate the full output matrix once, fill column-by-column.
  # This avoids the O(N*p^2) memory overhead of cbind-growing a data.frame
  # inside a loop.
  X <- matrix(NA_real_, nrow = N, ncol = p, dimnames = list(NULL, nm))
  for (j in seq_len(p)) {
    X[, j] <- runif(N,
                    min = model.factors["binf",  nm[j]],
                    max = model.factors["bsup",  nm[j]])
  }
  as.data.frame(X)
}
# end of function
################################################################################
#' @title Generate a random plan as a data frame. Columns are parameters. Values have triangle distribution
#' @description according to nominal, minimal and maximal values defined in a model.factors matrix
#' @param model.factors : matrix defining nominal, minimal (binf), maximal values (bsup) for a set of p parameters
#' @param N : size of sample
#' @return parameter matrix of dim = (N, p)
#' @importFrom triangle rtriangle
#' @export
param.rtriangle = function(model.factors, N)
{
  p  <- ncol(model.factors)
  nm <- colnames(model.factors)
  # Pre-allocate the full output matrix once, fill column-by-column.
  X <- matrix(NA_real_, nrow = N, ncol = p, dimnames = list(NULL, nm))
  for (j in seq_len(p)) {
    # rtriangle: a = lower limit, b = upper limit, c = mode
    X[, j] <- rtriangle(N,
                         a = model.factors["binf",    nm[j]],
                         b = model.factors["bsup",    nm[j]],
                         c = model.factors["nominal", nm[j]])
  }
  as.data.frame(X)
}
# end of function
################################################################################
#' @title Build the q.arg argument for the  FAST function (sensitivity analysis)
#' @description according to minimal and maximal values defined in a model.factors matrix 
#' @param model.factors : matrix defining minimal (binf) and maximal values (bsup) for a set of p parameters
#' @return a list of list
#' @export 
q.arg.fast.runif = function(model.factors){
  nm <- colnames(model.factors)
  p  <- length(nm)
  # Pre-allocate list of the correct length; assign by index.
  # Avoids O(p^2) repeated list copying from c(list, list(...)).
  q.arg.fast <- vector("list", p)
  for (j in seq_len(p)) {
    q.arg.fast[[j]] <- list(min = model.factors["binf", nm[j]],
                            max = model.factors["bsup", nm[j]])
  }
  q.arg.fast
}
# end of function
