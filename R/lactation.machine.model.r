################################ FUNCTIONS #####################################
# Contribution of Juliette Adrian, Master2 internship, january-jully 2013
#' @title The Lactation model with milking machine
#' @description \strong{Model description.}
#' This model is a model of lactating mammary glands of cattle described by Heather et al. (1983). This model was then inspired more complex models based on these principles.
#' This model simulates the dynamics of the production of cow's milk.
#' the system is represented by 6 state variables: change in hormone levels (H), the production and loss of milk secreting cells (CS), and removing the secretion of milk (M), the average quantity of milk contained in the animal (Mmean), the amount of milk removed (RM) and yield (Y).
#' The model has a time step dt = 0.001 for milking machines.
#' The model is defined by a few equations, with a total of twenty parameters for the described process.
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
#' @param rma : parameter of milk m (t) function
#' @param t1 : parameter of milk m (t) function
#' @param t2 : parameter of milk m (t) function
#' @param t3 : parameter of milk m (t) function
#' @param t4 : parameter of milk m (t) function
#' @param t5 : parameter of milk m (t) function
#' @param t6 : parameter of milk m (t) function
#' @param duration : duration of simulation
#' @param dt : time step
#' @param CSi : initial Number of secretory cells
#' @param Mi : initial Quantity of milk in animal (kg)
#' @return matrix with CS,M,Mmoy,RM
#' @export
lactation.machine.model <- function(cu,kdiv,kdl,kdh,km,ksl,kr,ks,ksm,mh,mm,p,mum,rma,t1,t2,t3,t4,t5,t6,duration,dt,CSi,Mi) {
  .lactation_engine(cu, kdiv, kdl, kdh, km, ksl, kr, ks, ksm, mh, mm, p, mum,
                    duration, dt, is_machine = TRUE,
                    rma = rma, t1 = t1, t2 = t2, t3 = t3, t4 = t4, t5 = t5, t6 = t6,
                    CSi = CSi, Mi = Mi)
}
# end of file

################################################################################
#' @title The Lactation model for use with lactation.machine.simule
#' @description see lactation.calf.model for model description.
#' @param param : a vector of parameters containning (cu,kdiv,kdl,kdh,km,ksl,kr,ks,ksm,mh,mm,p,mum,rma,t1,t2,t3,t4,t5,t6)(see lactation.model.machine)
#' @param duration : duration of simulation
#' @param dt : time step
#' @param CSi : initial Number of secretory cells
#' @param Mi : initial Quantity of milk in animal (kg)
#' @return data.frame with CS, M, Mmoy, RM, day, week
#' @export
lactation.machine.model2=function(param,duration,dt,CSi,Mi)
{
 # use lactation.model.machine function to run the model
return(lactation.machine.model(param["cu"],param["kdiv"],param["kdl"],param["kdh"],param["km"],param["ksl"],param["kr"],param["ks"],param["ksm"],param["mh"],param["mm"],param["p"],param["mum"],param["rma"],param["t1"],param["t2"],param["t3"],param["t4"],param["t5"],param["t6"],duration,dt,CSi,Mi))
}
# End of file
