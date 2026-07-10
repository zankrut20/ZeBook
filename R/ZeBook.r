#' @title Working with Dynamic Models for Agriculture and Environment
#' @description
#' \strong{ZeBook}
#' Working with Dynamic Models for Agriculture and Environment
#' (Working with Dynamic Crop Models)

#'
#' Linked to book \strong{Working with Dynamic Crop Models} (Elsevier), Third edition, 27 septembre 2018 by Wallach, Makowski, Jones and Brun. \url{https://www.modelia.org/moodle/course/view.php?id=61}
#'
#' A full description of the models is in the book in appendix of the book.
#'
#' Chapter numbers have changed between Second edition and Third Edition. Here the chapter numbers in the demo were changed to fit to Third edition. But all materials available in Second edition are still available in this version. 
#'
#' \strong{ACKNOWLEDGMENTS}
#' The project "Associate a level of error in predictions of models for agronomy" (CASDAR 2010-2013) and the French network "RMT modeling and agriculture", \url{https://www.modelia.org}) have contributed to the development of this R package. This project and network are lead by ACTA (French Technical Institute for Agriculture) and was funded by a grant from the Ministry of Agriculture and Fishing of France.
#'
#' \strong{Other contributions}
#' Juliette Adrian, Master2 internship (ACTA, \url{https://www.modelia.org/moodle/mod/resource/view.php?id=1027}), january-jully 2013.
#'
#' Sylvain Toulet, Master2 internship (INRAE, \url{https://www.modelia.org/moodle/mod/resource/view.php?id=965}), january-jully 2012.
#'                                     
#' @aliases ZeBook
#' @author  Francois Brun (ACTA)  \email{francois.brun@@acta.asso.fr}, David Makowski (INRAE), Daniel Wallach (INRAE), James W. Jones (U.of Florida),
#' @references Working with Dynamic Crop Models (Elsevier), Third edition. 2019
#' \url{https://www.modelia.org}
#' @keywords models agricultural agronomy crop environment methods tools evaluation uncertainty sensitivity parameter estimation bayesian assimilation
#' @importFrom graphics hist plot plot.default lines legend par abline barplot text mtext axis
#' @importFrom stats approx cor cov na.omit runif sd var
"_PACKAGE"
