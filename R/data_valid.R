#' @title Example Validation Data for 'rmodcal' Package
#' @description Example cohort with simulated case-cohort study where N = 10000 and 250 cases. This dataset is used to run the example for function  \code{\link{model.calibration}}.
#' @format A data frame with N = 10000 rows and 8 variables:
#' \describe{
#'   \item{X}{continuous covariate available in phase 1}
#'   \item{Z}{continuous covariate only available in phase 2; NA if missing}
#'   \item{Uz}{an ancillary predictor for Z2 available in phase 1} 
#'   \item{eventime}{time-to-event (event = either case or censored)}
#'   \item{ind.fail}{binary outcome status (1 = case; 0 = non-case)}
#'   \item{ind.ph2}{phase 2 inclusion status (1 = included as a case or selected control; 0 = otherwise)}
#'   \item{incl.prob}{inclusion probability }
#'   \item{wgt}{sampling weights (\code{=1/incl.prob})}
#'}
"data_valid"