#' Assessing Risk Model Calibration with Missing Covariates
#'
#' Assess risk model calibration using a ratio of observed to expected number of cases
#' where some risk model covariates are missing due to two-phase sampling designs such as case-cohort or nested case-control designs. 
#' Default is to use inverse probability weights (also known as sampling weights, survey weights, or designs weights).
#' Can use adjusted weights where the adjusted-weighted sum and the total sum are equal for pseudo risk values (using \code{\link{calib}}).
#' Can also compare with multiple imputation approach (e.g., \code{\link{mice}}). 
#' Standard errors are provided for any available approaches.
#' @param risk risk values computed from a risk model to validate with a two-phase designed validation cohort. Will be used to estimate the expected number of cases. Can be a numeric vector of length \code{n} without missing components or length \code{N} with NA for the individuals with missing covariates.
#' @param ind.case binary outcome status (1 = case; 0 = non-case). Will be used to compute the observed number of cases. This should be a vector of length \code{N}. 
#' @param ind.ph2 binary phase 2 inclusion indicators (1 = included; 0 = not included). This should be a vector of length \code{N}. 
#' @param incl.prob phase 2 inclusion probabilities. This should be a vector of length \code{N}. 
#' @param survtime follow-up times (i.e., time to event where event can be either case or loss to follow-up due to administrative censoring or other causes). This should be a vector of length \code{N}. 
#' @param tau risk projection time used to compute \code{risk}. The observed number of cases will be computed only up to this projection time. If NULL, risk is projected up to the maximum of \code{survtime}.
#' @param pseudo.risk pseudo risk values computed from a risk model to validate with a validation cohort where missing covariates are substituted by predicted values e.g., using \code{\link{glm}}.  Will be used as auxiliary statistics for weight adjustment using \code{\link{calib}}. This should be a numeric vector of length \code{N}.
#' @param imputed.risks multiple-imputed risk values e.g., by \code{\link{mice}}. This should be an \code{N} x \code{K} matrix where \code{K} is the number of imputations.
#' @param m.ifncc number of sampled controls per case if missing is by nested case-control samplign designs. Default is NULL.
#' @keywords ModelCal
#' @return estimates and standard errors of the ratio of observed and expected number of cases using inverse probability weights, pseudo risk-adjusted weights (if \code{pseudo.risk} is not NULL), and multiple imputations (if \code{imputed.risks} is not NULL). 
#' @export 
#' @seealso \code{\link{data_valid}}, \code{\link{riskmod}}
#' @examples #### Assessing Calibration of Pure Cumulative Risk Model using Example Two-phased Validation Data ####
#' tau = 30 # tau can range from 0 to 50 for the example data
#' data_valid.hat = data_valid
#' data_valid.hat$eventime = pmin(tau, data_valid$eventime)
#' # Risk Values in Phase 2
#' # library(survival)
#' rhat = predict(riskmod, newdata = data_valid.hat, type = "expected")
#' 
#' # Pseudo-risk Values in Phase 1 using Predicted Exposure Z
#' predmod = glm(Z ~ X + Uz, family = gaussian, data = data_valid, weights = wgt, subset = (ind.ph2==1))
#' data_valid.hat$Z = predict(predmod, data_valid, type = "response") 
#' prisk = predict(riskmod, newdata = data_valid.hat, type = "expected")
#' 
#' # Imputed Risk Values in Phase 1 using MICE
#' # library(mice)
#' data_valid.tilde = data_valid
#' data_valid.tilde$logT = log(data_valid.tilde$eventime)
#' ini = mice(data_valid.tilde, max = 0, print = FALSE)
#' ini$predictorMatrix[c("Z"), ] <- 0
#' ini$predictorMatrix["Z", c("X", "Uz", "logT", "ind.fail")] <- 1
#' miceobj = mice(data_valid.tilde, meth = ini$meth, pred = ini$predictorMatrix, seed = 2, m = 5, maxit = 5, print = F)
#' predmice = with(miceobj, predict(riskmod, newdata = data.frame(eventime = pmin(tau, eventime), ind.fail, X, Z), type = "expected"))
#' irisk = simplify2array(predmice$analyses)
#' 
#' # Assessing Model Calibration 
#' res = model.calibration(risk = rhat, 
#' ind.case = data_valid$ind.fail, 
#' ind.ph2 = data_valid$ind.ph2, 
#' incl.prob = data_valid$incl.prob,
#' survtime = data_valid$eventime,
#' tau = tau,
#' pseudo.risk = prisk,
#' imputed.risks = irisk)
#' 
#' print(round(res, 3))
#' 
model.calibration = function(risk, ind.case, ind.ph2, incl.prob, survtime, tau = NULL, pseudo.risk = NULL, imputed.risks = NULL, m.ifncc = NULL){
  
  if(!require(sampling)){
    install.packages("sampling")
    library(sampling)
  }
  
  n = length(ind.case)
  V = ind.ph2
  
  if(is.null(tau)){
    tau = max(survtime)
  }
  if(length(risk)!=n){
    temp = rep(0, n)
    temp[V==1] = risk
    risk = temp
  }
  if(any(is.na(risk))){
    temp = rep(0, n)
    temp[V==1] = risk[V==1]
    risk = temp
  }
  
  Oi = ind.case*(survtime<=tau)
  O = sum(Oi)
  
  #### Inclusion Probability Weighting ####
  w0 = 1/incl.prob
  E0i = V*w0*risk
  E0 = sum(V*w0*risk)
  
  if(!is.null(m.ifncc)){
    E0i_ncc = E0i[V == 1 & Oi == 0]
    omega = wcov.ncc(survtime = survtime, samplestat = ind.ph2 + ind.case, m = m.ifncc, psample = incl.prob)
    var_E0 = as.numeric(t(E0i_ncc) %*% omega %*% E0i_ncc) + n/(n-1)*sum(V*w0*(risk-mean(E0i))^2)
  }
  if(is.null(m.ifncc)){
    var_E0 = sum((1-1/w0)*E0i^2) + n/(n-1)*sum(V*w0*(risk-mean(E0i))^2)
  }
  cov_OE0 = cov(Oi, V*w0*risk)*(n-1)
  
  result = inclusion_weighting = c(O, E0, var_E0, cov_OE0)
  
  #### Weight Adjustment when Pseudo-risk is Given ####
  if(!is.null(pseudo.risk)){
    aux = cbind(1, pseudo.risk)
    g = wadj = rep(0, n)
    g[V==1] = calib(Xs = aux[V==1,], d = w0[V==1], total = colSums(aux), method = "raking")
    wadj[V==1] = g[V==1]*w0[V==1]
    
    Eadj = sum(V*wadj*risk)
    if(!is.null(m.ifncc)){
      var_Eadj = var.Eadj(risk[V==1], pseudo.risk, ind.ph2 = V, event = Oi, 
                          incl.prob = incl.prob, gcal = g[V==1], omega = omega)
    }
    if(is.null(m.ifncc)){
      var_Eadj = var.Eadj(risk[V==1], pseudo.risk, ind.ph2 = V, event = Oi, 
                          incl.prob = incl.prob, gcal = g[V==1])
    }
    cov_OEadj = cov(Oi, V*wadj*risk)*(n-1)
    
    result = rbind(result, weight_adjustment=c(O, Eadj, var_Eadj, cov_OEadj))
  }
  
  #### Multiple Imputation when Imputed-risks is Given ####
  if(!is.null(imputed.risks)){
    K = ncol(imputed.risks)
    Emi = sum(rowMeans(imputed.risks))
    var_Emi = mean(apply(imputed.risks, 2, var)*n)+(1+1/K)*var(colSums(imputed.risks)) # Rubin's formula
    cov_OEmi = cov(Oi, rowMeans(imputed.risks))*(n-1)
    
    result = rbind(result, multiple_imputation=c(O, Emi, var_Emi, cov_OEmi))
  }
  
  if(length(result) == 4){
    est = result[1]/result[2]
    se = ratiosd(result[1], result[2], result[3], result[4])
    ci_lower = est - qnorm(.975)*se
    ci_upper = est + qnorm(.975)*se
    result_summary = c(O = result[1], E = result[2], est, se, ci_lower, ci_upper)
    result_summary = matrix(result_summary, ncol = 6, 
                            dimnames = list("inclusion_weighting", c("O", "E", "est", "se", "ci_lower", "ci_upper")))
  }
  if(length(result) != 4){
    est = apply(result, 1, function(x){x[1]/x[2]})
    se = apply(result, 1, function(x){ratiosd(x[1], x[2], x[3], x[4])})
    ci_lower = est - qnorm(.975)*se
    ci_upper = est + qnorm(.975)*se
    result_summary = cbind(O = result[,1], E = result[,2], est, se, ci_lower, ci_upper)
  }
  
  return(result_summary)
}


