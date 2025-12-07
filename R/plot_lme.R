#' Plot Model-Predicted Trajectories from a Specific LME Model
#'
#' This function generates a visualization of the predicted median multiplicative
#' change over time, along with 95% confidence intervals (CIs), for two groups
#' ("Heparin" and "Standard"). It specifically requires an LME object (likely from
#' \code{nlme::lme} or similar) where the model includes a **time interaction term**
#' and the predictions are made on a **logarithmic scale** (e.g., log-transformed data).
#'
#' The function assumes a model structure where the last two fixed-effect coefficients
#' represent the slope of the reference group (Heparin) and the difference in slopes
#' between the two groups, respectively.
#'
#' @param object An LME model object (e.g., from \code{lme()}) that contains the
#'   fixed effects coefficients in \code{object$coefficients$fixed} and the
#'   fixed effects covariance matrix in \code{object$varFix}.
#' @param r A numeric vector representing the time points (e.g., days) over which
#'   predictions should be calculated. Defaults to \code{1:240}.
#' @param xname A character string for the x-axis label. Defaults to
#'   \code{"Days since randomisation"}.
#' @param yname A character string for the y-axis label. Defaults to
#'   \code{"median multiplicative change in DDimerAUC"}.
#'
#' @return The function prints and implicitly returns a \code{ggplot} object
#'   showing the two predicted trajectories (lines) and their corresponding 95%
#'   confidence bands (ribbons).
#'
#' @details
#' The predictions and confidence intervals are calculated using the following steps:
#'
#' 1.  **Coefficient Extraction**: The function specifically targets the last two
#'     fixed effects:
#'     * $\beta_{H}$ (Reference slope/Heparin slope): \code{object$coefficients$fixed[N_var-1]}
#'     * $\beta_{S-H}$ (Slope difference/Standard vs Heparin): \code{object$coefficients$fixed[N_var]}
#' 2.  **Predicted Values (Log Scale)**:
#'     * Heparin Group (H): $\text{pred}_{H} = \beta_{H} \times \text{time}$
#'     * Standard Group (S): $\text{pred}_{S} = (\beta_{H} + \beta_{S-H}) \times \text{time}$
#' 3.  **Standard Error (SE)**: The standard errors (\code{SE_H} and \code{SE_S})
#'     are calculated using the fixed effects covariance matrix (\code{object$varFix})
#'     and the standard formula for the variance of a linear combination of estimators.
#'     
#' 4.  **Transformation (Percent Change)**: The predictions and CIs are transformed
#'     from the log scale to a percentage multiplicative change scale using the
#'     formula: $\text{Percent Change} = 100 \times (\exp(\text{Prediction}_{\text{log}}) - 1)$.
#'     The 95\% CIs use the factor $1.96 \times \text{SE}$.
#'
#' **Note on Time Vector**: If \code{r} has 240 elements, the x-axis labels are
#' divided by 24 (e.g., to convert 24-hour periods to days/weeks).
#'
#' @importFrom ggplot2 ggplot aes scale_x_continuous scale_y_continuous theme
#'   element_text theme_minimal labs geom_line geom_ribbon
#' @importFrom stats qt
#' @importFrom rlang sym
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' # Assuming 'lme_model_fit' is a fitted LME object with the required structure
#'
#' # Plotting the results over 240 time points (e.g., 10 days if divided by 24)
#' plot_lme(
#'   object = lme_model_fit,
#'   r = 1:240,
#'   xname = "Time (24h intervals)",
#'   yname = "DDimer AUC Change (%)"
#' )
#'
#' # Plotting results over 10 days
#' plot_lme(
#'   object = lme_model_fit,
#'   r = 1:10,
#'   xname = "Time (Days)",
#'   yname = "Relative Change (%)"
#' )
#' }
#'
plot_lme <- function(object,r=1:240,xname="Days since randomisation",yname="median multiplicative change in DDimerAUC"){
  
  time=r
  N <- length(r)
  N_var <- dim(object$varFix)[1]
  covar1 <- object$varFix[c(N_var-1),c(N_var-1)]
  covar2 <- object$varFix[(N_var-1):N_var,(N_var-1):N_var]
  mat1 <- cbind(time)
  mat2 <- cbind(time,time)
  SE_H <- sqrt(diag(mat1%*%covar1%*%t(mat1)))
  SE_S <- sqrt(diag(mat2%*%covar2%*%t(mat2)))
  # predictions are on log scale
  pred_H <- object$coefficients$fixed[N_var-1]*time
  pred_S <- object$coefficients$fixed[N_var-1]*time+ object$coefficients$fixed[N_var]*time
  lower_H <- 100*exp(pred_H-1.96*SE_H)-100
  upper_H <- 100*exp(pred_H+1.96*SE_H)-100
  lower_S <- 100*exp(pred_S-1.96*SE_S)-100
  upper_S <- 100*exp(pred_S+1.96*SE_S)-100
  pred_H <- 100*exp(pred_H)-100
  pred_S <- 100*exp(pred_S)-100
  if(length(r)==240) r <- r/24
  df <- data.frame(time=rep(r,2),group=factor(c(rep("Heparin",N),rep("Standard",N))),pred=c(pred_H,pred_S),lower=c(lower_H,lower_S),upper=c(upper_S,upper_S))
  
  p <- ggplot(data = df, aes(x = time, y = pred, group = group,colour=group,fill=group))+scale_x_continuous(name=xname,breaks=c(0,2,4,6,8,10)) +scale_y_continuous(name=yname,limits=c(min(df$lower),max(df$upper))) + theme(axis.text=element_text(size=12),
                                                                                                                                                                                                                             axis.title=element_text(size=14,face="bold"),strip.text.x = element_text(
                                                                                                                                                                                                                               size = 12, face = "bold"
                                                                                                                                                                                                                             ))+geom_line(linetype=3,linewidth=3)+ geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.3)
  p
}