#' @title Computes integral of a very peaked function
#' @description   # computes the logarithm of the integral of exp(logfun) from 0 to Inf under the following assumptions:
# . exp(logfun) has a single, sharply peaked maximum
# . exp(logfun) is increasing to the left of the peak and decreasing to the right of the peak
# . exp(logfun) can be zero or positive at zero
# . exp(logfun) tends to zero at infinity
#' @param logfun the logarithm of the function to integrate
#' @param xx the initial set of points on which to evaluate the function
#' @param xcutoff when the maximum has been found among the xx, this parameter sets the width of the interval to find the maximum in
#' @param ycutoff set the threshold below which (on a log scale) the function is deemed negligible, i.e. that it does not contribute to the integral)
#' @param ymaxthreshold sets the deviation allowed in finding the maximum among the xx
#' @return the result of the integration
#' @references Haegeman, B. & Etienne, R.S. (2016). A general sampling formula for community abundance data. Methods in Ecology & Evolution. In review.
#' @export

integral_peak <- function(logfun, xx = seq(-100,10,2), xcutoff = 2, ycutoff = 40, ymaxthreshold = 1E-12)
{
   # 1/ determine integrand peak
   yy <- xx + logfun(exp(xx));
   yy[which(is.na(yy) | is.nan(yy))] <- -Inf;
   yymax <- max(yy);
   if(yymax == -Inf)
   {
      logQ <- -Inf;
      return(logQ);
   }
   iimax <- which(yy >= (yymax - ymaxthreshold));
   xlft <- xx[iimax[1]] - xcutoff;
   xrgt <- xx[iimax[length(iimax)]] + xcutoff;
   optfun <- function(x) x + logfun(exp(x));
   optres <- stats::optimize(f = optfun, interval = c(xlft,xrgt), maximum = TRUE, tol = 1e-10);
   xmax <- optres$maximum;
   ymax <- optres$objective;

   # 2/ determine peak width
   iilft <- which((xx < xmax) & (yy < (ymax - ycutoff)));
   if(length(iilft) == 0)
   {
      xlft <- xx[1] - xcutoff;
   } else
   {
      ilft <- iilft[length(iilft)];
      xlft <- xx[ilft];
   }
   iirgt <- which((xx > xmax) & (yy < (ymax - ycutoff)));
   if(length(iirgt) == 0)
   {
      xrgt <- xx[length(xx)] + xcutoff;
   } else
   {
      irgt <- iirgt[1];
      xrgt <- xx[irgt];
   }

   # 3/ compute integral
   intfun <- function(x) exp((x + logfun(exp(x))) - ymax);
   intres <- stats::integrate(f = intfun, lower = xlft, upper = xrgt, rel.tol = 1e-10, abs.tol = 1e-10);
   corrfact <- intres$value;
   logQ <- ymax + log(corrfact);
   return(logQ);
}
