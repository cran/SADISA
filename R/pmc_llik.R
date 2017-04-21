pmc_lprj <- function(pars,qq,k)
{
   logfun <- function(x) pmc_prj_int(x,pars,qq,k);
   lprj <- integral_peak(logfun);
   return(lprj);
}

pmc_prj_int <- function(x,pars,qq,k)
{
   th = pars[1];
   ii = pars[2];
   y <- rep(0,length(x));
   for(cnt in 1:length(x))
   {
      if(k > 1)
      {
         y[cnt] <- sum(log(ii * x[cnt] + (1:(k - 1))))
      } else
      {
         y[cnt] <- 0
      }

      y[cnt] <- y[cnt] + log(ii) + ii * x[cnt] * log(1 - qq) + th * log(x[cnt]) - th * x[cnt];
   }
   y <- y + k * log(qq) - lgamma(k + 1) + th * log(th) - lgamma(th);
   return(y);
}
