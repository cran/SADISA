pm_estot <- function(pars,qq)
{
   th <- pars[1];
   ii <- pars[2];
   estot <- th * log(th - ii * log(1 - qq)) - th * log(th);
   return(estot);
}

pm_lesk <- function(pars,qq,k)
{
   logfun <- function(x) pm_lesk_int(x,pars,qq,k);
   lesk <- integral_peak(logfun);
   return(lesk);
}

pm_lesk_int <- function(x,pars,qq,k)
{
   th <- pars[1];
   ii <- pars[2];
   y <- rep(0,length(x));
   for(cnt in 1:length(x))
   {
      if(k > 1)
      {
         y[cnt] <- sum(log(ii * x[cnt] + (1:(k - 1))));
      } else
      {
         y[cnt] <- 0;
      }
      y[cnt] <- y[cnt] + log(ii) + ii * x[cnt] * log(1 - qq) - th * x[cnt];
   }
   y <- y + k * log(qq) - lgamma(k + 1) + log(th);
   return(y);
}
