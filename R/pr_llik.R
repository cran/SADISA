pr_lesk <- function(pars,qq,k)
{
   logfun <- function(x) pr_lesk_int(x,pars,qq,k);
   lesk <- integral_peak(logfun);
   return(lesk);
}

pr_estot <- function(pars,qq)
{
   th <- pars[1];
   ph <- pars[2];
   ii <- pars[3];
   estot <- th * log((th * ph - (th + ph) * ii * log(1 - qq))/(th * ph - th * ii * log(1 - qq)));
   return(estot);
}

pr_lesk_int <- function(x,pars,qq,k)
{
   th <- pars[1];
   ph <- pars[2];
   ii <- pars[3];
   y <- rep(0,length(x));
   for(cnt in 1:length(x))
   {
      if(x[cnt] == 0)
      {
         y[cnt] = -Inf;
      } else
      {
         if(k > 1)
         {
            y[cnt] <- sum(log(ii * x[cnt] + (1:(k - 1))));
         } else
         {
            y[cnt] <- 0;
         }
         y[cnt] <- y[cnt] + ii * x[cnt] * log(1 - qq) + log(ii) + log(th);
         if(ph^2 * x[cnt]/(th + ph) > 100)
         {
             y[cnt] <- y[cnt] - th * ph * x[cnt]/(th + ph);
         } else
         {
             y[cnt] <- y[cnt] - ph * x[cnt] + log(expm1(ph^2 * x[cnt]/(th + ph)));
         }
      }
   }
   y <- y + k * log(qq) - lgamma(k + 1);
   return(y);
}
