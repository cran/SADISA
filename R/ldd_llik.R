ldd_ejtot <- function(pars,qq)
{
   logfun <- function(x) ldd_ejtot_int(x,pars,qq);
   lejtot <- integral_peak(logfun);
   ejtot <- exp(lejtot);
   return(ejtot);
}

ldd_ejtot_int <- function(x,pars,qq)
{
   th <- pars[1];
   al <- pars[2];
   ii <- pars[3];
   y <- rep(0,length(x));
   for(cnt in 1:length(x))
   {
      if(ii * x[cnt] > al)
      {
         y[cnt] <- log(th) - th * x[cnt] + log(ii) + log(ii * x[cnt] - al) - log(ii * x[cnt] - (1 - qq)^(ii * x[cnt] - al) * al);
      } else
      {
         y[cnt] <- log(th) - th * x[cnt] + log(ii) + (al - ii*x[cnt]) * log(1 - qq) + log(al - ii * x[cnt]) - log(al - (1 - qq)^(-ii * x[cnt] + al) * ii * x[cnt]);
      }
   }
   y <- y + log(qq) - log(1 - qq);
   return(y);
}

ldd_lesk <- function(pars,qq,k)
{
   logfun <- function(x) ldd_lesk_int(x,pars,qq,k);
   lesk <- integral_peak(logfun);
   return(lesk);
}

ldd_estot <- function(pars,qq)
{
   logfun <- function(x) ldd_lestot_int(x,pars,qq);
   lestot <- integral_peak(logfun);
   estot <- exp(lestot);
   return(estot);
}

ldd_lesk_int <- function(x,pars,qq,k)
{
   th <- pars[1];
   al <- pars[2];
   ii <- pars[3];
   y <- rep(0,length(x));
   for(cnt in 1:length(x))
   {
      if(ii * x[cnt] > al)
      {
         y[cnt] = log(th) - th*x[cnt] + log(ii) + sum(log(ii * x[cnt]-al + (0:(k - 1)))) - (-ii * x[cnt] + al) * log(1 - qq) - log(ii * x[cnt] - (1 - qq)^(ii * x[cnt] - al) * al);
      } else
      {
         if(k > 1)
         {
            y[cnt] <- sum(log(ii * x[cnt] - al + (1:(k - 1))));
         } else
         {
            y[cnt] <- 0;
         }
         y[cnt] <- y[cnt] + log(th) - th * x[cnt] + log(ii) + log(al - ii * x[cnt]) - log(al - (1 - qq)^(-ii * x[cnt] + al) * ii * x[cnt]);
      }
   }
   y <- y + k * log(qq) - lgamma(k + 1);
   return(y);
}

ldd_lestot_int <- function(x,pars,qq)
{
   th <- pars[1];
   al <- pars[2];
   ii <- pars[3];
   y <- rep(0,length(x));
   for(cnt in 1:length(x))
   {
      if(ii * x[cnt] > al)
      {
         y[cnt] <- log(th) - th * x[cnt] + log(ii) + log(1 - (1 - qq)^(ii * x[cnt] - al)) - log(ii * x[cnt] - (1 - qq)^(ii * x[cnt] - al) * al);
      } else
      {
         y[cnt] <- log(th) - th * x[cnt] + log(ii) + log(1 - (1 - qq)^(-ii * x[cnt] + al)) - log(al - (1 - qq)^(-ii * x[cnt] + al) * ii * x[cnt]);
      }
   }
   return(y);
}
