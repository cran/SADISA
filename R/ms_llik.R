ms_llik <- function(pars,sf)
{
   N <- dim(sf)[1];
   Mp <- dim(sf)[2];
   M <- Mp - 1;
   th <- pars[1];
   ii <- pars[2:length(pars)];
   lth <- log(th);
   lii <- log(ii);
   if(lth < -20 || lth > 20 || sum(lii < -20) || sum(lii > 20))
   {
      llik <- -Inf;
      return(llik);
   }
   jj <- t(sf[,M + 1]) %*% sf[,1:M];
   qq <- jj/(ii + jj);
   estot <- ms_estot(pars,qq);
   llik <- -estot;
   for(cnt in 1:N)
   {
      k <- sf[cnt,1:M];
      lesk <- ms_lesk(pars,qq,k);
      sk <- sf[cnt,Mp];
      llik <- llik + sk * lesk - lgamma(sk + 1);
   }
   return(llik);
}

ms_estot <- function(pars,qq)
{
   logfun <- function(x) ms_estot_int(x,pars,qq);
   lestot <- integral_peak(logfun);
   estot <- exp(lestot);
   return(estot);
}

ms_lesk <- function(pars,qq,k)
{
   logfun <- function(x) ms_lesk_int(x,pars,qq,k);
   lesk <- integral_peak(logfun);
   return(lesk);
}

ms_estot_int <- function(x,pars,qq)
{
   th <- pars[1];
   ii <- pars[2:length(pars)];
   al <- -sum(ii * log(1 - qq));
   y <- rep(0,length(x));
   for(cnt in 1:length(x))
   {
      if(x[cnt] > 0)
      {
         y[cnt] <- log(1 - exp(-al * x[cnt])) + log(th) - th * x[cnt] - log(x[cnt]);
      } else
      {
         y[cnt] <- log(al) + log(th);
      }
   }
   return(y);
}

ms_lesk_int <- function(x,pars,qq,k)
{
   th <- pars[1];
   ii <- pars[2:length(pars)];
   y <- rep(0,length(x));
   for(cnt in 1:length(x))
   {
      if(x[cnt] > 0)
      {
         y[cnt] <- sum(lgamma(ii * x[cnt] + k)) - sum(lgamma(ii * x[cnt])) - sum(lgamma(k + 1)) + sum((ii * x[cnt]) * log(1 - qq)) + sum(k * log(qq))  + log(th) - th * x[cnt] - log(x[cnt]);
      } else
      {
         if(sum(k > 0) == 1)
         {
            y[cnt] <- log(ii[k > 0]) + k[k > 0] * log(qq[k > 0]) - log(k[k > 0]) + log(th);
         } else
         {
            y[cnt] <- -Inf;
         }
      }
   }
   return(y);
}
