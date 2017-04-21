mdd_lesk <- function(pars,qq,k)
{
   logfun <- function(x) mdd_lesk_int(x,pars,qq,k)
   lesk <- integral_peak(logfun);
   return(lesk);
}

mdd_estot <- function(pars,qq)
{
   logfun <- function(x) mdd_lestot_int(x,pars,qq)
   lestot <- integral_peak(logfun);
   estot <- exp(lestot);
   return(estot);
}

mdd_lesk_int <- function(x,pars,qq,k)
{
   th <- pars[1];
   al <- pars[2];
   ii <- pars[3];
   y <- rep(0,length(x));
   if(al > 0)
   {
      be <- 1/(1 - al) #exponent of substitution
      for(cnt in 1:length(x))
      {
         if(k > 1)
         {
            y[cnt] <- sum(log(ii * x[cnt]^be + (1:(k - 1))));
         } else
         {
            y[cnt] <- 0;
         }
         y[cnt] <- y[cnt] + (1 - al) * log(th) - lgamma(2 - al) + k * log(qq) - lgamma(k + 1) + ii * x[cnt]^be * log(1 - qq) + log(ii) - th * x[cnt]^be;
      }
   } else
   {
      for(cnt in 1:length(x))
      {
         if(k > 1)
         {
            y[cnt] <- sum(log(ii * x[cnt] + (1:(k - 1))));
         } else
         {
            y[cnt] <- 0;
         }
         y[cnt] <- y[cnt] + (1 - al) * log(th) - lgamma(1 - al) + k * log(qq) - lgamma(k + 1) + ii * x[cnt] * log(1 - qq) + log(ii) - al * log(x[cnt]) - th * x[cnt];
      }
   }
   return(y);
}

mdd_lestot_int <- function(x,pars,qq)
{
   th <- pars[1];
   al <- pars[2];
   ii <- pars[3];
   y <- rep(0,length(x));
   if(al > 0)
   {
      be <- 1/(1 - al) #exponent of substitution
      for(cnt in 1:length(x))
      {
         if(x[cnt] > 0)
         {
            y[cnt] <- (1 - al) * log(th) - lgamma(2 - al) + log(-expm1(ii * log(1 - qq) * x[cnt]^be)) - log(x[cnt]^be) - th * x[cnt]^be;
         } else
         {
            y[cnt] <- (1 - al) * log(th) - lgamma(2 - al) + log(-ii * log(1 - qq));
         }
      }
   } else
   if(al < 0)
   {
      for(cnt in 1:length(x))
      {
         if(x[cnt] > 0)
         {
            y[cnt] <- (1 - al) * log(th) - lgamma(1 - al) + log(-expm1(ii * log(1 - qq) * x[cnt])) - (1 + al) * log(x[cnt]) - th * x[cnt];
         } else
         {
            y[cnt] <- -Inf;
         }
      }
   } else
   {
      for(cnt in 1:length(x))
      {
         if(x[cnt] > 0)
         {
            y[cnt] <- log(th) + log(-expm1(ii * log(1 - qq) * x[cnt])) - log(x[cnt]) - th * x[cnt];
         } else
         {
            y[cnt] <- log(th) + log(-ii * log(1 - qq));
         }
      }
   }
   return(y);
}
