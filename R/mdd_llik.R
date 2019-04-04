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
            #if(any(is.nan(y))) { print(k); print(x); print(y)}
         } else
         {
            y[cnt] <- 0;
         }
         y[cnt] <- y[cnt] + ii * x[cnt]^be * log(1 - qq) - th * x[cnt]^be;
      }
      y <- y - lgamma(2 - al);
   } else # al <= 0
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
         y[cnt] <- y[cnt] + ii * x[cnt] * log(1 - qq) - al * log(x[cnt]) - th * x[cnt];
      }
      y <- y - lgamma(1 - al);
   }
   y <- y + (1 - al) * log(th) + k * log(qq) - lgamma(k + 1) + log(ii)
   return(y);
}

mdd_lestot_int <- function(x,pars,qq)
{
   th <- pars[1];
   al <- pars[2];
   ii <- pars[3:length(pars)];
   iiqq <- -sum(ii * log(1 - qq));
   y <- rep(0,length(x));
   if(al > 0)
   {
      be <- 1/(1 - al) #exponent of substitution
      for(cnt in 1:length(x))
      {
         if(x[cnt] > 0)
         {
            y[cnt] <- logminexpm1approx(iiqq,x[cnt],be) - be * log(x[cnt]) - th * x[cnt]^be;
            #y[cnt] <- log(-expm1(-iiqq * xcnt^be)) - be * log(x[cnt]) - th * x[cnt]^be;
         } else
         {
            y[cnt] <- log(iiqq);
         }
      }
      y <- y + (1 - al) * log(th) - lgamma(2 - al);
   } else
   if(al < 0)
   {
      for(cnt in 1:length(x))
      {
         if(x[cnt] > 0)
         {
            y[cnt] <- logminexpm1approx(iiqq,x[cnt],1) - (1 + al) * log(x[cnt]) - th * x[cnt];
           #y[cnt] <- log(-expm1(-iiqq * x[cnt])) - (1 + al) * log(x[cnt]) - th * x[cnt];
         } else
         {
            y[cnt] <- -Inf;
         }
      }
      y <- y + (1 - al) * log(th) - lgamma(1 - al);
   } else #al = 0 (pm)
   {
      for(cnt in 1:length(x))
      {
         if(x[cnt] > 0)
         {
            y[cnt] <- logminexpm1approx(iiqq,x[cnt],1) - log(x[cnt]) - th * x[cnt];
            #y[cnt] <- log(-expm1(-iiqq * x[cnt])) - log(x[cnt]) - th * x[cnt];
         } else
         {
            y[cnt] <- log(iiqq);
         }
      }
      y <- y + log(th);
   }
   return(y);
}
