ms_llik <- function(pars,sf,model)
{
   N <- dim(sf)[1];
   Mp <- dim(sf)[2];
   M <- Mp - 1;
   if(model[1] == 'pm')
   {
      ii <- pars[2:length(pars)];
      ms_estot <- pmdlms_estot;
      checkpars <- log(pars);
   } else if(model[1] == 'pr')
   {
      ii <- pars[3:length(pars)];
      ms_estot <- prdlms_estot;
      checkpars <- log(pars);
   } else if(model[1] == 'dd')
   {
      ii <- pars[3:length(pars)];
      ms_estot <- dddlms_estot;
      checkpars <- log(c(pars[1],ii));
      if(pars[2] < -10 | pars[2] > 0.99)
      {
         warning('alpha value is larger than 1 or very negative.')
         llik <- -Inf;
         return(llik);
      }
   }
   if(any(is.na(checkpars) | is.nan(checkpars)))
   {
      warning('One or more parameter values are not numbers.')
      llik <- -Inf;
      return(llik)
   }
   if(any(checkpars < -20 | checkpars > 20))
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
      lesk <- ms_lesk(pars,qq,k,model);
      sk <- sf[cnt,Mp];
      llik <- llik + sk * lesk - lgamma(sk + 1);
   }
   return(llik);
}

pmdlms_estot <- function(pars,qq)
{
   logfun <- function(x) pmdlms_estot_int(x,pars,qq);
   lestot <- integral_peak(logfun);
   estot <- exp(lestot);
   return(estot);
}

prdlms_estot <- function(pars,qq)
{
   th <- pars[1];
   ph <- pars[2];
   ii <- pars[3:length(pars)]
   aux <- sum(ii * log(1 - qq));
   estot <- th * log((th * ph - (th + ph) * aux)/(th * ph - th * aux));
   return(estot);
}

dddlms_estot <- function(pars,qq)
{
   logfun <- function(x) dddlms_estot_int(x,pars,qq);
   lestot <- integral_peak(logfun);
   estot <- exp(lestot);
   return(estot);
}

ms_lesk <- function(pars,qq,k,model)
{
   if(model[1] == 'pm' && model[2] == 'dl')
   {
      ms_lesk_int <- pmdlms_lesk_int;
   }
   if(model[1] == 'pr' && model[2] == 'dl')
   {
      ms_lesk_int <- prdlms_lesk_int;
   }
   if(model[1] == 'dd' && model[2] == 'dl')
   {
      ms_lesk_int <- dddlms_lesk_int;
   }
   logfun <- function(x) ms_lesk_int(x,pars,qq,k);
   lesk <- integral_peak(logfun);
   return(lesk);
}

pmdlms_estot_int <- function(x,pars,qq)
{
   th <- pars[1];
   ii <- pars[2:length(pars)];
   iiqq <- -sum(ii * log(1 - qq));
   y <- rep(0,length(x));
   for(cnt in 1:length(x))
   {
      if(x[cnt] > 0)
      {
         y[cnt] <- log(-expm1(-iiqq * x[cnt])) + log(th) - th * x[cnt] - log(x[cnt]);
      } else
      {
         y[cnt] <- log(iiqq) + log(th);
      }
   }
   return(y);
}

prdlms_estot_int <- function(x,pars,qq)
{
   th <- pars[1];
   ph <- pars[2];
   ii <- pars[3:length(pars)];
   iiqq <- -sum(ii * log(1 - qq));
   y <- rep(0,length(x));
   for(cnt in 1:length(x))
   {
      if(x[cnt] == 0)
      {
         y[cnt] <- -Inf;
      } else
      {
         if(ph^2 * x[cnt]/(th + ph) > 100)
         {
            aux <- -th * ph * x[cnt]/(th + ph);
         } else
         {
            aux = -ph * x[cnt] + log(expm1(ph^2 * x[cnt]/(th + ph)));
         }
         y[cnt] <- aux + log(th) + log(1 - exp(-iiqq * x[cnt])) - log(x[cnt]);
      }
   }
   return(y);
}

dddlms_estot_int <- function(x,pars,qq)
{
   th <- pars[1];
   al <- pars[2];
   ii <- pars[3:length(pars)];
   iiqq <- -sum(ii * log(1 - qq));
   y <- rep(0,length(x));
   if(al > 0)
   {
      be <- 1/(1 - al); # exponent of substitution
      for(cnt in 1:length(x))
      {
         if(x[cnt] > 0)
         {
            y[cnt] <- (1 - al) * log(th) - lgamma(2 - al) + log(-expm1(-iiqq * x[cnt]^be)) -log(x[cnt]^be) - th * x[cnt]^be;
         } else
         {
            y[cnt] <- (1 - al) * log(th) - lgamma(2-al) + log(iiqq);
         }
      }
   } else if(al < 0)
   {
      for(cnt in 1:length(x))
      {
         if(x[cnt] > 0)
         {
            y[cnt] <- (1 - al) * log(th) -lgamma(1 - al) + log(-expm1(-iiqq*x[cnt])) - (1 + al) * log(x[cnt]) - th * x[cnt];
         } else
         {
            y[cnt] <- -Inf;
         }
      }
   } else # al=0
   {
      y <- pmdlms_estot_int(x,pars[-2],qq)
      #for(cnt in 1:length(x))
      #{
      #   if(x[cnt] > 0)
      #   {
      #      y[cnt] <- log(th) + log(-expm1(-iiqq * x[cnt])) - log(x[cnt]) - th * x[cnt];
      #   } else
      #   {
      #      y[cnt] <- log(th) + log(iiqq);
      #   }
      #}
   }
   return(y)
}

pmdlms_lesk_int <- function(x,pars,qq,k)
{
   th <- pars[1];
   ii <- pars[2:length(pars)];
   if(min(ii) < 0)
   {
      stop('ii must be non-negative.')
   }
   if(min(x) < 0)
   {
      stop('x must be non-negative.')
   }

   y <- rep(0,length(x));
   for(cnt in 1:length(x))
   {
      if(x[cnt] > 0)
      {
         y[cnt] <- sum(difgamln(ii * x[cnt],k)) - sum(lgamma(k + 1)) + sum((ii * x[cnt]) * log(1 - qq)) + sum(k * log(qq))  + log(th) - th * x[cnt] - log(x[cnt]);
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

prdlms_lesk_int <- function(x,pars,qq,k)
{
   th <- pars[1];
   ph <- pars[2];
   ii <- pars[3:length(pars)];
   y <- rep(0,length(x));
   for(cnt in 1:length(x))
   {
      if(x[cnt] == 0)
      {
         y[cnt] <- -Inf;
      } else
      {
         if(ph^2 * x[cnt]/(th + ph) > 100)
         {
            aux <- -th * ph * x[cnt]/(th + ph);
         } else
         {
            aux <- -ph * x[cnt] + log(expm1(ph^2 * x[cnt]/(th + ph)));
         }
         y[cnt] <- aux + sum(ii * log(1 - qq)) * x[cnt] + sum(difgamln(ii * x[cnt],k)) - log(x[cnt]);
      }
   }
   y <- y + log(th) + sum(k * log(qq)) - sum(lgamma(k + 1));
   return(y)
}

dddlms_lesk_int <- function(x,pars,qq,k)
{
   th <- pars[1];
   al <- pars[2];
   ii <- pars[3:length(pars)];
   iiqq <- -sum(ii * log(1 - qq));
   y <- rep(0,length(x));
   if(al > 0)
   {
      be <- 1/(1 - al); # exponent of substitution
      for(cnt in 1:length(x))
      {
         y[cnt] <- -lgamma(2 - al) - iiqq * x[cnt]^be + sum(log(ii[k > 0])) + (sum(k > 0) - 1) * be * log(x[cnt]) + sum(difgamln(ii * x[cnt]^be + 1,k - (k > 0))) - th * x[cnt]^be;
      }
   } else if(al < 0)
   {
      for(cnt in 1:length(x))
      {
         y[cnt] <- -lgamma(1 - al) - iiqq * x[cnt] + sum(log(ii[k > 0])) + (sum(k > 0) - 1 - al) * log(x[cnt]) + sum(difgamln(ii * x[cnt] + 1,k - (k > 0))) - th * x[cnt];
      }
   } else # al=0
   {
      y <- pmdlms_lesk_int(x,pars[-2],qq,k);
      return(y);
   }
   y <- y + (1 - al) * log(th) + sum(k * log(qq)) - sum(lgamma(k + 1));
   return(y);
}
