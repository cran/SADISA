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
#' @references Haegeman, B. & R.S. Etienne (2017). A general sampling formula for community structure data. Methods in Ecology & Evolution. In press.
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
   intfun <- function(x)
   {
      #if(any(is.nan(logfun(exp(x)))))
      #{
      #   print(exp(x))
      #   print(logfun(exp(x)))
      #   print(logfun)
      #}
      return(exp((x + logfun(exp(x))) - ymax))
   }
   intres <- stats::integrate(f = intfun, lower = xlft, upper = xrgt, rel.tol = 1e-10, abs.tol = 1e-10);
   corrfact <- intres$value;
   logQ <- ymax + log(corrfact);
   return(logQ);
}

#' @title Converts different formats to represent multiple sample data
#' @description Converts the full abundance matrix into species frequencies
#' If S is the number of species and M is the number of samples, then fa is
#' the full abundance matrix of dimension S by M. The  for example
#' fa = [0 1 0;3 2 1;0 1 0] leads to sf = [0 1 0 2;3 2 1 1];
#' @param fa the full abundance matrix with species in rows and samples in columns
#' @return the sample frequency matrix
#' @references Haegeman, B. & R.S. Etienne (2017). A general sampling formula for community structure data. Methods in Ecology & Evolution. In press.
#' @export

convert_fa2sf <- function(fa)
{
   dfa <- dim(fa);
   S <- dfa[1];
   M <- dfa[2];
   aux <- unique(fa);
   da <- dim(aux)[1]
   freq <- rep(0,da);
   for(cnt in 1:da)
   {
      ref <- t(t(rep(1,S))) %*% (aux[cnt,]);
      eqs <- rowSums(fa == ref);
      freq[cnt] = sum(eqs == M);
   }
   sf <- cbind(aux,freq);
   colnames(sf) <- NULL;
   return(sf)
}

convert_sf2fa <- function(sf)
{
   notzero <- which(sf != 0)
   fa <- NULL
   for(i in length(notzero):1)
   {
      fa <- c(fa,rep(notzero[i],sf[notzero[i]]))
   }
   return(fa)
}

difgamln <- function(a,n)
{
# computes gammaln(a+n)-gammaln(a) as a sum of logarithms
# for large a the result is numerically more precise than
# directly evaluating the difference of gammaln
# . a = real number, or vector of real numbers
# . n = natural number, or vector of natural numbers
   if(length(a) != length(n))
   {
      stop('Input arguments should have same length.')
   }
   if(min(a) < 0)
   {
      stop('The first argument of difgamln must be non-negative.')
   }
   if(min(n) < 0)
   {
      stop('The second argument of difgamln must be non-negative.')
   }
   cvec <- rep(0,length(a));
   for(ctr in 1:length(a))
   {
      if(a[ctr] < 1e5 && n[ctr] > 100)
      {
         cvec[ctr] <- lgamma(a[ctr] + n[ctr]) - lgamma(a[ctr]);
      } else
      {
         if(n[ctr] == 0)
         {
            cvec[ctr] <- 0;
         } else #n[ctr] > 0
         {
            cvec[ctr] <- sum(log(a[ctr] + (0:(n[ctr] - 1))));
         }
         #if(is.nan(cvec[ctr])) {print(a[ctr]); print(n[ctr])}
         #if(cvec[ctr] == -Inf) {print(a[ctr]); print(n[ctr])}
      }
   }
   return(cvec)
}

checkiijj <- function(ii,jj,mult)
{
   if(is.list(ii))
   {
      if(is.list(ii[[1]]))
      {
         stop('I parameters not correctly specified.')
      }
      ii <- unlist(ii);
   }
   if(is.list(jj))
   {
      if(is.list(jj[[1]]))
      {
         stop('Sample sizes not correctly specified.')
      }
      jj <- unlist(jj);
   }
   if(length(ii) <= 1 && mult != 'single')
   {
      stop('You need to specify more than one I value.')
   }
   if(length(jj) <= 1 && mult != 'single')
   {
      stop('You need to specify more than one sample size.')
   }
   if(length(ii) > 1 && mult == 'single')
   {
      stop('You need to specify only one I value.')
   }
   if(length(jj) > 1 && mult == 'single')
   {
      stop('You need to specify only one sample size.')
   }
   return(list(ii = ii,jj = jj))
}

dec2bin <- function(y,ly)
{
   stopifnot(length(y) == 1, mode(y) == 'numeric')
   q1 <- (y / 2) %/% 1
   r <- y - q1 * 2
   res <- c(r)
   while(q1 >= 1)
   {
      q2 <- (q1 / 2) %/% 1
      r <- q1 - q2 * 2
      q1 <- q2
      res <- c(r, res)
   }
   res <- c(rep(0,ly - length(res)),res)
   return(res)
}

dec2binmat <- function(y)
{
   numrows <- 2^y
   res <- matrix(0,numrows,y)
   for(i in 0:(numrows-1))
   {
      res[i + 1,] <- dec2bin(i,y)
   }
   return(res)
}

is.nonnegativewholenumber <- function(x, tol = .Machine$double.eps^0.5)
{
   ans <- abs(x - DDD::roundn(x)) < tol & x >= 0
   if(any(is.na(ans)))
   {
      ans <- FALSE
   }
   return(all(ans))
}

logminexpm1approx <- function(iiqq,xcnt,be)
{
   if(abs(-iiqq * xcnt^be) < 1E-100)
   {
      ans <- log(iiqq) + be * log(xcnt)
   } else
   {
      ans <- log(-expm1(-iiqq * xcnt^be))
   }
   return(ans)
}
