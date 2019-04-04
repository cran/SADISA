#' @title Simulates species abundance data
#' @description Simulates species abundance data using the independent-species approach
#' @param th the theta parameter
#' @param ii the I parameters. When a vector, it is assumed that each value describes a sample or a guild
#' depending on whether mult == 'ms' or mult == 'mg'. When mult = 'both', a list of lists must be specified, with
#' each list element relates to a sample and contains a list of values across guilds.
#' @param jj the sample sizes for each sample and each guild. Must have the same structure as ii
#' @param mult When set to 'single', the loglikelihood of a single abundance vector will be computed
#' When set to 'mg' the loglikelihood for multiple guilds is computed.
#' When set to 'ms' the loglikelihood for multiple samples from the same metacommunity is computed.
#' When set to 'both' the loglikelihood for multiple guilds within multiple samples is computed.
#' @return abund abundance vector, a list of abundance vectors, or a list of lists of abundance vectors,
#' or a list of lists of lists of abundance vectors
#' The first layer of the lists corresponds to different simulations
#' When mult = 'mg', each list contains a list of abundance vectors for different guilds.
#' When mult = 'ms', each list contains a list of abundance vectors for different samples
#' from the same metacommunity. In this case the vectors should have equal lengths and may contain
#' zeros because there may be species that occur in multiple samples and species that do not occur
#' in some of the samples.
#' When mult = 'both', each list will be a list of lists of multiple guilds within a sample
#' @param nsim Number of simulations to perform
#' @details Not all combinations of metacommunity model and local community model have been implemented yet.
#' because this requires checking for numerical stability of the integration. The currently available model combinations are c('pm','dl').
#' @keywords model species-abundance-distribution
#' @references Haegeman, B. & R.S. Etienne (2017). A general sampling formula for community structure data. Methods in Ecology & Evolution 8: 1506-1519. doi: 10.1111/2041-210X.12807
#' @export
pmdl_sim <- function(th,ii,jj,mult = 'single',nsim = 1)
{
   if(!is.numeric(th) | length(th) > 1)
   {
      stop('th should be a single scalar.')
   }
   out <- list();
   for(i in 1:nsim)
   {
      out[[i]] <- list();
      if(mult == 'single')
      {
         ff <- checkiijj(ii,jj,mult); ii <- ff$ii; jj <- ff$jj; rm(ff);
         qq <- jj/(ii + jj);
         out[[i]] <- pmdlsing_sim(th,ii,qq);
      } else if(mult == 'ms')
      {
         ff <- checkiijj(ii,jj,mult); ii <- ff$ii; jj <- ff$jj; rm(ff);
         qq <- jj/(ii + jj);
         out[[i]] <- pmdlmult_sim(th,ii,qq);
      } else if(mult == 'mg')
      {
         loopguilds <- length(ii)
         ff <- checkiijj(ii,jj,mult); ii <- ff$ii; jj <- ff$jj; rm(ff);
         qq <- jj/(ii + jj);
         for(g in 1:loopguilds)
         {
            out[[i]][[g]] <- pmdlsing_sim(th,ii[g],qq[g]);
         }
      } else if(mult == 'both')
      {
         if(!is.list(ii) && !is.list(ii[[1]]))
         {
            stop('The I parameters must be in a list of lists.')
         }
         if(!is.list(jj) && !is.list(jj[[1]]))
         {
            stop('The sample sizes must be in a list of lists.')
         }
         numsam <- length(ii);
         if(numsam != length(jj))
         {
            stop('The dimensions of the ii list and the jj list do not match.')
         }
         loopguilds <- length(ii[[1]]);
         for(g in 1:loopguilds)
         {
            ii2 <- NULL;
            jj2 <- NULL;
            for(s in 1:numsam)
            {
              ii2 <- c(ii2,ii[[s]][[g]]);
              jj2 <- c(jj2,jj[[s]][[g]]);
            }
            qq <- jj2/(ii2 + jj2);
            out[[i]][[g]] <- pmdlmult_sim(th,ii2,qq);
         }
         out2 <- out;
         out[[i]] <- list();
         for(s in 1:numsam)
         {
            out[[i]][[s]] <- list();
            for(g in 1:loopguilds)
            {
               out[[i]][[s]][[g]] <- out2[[i]][[g]][[s]]
            }
         }
      } else
      {
         stop('"mult" should be either "single", "ms", "mg", or "both"')
      }
   }
   if(nsim == 1)
   {
      out <- out[[1]]
   }
   return(out);
}

pmdlsing_sim <- function(th,ii,qq)
{
# simulation of sample abundance data
# pm: point mutation speciation model
# dl: dispersal-limited sampling model
# . th = metacomm diversity (scalar)
# . ii = coupling strength (scalar)
# . qq = sampling effort (scalar)
# . nsim = number of simulations

# compute species abundance distribution
  jj <- max(c(10,DDD::roundn(10 * qq * ii/(1 - qq))));
  k1 <- unique(DDD::roundn(10^(seq(0,log10(jj),length.out = 1000))));
  aux <- rep(0,length(k1));
  for(cnt in 1:length(k1))
  {
     aux[cnt] <- pmdlmu_esk(th,ii,qq,k1[cnt]);
  }
  k2 <- 1:jj;
  kesk <- stats::spline(x = log2(k1),y = k1 * aux,xout = log2(k2))$y;
  esk <- kesk/k2;
  esk[which(esk < .Machine$double.eps)] <- .Machine$double.eps;

  # sample community abundances
  simsk <- stats::rpois(length(esk),esk);
  simfa <- convert_sf2fa(simsk);
  return(simfa);
}

pmdlmu_esk <- function(th,ii,qq,k)
{
  ff <- function(x)
  {
     return(exp(pmdlmu_esk_int(x,th,ii,qq,k)))
  }
  esk <- stats::integrate(f = ff,lower = 0,upper = 1,rel.tol = 1e-9,abs.tol = 0)$value;
  return(esk)
}

pmdlmu_esk_int <- function(xx,th,ii,qq,k)
{
   yy <- rep(0,length(xx));
   for(cnt in 1:length(xx))
   {
      x <- xx[cnt];
      if(x == 0)
      {
         if(sum(k > 0) == 1)
         {
            kgt0 <- which(k > 0);
            y <- log(ii[kgt0]) + k[kgt0] * log(qq[kgt0]) - log(k[kgt0]) + log(th);
         } else
         {
            y <- -Inf;
         }
      } else
      {
         y <- sum(lgamma(ii * x + k)) - sum(lgamma(ii * x)) - sum(lgamma(k + 1)) + sum((ii * x) * log(1 - qq)) + sum(k * log(qq)) + log(th) - th * x - log(x);
      }
      yy[cnt] <- y;
   }
   return(yy);
}

pmdlmult_sim <- function(th,ii,qq)
{
   # simulation of sample abundance data
   # pm: point mutation speciation model
   # dl: dispersal-limited sampling model
   # mu: multiple samples (from one metacomm)
   # . th = metacomm diversity (scalar)
   # . ii = coupling strength (vector of length M)
   # . qq = sampling effort (vector of length M)

   # sample metacommunity abundances
   ff <- function(x)
   {
      return(exp(pm_estot_int(x,th,ii,qq)))
   };
   es0 <- stats::integrate(f = ff,lower = 0,upper = 1,rel.tol = 1e-9,abs.tol = 0)$value;
   nx <- stats::rpois(n = length(es0),es0);
   xs <- rep(0,nx);
   for(cnt in 1:nx)
   {
      prob <- stats::runif(1);
      ff2 <- function(x)
      {
         intresult <- rep(0,length(x))
         for(i in 1:length(x))
         {
            intresult[i] <- stats::integrate(ff,0,x[i])$value - prob * es0
         }
         return(intresult)
      }
      xs[cnt] <- stats::uniroot(f = ff2,interval = c(0,1),tol = .Machine$double.eps)$root;
   }

   # sample local community abundances
   M <- length(ii);
   fasim <- matrix(0,nrow = nx, ncol = M);
   for(ctr in 1:nx)
   {
      nn <- ii * xs[ctr]
      kk <- stats::rnbinom(n = nn,size = nn, prob = 1 - qq);
      while(sum(kk) == 0)
      {
         kk <- stats::rnbinom(n = nn,size = nn, prob = 1 - qq);
      }
      fasim[ctr,] <- kk;
   }
   #sfsim <- convert_fa2sf(fasim);
   fasim <- t(fasim);
   ff <- list();
   for(sam in 1:M)
   {
      ff[[sam]] <- fasim[sam,]
   }
   fasim <- ff;
   rm(ff);
   return(fasim);
}

pmdlmu_es0_int <- function(xx,th,ii,qq)
{
  al <- -sum(ii * log(1 - qq));
  yy <- rep(0,length(xx));
  for(cnt in 1:length(xx))
  {
     x <- xx[cnt];
     if(x == 0)
     {
        y <- log(al) + log(th);
     } else {
        y <- log(1 - exp(-al * x)) + log(th) - th * x - log(x);
     }
     yy[cnt] <- y;
  }
  return(yy);
}
