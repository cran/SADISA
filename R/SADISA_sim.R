#' @title Simulates species abundance data
#' @description Simulates species abundance data using the independent-species approach
#' @param parsmc The model parameters. For the point mutation (pm) model this is theta and I. For the protracted
#' model (pr) this is theta, phi and I. For the density-dependent model (dd) - which can also be interpreted as
#' the per-species speciation model, this is theta and alpha.
#' @param ii The I parameter. When I is a vector, it is assumed that each value
#' describes a sample or a guild depending on whether mult == 'ms' or mult == 'mg'. When mult = 'both',
#' a list of lists must be specified, with each list element relates to a sample and contains a list of values across guilds.
#' @param jj the sample sizes for each sample and each guild. Must have the same structure as ii
#' @param model the chosen combination of metacommunity model and local community model
#' as a vector, e.g. c('pm','dl') for a model with point mutation in the metacommunity and
#' dispersal limitation.
#' The choices for the metacommunity model are: 'pm' (point mutation), 'rf' (random fission),
#' 'pr' (protracted speciation), 'dd' (density-dependence).
#' The choices for the local community model are: 'dl' (dispersal limitation), 'dd' (density-dependence).
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
SADISA_sim <- function(parsmc,ii,jj,model = c('pm','dl'),mult = 'single',nsim = 1)
{
   if(!((model[1] == 'pm' | model[1] == 'pr' | model[1] == 'dd') & model[2] == 'dl'))
   {
      stop('Simulations for this metacommunity model is not implemented yet.');
   }
   if(!is.numeric(parsmc))
   {
      stop('Parameters should be scalars.')
   }
   out <- list();
   for(i in 1:nsim)
   {
      out[[i]] <- list();
      if(mult == 'single')
      {
         ff <- checkiijj(ii,jj,mult); ii <- ff$ii; jj <- ff$jj; rm(ff);
         qq <- jj/(ii + jj);
         out[[i]] <- ms_sim(parsmc,ii,qq,model);
      } else if(mult == 'ms')
      {
         ff <- checkiijj(ii,jj,mult); ii <- ff$ii; jj <- ff$jj; rm(ff);
         qq <- jj/(ii + jj);
         out[[i]] <- ms_sim(parsmc,ii,qq,model)
      } else if(mult == 'mg')
      {
         loopguilds <- length(ii)
         ff <- checkiijj(ii,jj,mult); ii <- ff$ii; jj <- ff$jj; rm(ff);
         qq <- jj/(ii + jj);
         for(g in 1:loopguilds)
         {
            out[[i]][[g]] <- ms_sim(parsmc,ii[g],qq[g],model);
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
            #cat(paste('Simulation ', i,' and guild ',g,sep = '')); cat('\n');
            out[[i]][[g]] <- ms_sim(parsmc,ii2,qq,model);
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

ms_sim <- function(parsmc,ii,qq,model)
{
   if(model[1] == 'pm')
   {
      if(parsmc < 0 | min(unlist(ii)) < 0) stop('Parameters should be positive');
      model_es0_int <- pm_estot_int;
   } else
   if(model[1] == 'pr')
   {
      if(parsmc < 0 | min(unlist(ii)) < 0) stop('Parameters should be positive');
      model_es0_int <- pr_estot_int;
   } else
   if(model[1] == 'dd')
   {
      if(parsmc[1] < 0 | parsmc[2] > 1 | min(unlist(ii)) < 0) stop('Parameters, except alpha, should be positive; alpha cannot be larger than 1.');
      model_es0_int <- dd_estot_int;
   }

   # sample metacommunity abundances
   ff <- function(x)
   {
      return(exp(model_es0_int(x,parsmc,ii,qq)))
   }
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
            intresult[i] <- stats::integrate(ff,0,x[i])$value - prob * es0;
         }
         return(intresult);
      }
      #cat(paste('Compute xs[',cnt,']\n',sep = ''))
      xs[cnt] <- stats::uniroot(f = ff2,interval = c(0,1),tol = .Machine$double.eps)$root;
   }
   if(model[1] == 'dd' && parsmc[2] > 0) # this is alpha
   {
      xs <- xs^(1/(1 - parsmc[2]));
   }

   # sample local community abundances
   M <- length(ii);
   fasim <- matrix(0,nrow = nx, ncol = M);

   configs <- dec2binmat(M)
   for(ctr in 1:nx)
   {
      nn <- ii * xs[ctr]
      kk <- stats::rnbinom(n = M,size = nn, prob = 1 - qq);
      fasim[ctr,] <- kk;
      if(sum(kk) == 0)
      {
         p0 <- stats::dnbinom(x = 0, size = nn, prob = 1 - qq)
         p <- rep(0,2^M - 1)
         for(i in 2:(2^M))
         {
            p[i - 1] <- prod(p0^(1 - configs[i,]) * (1 - p0)^configs[i,])
         }
         config_possible <- which(c(0,p) != 0)
         if(length(config_possible) == 0)
         {
            totsam <- rowSums(configs)
            config_possible <- which(totsam == 1)
            p <- rep(1,length(config_possible))
         } else
         {
            p <- p[config_possible - 1]
         }
         config_sampled <- configs[DDD::sample2(config_possible,size = 1,prob = p),]
         samples_present <- which(config_sampled == 1)
         for(i in samples_present)
         {
            fele <- stats::runif(n = 1, min = stats::dnbinom(0, size = nn[i], prob = 1 - qq[i]), max = 1)
            fasim[ctr,i] <- stats::qnbinom(fele, size = nn[i], prob = 1 - qq[i])
            if(fasim[ctr,i] <= 0 | fasim[ctr,i] == Inf | is.na(fasim[ctr,i]) | is.nan(fasim[ctr,i]))
            {
               nmax <- 1000000
               probs <- stats::dnbinom(1:nmax, size = nn[i], prob = 1 - qq[i])
               which0 <- which(probs <= 0)
               if(length(which0) > 0)
               {
                  nmax <- min(which0) - 1
               }
               fasim[ctr,i] <- DDD::sample2(x = 1:nmax, size = 1, prob = probs[1:nmax])
            }
         }
      }
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

pm_estot_int <- function(x,parsmc,ii,qq)
{
   return(pmdlms_estot_int(x,c(parsmc,ii),qq));
}

dd_estot_int <- function(x,parsmc,ii,qq)
{
   return(mdd_lestot_int(x,c(parsmc,ii),qq))
}

pr_estot_int <- function(x,parsmc,ii,qq)
{
   th <- parsmc[1];
   ph <- parsmc[2];
   iiqq <- -sum(ii * log(1 - qq));
   y <- rep(0,length(x));
   for(cnt in 1:length(x))
   {
      if(x[cnt] > 0)
      {
         if(ph^2 *x[cnt]/(th + ph) > 100)
         {
            y[cnt] <- -th * ph * x[cnt]/(th + ph);
         } else
         {
            y[cnt] <- -ph * x[cnt] + log(expm1(ph^2 * x[cnt]/(th + ph)));
         }
         y[cnt] <- y[cnt] + log(-expm1(-iiqq * x[cnt])) - log(x[cnt]);
      } else
      {
         y[cnt] <- -Inf;
      }
   }
   y <- y + log(th);
   return(y);
}
