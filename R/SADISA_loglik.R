#' @title Computes loglikelihood for requested model
#' @description Computes loglikelihood for requested model using independent-species approach
#' @param abund abundance vector or a list of abundance vectors.
#' When a list is provided and mult = 'mg' (the default), it is assumed that the different vectors
#' apply to different guilds. When mult = 'ms' then the different vectors apply to multiple samples
#' from the same metacommunity. In this case the vectors should have equal lengths and may contain
#' zeros because there may be species that occur in multiple samples and species that do not occur
#' in some of the samples. When mult= 'both', abund should be a list of lists, each list representing
#' multiple guilds within a sample
#' @param pars a vector of model parameters or a list of vectors of model parameters.
#' When a list is provided and mult = 'mg' (the default), it is assumed that the different vectors
#' apply to different guilds. Otherwise, it is assumed that they apply to multiple samples.
#' @param model the chosen combination of metacommunity model and local community model
#' as a vector, e.g. c('pm','dl') for a model with point mutation in the metacommunity and
#' dispersal limitation.
#' The choices for the metacommunity model are: 'pm' (point mutation), 'rf' (random fission),
#' 'pr' (protracted speciation), 'dd' (density-dependence).
#' The choices for the local community model are: 'dl' (dispersal limitation), 'dd' (density-dependence).
#' @param mult When set to 'single' (the default), the loglikelihood for a single sample is computed
#' When set to 'mg' the loglikelihood for multiple guilds is computed.
#' When set to 'ms' the loglikelihood for multiple samples from the same metacommunity is computed.
#' When set to 'both' the loglikelihood for multiple guilds within multiple samples is computed.
#' @return loglikelihood
#' @details Not all combinations of metacommunity model and local community model have been implemented yet.
#' because this requires checking for numerical stability of the integration. The currently available model combinations are, for a single sample, c('pm','dl'), c('pm','rf'), c('dd','dl'),
#' c('pr','dl'), c('pm','dd'), and for multiple samples, c('pm','dl').
#' @keywords model species-abundance-distribution
#' @references Haegeman, B. & R.S. Etienne (2017). A general sampling formula for community structure data. Methods in Ecology & Evolution 8: 1506-1519. doi: 10.1111/2041-210X.12807
#' @examples
#' data(datasets);
#' abund_bci <- datasets$dset1.abunvec[[1]];
#' data(fitresults);
#' data.paropt <- fitresults$fit1a.parsopt[[1]];
#' result <- SADISA_loglik(abund = abund_bci,pars = data.paropt,model = c('pm','dl'));
#' cat('The difference between result and the value in fitresults.RData is:',
#' result - fitresults$fit1a.llikopt[[1]]);
#' @export

SADISA_loglik <- function(
   abund,
   pars,
   model,
   mult = 'single'
   )
{
   if(!is.list(abund))
   {
      abund <- list(abund);
   }
   if(!is.list(pars))
   {
      pars <- list(pars);
   }
   if((length(abund) != 1) && (length(abund) != length(pars)))
   {
      stop('The list of abundances has a different dimension than the list of parameters.');
   }
   loglik = 0;
   if(mult == 'mg' | mult == 'single')
   {
      for(i in 1:length(abund))
      {
         nn <- abund[[i]];
         if(min(nn) == 0)
         {
            nn <- nn[-which(nn == 0)];
         }
         nu <- sort(unique(nn));
         if(!is.nonnegativewholenumber(nu))
         {
            stop('The abundances should be non-negative integers.')
         }
         ss <- pracma::histc(nn,nu)$cnt;
         loglik <- loglik + model_llik(model = model, pars = pars[[i]], nn = nn, nu = nu, ss = ss);
      }
   } else # multiple samples
   {
      abund2 <- abund
      numsam <- length(abund2)
      pars2 <- pars
      if(mult == 'both' && is.list(abund2[[1]]))
      {
         numguilds <- length(abund2[[1]])
      } else
      {
         numguilds <- 1
      }
      if(!((model[1] == 'pm' | model[1] == 'pr' | model[1] == 'dd') & model[2] == 'dl'))
      {
         stop('Multiple samples is not implemented for this model');
      } else
      {
         for(g in 1:numguilds)
         {
            if(numguilds == 1)
            {
               abund <- abund2
               pars <- pars2
               sumabund <- rep(0,length(abund2[[1]]))
               for(j in 1:numsam)
               {
                  sumabund <- sumabund + abund2[[j]]
               }
               if(min(sumabund) == 0)
               {
                  stop('There are species that have zero abundance across all samples.');
               }
            } else
            {
               abund <- list()
               pars <- list()
               sumabund <- rep(0,length(abund2[[1]][[g]]))
               for(j in 1:numsam)
               {
                  abund[[j]] <- abund2[[j]][[g]]
                  pars[[j]] <- pars2[[j]][[g]]
                  sumabund <- sumabund + abund[[j]]
               }
               if(min(sumabund) == 0)
               {
                  stop('There are species that have zero abundance across all samples.');
               }
            }
            nn <- NULL;
            for(j in 1:numsam)
            {
               if(length(abund[[j]]) != length(abund[[1]]) )
               {
                  stop(paste0('Sample ',j, ' for guild ', g, ' does not have the same length as sample 1.\n'))
               }
               nn <- cbind(nn,abund[[j]]);
            }
            nu <- unique(nn);
            if(!is.nonnegativewholenumber(nu))
            {
               stop('The abundances should be non-negative integers.')
            }
            numuni <- dim(nu)[1];
            sf <- cbind(nu,rep(0,numuni));
            for(cnt in 1:numuni)
            {
               for(j in 1:dim(nn)[1])
               {
                  sf[cnt,numsam + 1] <- sf[cnt,numsam + 1] + prod(nn[j,] == nu[cnt,]);
               }
            }
            if(model[1] == 'pm' | model[1] == 'rf')
            {
               parsnew <- pars[[1]][1];
               for(cc in 1:numsam)
               {
                  parsnew <- c(parsnew,pars[[cc]][2]);
               }
            } else
            if(model[1] == 'pr' | model[1] == 'dd')
            {
               parsnew <- pars[[1]][1:2];
               for(cc in 1:numsam)
               {
                  parsnew <- c(parsnew,pars[[cc]][3]);
               }
            }
            #loglik <- loglik + pmdlms_llik(pars = parsnew,sf = sf
            loglik <- loglik + ms_llik(pars = parsnew,sf = sf,model = model);
         }
      }
   }
   if(loglik > 0)
   {
      stop('The loglikelihood is larger than 0. This is probably due to a numerical problem.')
   }
   return(loglik);
}

model_llik <- function(model,pars,nn,nu,ss)
{
   llik <- 0;
   if(model[1] == 'pm' & model[2] == 'dl')
   {
      model_estot <- pm_estot;
      model_lesk <- pm_lesk;
   } else
   if(model[1] == 'pmc' & model[2] == 'dl')
   {
      j <- sum(nn);
      qq <- j/(pars[length(pars)] + j);
      lprj <- pmc_lprj(pars = pars,qq = qq,k = j);
      llik <- llik - lprj;
      model_estot <- pm_estot;
      model_lesk <- pm_lesk;
   } else
   if(model[1] == 'rf' & model[2] == 'dl')
   {
      model_estot <- rf_estot;
      model_lesk <- rf_lesk;
   } else
   if(model[1] == 'dd' & model[2] == 'dl')
   {
      model_estot <- mdd_estot;
      model_lesk <- mdd_lesk;
   } else
   if(model[1] == 'pr' & model[2] == 'dl')
   {
      model_estot <- pr_estot;
      model_lesk <- pr_lesk;
   } else
   if(model[1] == 'pm' & model[2] == 'dd')
   {
      model_estot <- ldd_estot;
      model_lesk <- ldd_lesk;
   } else
   {
      stop('The chosen combination of metacommunity model and local community model has not yet been implemented.');
      llik <- NA;
      return(llik);
   }
   if(model[1] == 'dd' | model[2] == 'dd')
   {
      al <- pars[2];
      if(al < -10 | al > 0.99)
      {
         warning('alpha is larger than 1 or very negative.')
         llik <- -Inf;
         return(llik);
      }
   }
   logpars <- log(pars[pars > 0]);
   if(min(logpars) < -20 | max(logpars) > 20)
   {
      llik <- -Inf;
      return(llik);
   }
   j <- sum(nn);
   qq <- j/(pars[length(pars)] + j);
   if(model[2] == 'dd')
   {
      funzero <- function(x) ldd_ejtot(pars = pars,qq = x) - j;
      reszero <- stats::uniroot(f = funzero,interval = c(1e-5,1-1e-5), tol = 1e-10);
      qq <- reszero$root;
   } else
   if(model[2] == 'ss')
   {
      qq <- j;
   }
   estot <- model_estot(pars,qq);
   lesk <- rep(0,length(nu));
   for(cnt in 1:length(nu))
   {
      lesk[cnt] <- model_lesk(pars,qq,k = nu[cnt]);
   }
   llik <- llik - estot + sum(ss * lesk) - sum(lgamma(ss + 1));
   return(llik);
}
