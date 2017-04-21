#' @title Computes loglikelihood for requested model
#' @description Computes loglikelihood for requested model using independent-species approach
#' @param abund abundance vector or a list of abundance vectors.
#' When a list is provided and mult = 'mg' (the default), it is assumed that the different vectors
#' apply to different guilds. When mult = 'ms' then the different vectors apply to multiple samples.
#' from the same metacommunity. In this case the vectors should have equal lengths and may contain
#' zeros because there may be species that occur in multiple samples and species that do not occur
#' in some of the samples.
#' @param pars a vector of model parameters or a list of vectors of model parameters.
#' When a list is provided and mult = 'mg' (the default), it is assumed that the different vectors
#' apply to different guilds. Otherwise, it is assumed that they apply to multiple samples.
#' @param model the chosen combination of metacommunity model and local community model
#' as a vector, e.g. c('pm','dl') for a model with point mutation in the metacommunity and
#' dispersal limitation.
#' The choices for the metacommunity model are: 'pm' (point mutation), 'rf' (random fission),
#' 'pr' (protracted speciation), 'dd' (density-dependence).
#' The choices for the local community model are: 'dl' (dispersal limitation), 'dd' (density-dependence).
#' @param mult When set to 'mg' (the default) the loglikelihood for multiple guilds is computed.
#' When set to 'ms' the loglikelihood for multiple samples from the same metacommunity is computed.
#' @return loglikelihood
#' @details Not all combinations of metacommunity model and local community model have been implemented yet.
#' because this requires checking for numerical stability of the integration. The currently available model combinations are, for a single sample, c('pm','dl'), c('pm','rf'), c('dd','dl'),
#' c('pr','dl'), c('pm','dd'), and for multiple samples, c('pm','dl').
#' @keywords model species-abundance-distribution
#' @references Haegeman, B. & Etienne, R.S. (2016). A general sampling formula for community abundance data. Methods in Ecology & Evolution. In review.
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
   mult = 'mg'
   )
{
   if(!is.list(abund))
   {
      abund <- list(abund);
      pars <- list(pars);
   }
   loglik = 0;
   if(mult == 'mg')
   {
      for(i in 1:length(abund))
      {
         nn <- abund[[i]];
         nu <- sort(unique(nn));
         ss <- pracma::histc(nn,nu)$cnt;
         loglik <- loglik + model_llik(model = model, pars = pars[[i]], nn = nn, nu = nu, ss = ss);
      }
   } else  # multiple samples
   if(model[1] == 'pm' && model[2] == 'dl')
   {
      numsam <- length(abund);
      nn <- NULL;
      for(j in 1:numsam)
      {
         nn <- cbind(nn,abund[[j]]);
      }
      nu <- unique(nn);
      numuni <- dim(nu)[1];
      sf <- cbind(nu,rep(0,numuni));
      for(cnt in 1:numuni)
      {
         for(j in 1:dim(nn)[1])
         {
            sf[cnt,numsam + 1] <- sf[cnt,numsam + 1] + prod(nn[j,] == nu[cnt,]);
         }
      }

      pars2 <- pars[[1]][1];
      for(cc in 1:numsam)
      {
         pars2 <- c(pars2,pars[[cc]][2]);
      }
      loglik <- ms_llik(pars = pars2,sf = sf);
   } else
   {
      warning('Multiple samples is not implemented for this model');
      loglik <- NA;
      return(loglik);
   }
   return(loglik);
}

model_llik <- function(model,pars,nn,nu,ss)
{
   llik <- 0;
   if(model[1] == 'pm' && model[2] == 'dl')
   {
      model_estot <- pm_estot;
      model_lesk <- pm_lesk;
   } else
   if(model[1] == 'pmc' && model[2] == 'dl')
   {
      j <- sum(nn);
      qq <- j/(pars[length(pars)] + j);
      lprj <- pmc_lprj(pars = pars,qq = qq,k = j);
      llik <- llik - lprj;
      model_estot <- pm_estot;
      model_lesk <- pm_lesk;
   } else
   if(model[1] == 'rf' && model[2] == 'dl')
   {
      model_estot <- rf_estot;
      model_lesk <- rf_lesk;
   } else
   if(model[1] == 'dd' && model[2] == 'dl')
   {
      model_estot <- mdd_estot;
      model_lesk <- mdd_lesk;
   } else
   if(model[1] == 'pr' && model[2] == 'dl')
   {
      model_estot <- pr_estot;
      model_lesk <- pr_lesk;
   } else
   if(model[1] == 'pm' && model[2] == 'dd')
   {
      model_estot <- ldd_estot;
      model_lesk <- ldd_lesk;
   } else
   {
      cat('\nThe chosen combination of metacommunity model and local community model has not yet been implemented.\n');
      llik <- NA;
      return(llik);
   }
   if(model[1] == 'dd' || model[2] == 'dd')
   {
      al <- pars[2];
      if(al < -10 | al > 0.999)
      {
         llik <- -Inf;
         return(llik);
      }
   }
   logpars <- suppressWarnings(log(pars));
   logpars2 <- logpars[!is.nan(logpars) & logpars != -Inf];
   if(min(logpars2) < -20 || max(logpars2) > 20)
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
      qq = j;
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

