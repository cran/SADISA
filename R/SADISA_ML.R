#' @title Performs maximum likelihood parameter estimation for requested model
#' @description Computes maximum loglikelihood and corresponding parameters for the requested model using the independent-species approach.
#' For optimization it uses various auxiliary functions in the DDD package.
#' @param abund abundance vector or a list of abundance vectors.
#' When a list is provided and mult = 'mg' (the default), it is assumed that the different vectors
#' apply to different guilds. When mult = 'ms' then the different vectors apply to multiple samples.
#' from the same metacommunity. In this case the vectors should have equal lengths and may contain
#' zeros because there may be species that occur in multiple samples and species that do not occur
#' in some of the samples.
#' @param initpars a vector, or - when there are multiple samples or multiple guilds - a matrix of parameter values
#' @param labelpars a vector, or - when there are multiple samples or multiple guilds - a matrix indicating whether the parameters
#' in initpars must remain fixed (0) during optimization, optimized (1), or - in the case of multiple samples or guilds - set equal
#' to the parameter of the first sample/guild (2).
#' @param model the chosen combination of metacommunity model and local community model
#' as a vector, e.g. c('pm','dl') for a model with point mutation in the metacommunity and
#' dispersal limitation.
#' The choices for the metacommunity model are: 'pm' (point mutation), 'rf' (random fission),
#' 'pr' (protracted speciation), 'dd' (density-dependence).
#' The choices for the local community model are: 'dl' (dispersal limitation), 'dd' (density-dependence).
#' @param mult When set to 'mg' (the default) the loglikelihood for multiple guilds is computed.
#' When set to 'ms' the loglikelihood for multiple samples from the same metacommunity is computed.
#' @details Not all combinations of metacommunity model and local community model have been implemented yet.
#' because this requires checking for numerical stability of the integration. The currently available model combinations are, for a single sample, c('pm','dl'), c('pm','rf'), c('dd','dl'),
#' c('pr','dl'), c('pm','dd'), and for multiple samples, c('pm','dl').
#' @param tol a vector containing three numbers for the relative tolerance in the parameters, the relative tolerance in the function, and the absolute tolerance in the parameters.
#' @param maxiter sets the maximum number of iterations
#' @param optimmethod sets the optimization method to be used, either subplex (default) or an alternative implementation of simplex.
#' @keywords model species-abundance-distribution
#' @references Haegeman, B. & Etienne, R.S. (2016). A general sampling formula for community abundance data. Methods in Ecology & Evolution. In review.
#' @examples
#' utils::data(datasets);
#' utils::data(fitresults);
#' result <- SADISA_ML(
#'    abund = datasets$dset1.abunvec[[1]],
#'    initpars = fitresults$fit1a.parsopt[[1]],
#'    labelpars = c(1,1),
#'    model = c('pm','dl'),
#'    tol = c(1E-1, 1E-1, 1E-1)
#'    );
#' # Note that tolerances should be set much lower than 1E-1 to get the best results.
#' @export
#'
SADISA_ML <- function(
   abund,
   initpars,
   labelpars,
   model = c('pm','dl'),
   mult = 'mg',
   tol = c(1E-6, 1E-6, 1E-6),
   maxiter = 1000 * round((1.25)^length(which(labelpars == 1))),
   optimmethod = 'subplex'
   )
{
   if(is.list(initpars))
   {
      ff <- NULL;
      for(i in 1:length(initpars))
      {
         ff <- rbind(ff,initpars[[i]]);
      }
      initpars <- ff;
      rm(ff);
   }
   if(is.list(labelpars))
   {
      ff <- NULL;
      for(i in 1:length(labelpars))
      {
         ff <- rbind(ff,labelpars[[i]]);
      }
      labelpars <- ff;
      rm(ff);
   }
   trpars <- initpars/(1 + initpars);
   trpars[which(initpars == Inf, arr.ind = TRUE)] <- 1;
   trparsopt <- trpars[which(labelpars == 1,arr.ind = TRUE)];
   trparsfix <- trpars[which(labelpars == 0,arr.ind = TRUE)];
   initloglik <- SADISA_loglik_choosepar(trparsopt = trparsopt,trparsfix = trparsfix,labelpars = labelpars,abund = abund,model = model,mult = mult);
   cat("The loglikelihood for the initial parameter values is ",initloglik,".\n",sep = '');
   utils::flush.console();
   if(initloglik == -Inf)
   {
      cat("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n");
      out <- NA;
      return(out);
   } else {
      optimpars <- c(tol,maxiter);
      out <- DDD::optimizer(optimmethod = optimmethod,optimpars = optimpars,fun = SADISA_loglik_choosepar,trparsopt = trparsopt,trparsfix = trparsfix,labelpars = labelpars,abund = abund,model = model,mult = mult);
   }
   if(out$conv != 0)
   {
      cat("Optimization has not converged. Try again with different initial values.\n");
      out <- NA;
      return(out);
   }
   MLtrpars <- as.numeric(unlist(out$par));
   MLpars <- MLtrpars/(1 - MLtrpars);
   ML <- as.numeric(unlist(out$fvalues));
   pars <- initpars;
   pars[which(labelpars == 1)] <- MLpars;
   ff <- which(labelpars == 2,arr.ind = TRUE);
   if(!is.null(dim(ff)))
   {
      pars[ff] <- pars[1,ff[,2]];
   }
   rm(ff);
   out <- list(pars = pars, loglik = ML, conv = unlist(out$conv));
   cat('\nParameters after likelihood maximization:\n');
   print(pars);
   cat('\nMaximum loglikelihood:\n',ML,'\n\n');
   return(out);
}

SADISA_loglik_choosepar <- function(trparsopt,trparsfix,labelpars,abund,model,mult)
{
   if(!is.list(abund))
   {
      abund <- list(abund);
   }
   trpars1 <- 0 * labelpars;
   if(is.null(dim(trpars1)))
   {
      trpars1 <- as.matrix(t(trpars1),byrow = TRUE);
   }
   trpars1[which(labelpars == 1, arr.ind = TRUE)] <- trparsopt;
   trpars1[which(labelpars == 0, arr.ind = TRUE)] <- trparsfix;
   ff <- which(labelpars == 2,arr.ind = TRUE);
   if(!is.null(dim(ff)))
   {
      trpars1[ff] <- trpars1[1,ff[,2]];
   }
   if(max(trpars1) > 1 || min(trpars1) < 0)
   {
      loglik = -Inf;
      return(loglik);
   } else {
      pars1 <- trpars1/(1 - trpars1);
      pars <- list();
      for(i in 1:length(abund))
      {
         pars[[i]] <- pars1[i,];
      }
      loglik = SADISA_loglik(abund = abund, pars = pars, model = model, mult = mult);
      if(is.nan(loglik) || is.na(loglik))
      {
         cat("There are parameter values used which cause numerical problems.\n")
         loglik <- -Inf;
         return(loglik);
      }
   }
   return(loglik);
}
