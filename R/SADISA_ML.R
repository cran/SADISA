#' @title Performs maximum likelihood parameter estimation for requested model
#' @description Computes maximum loglikelihood and corresponding parameters for the requested model using the independent-species approach.
#' For optimization it uses various auxiliary functions in the DDD package.
#' @param abund abundance vector or a list of abundance vectors.
#' When a list is provided and mult = 'mg' (the default), it is assumed that the different vectors
#' apply to different guilds. When mult = 'ms' then the different vectors apply to multiple samples.
#' from the same metacommunity. In this case the vectors should have equal lengths and may contain
#' zeros because there may be species that occur in multiple samples and species that do not occur
#' in some of the samples.
#' @param initpars a vector of initial values of the parameters to be optimized and fixed. See \code{labelpars}
#' for more explanation.
#' @param idpars a vector stating whether the parameters in \code{initpars} should be optimized (1) or
#' remain fixed (0).
#' @param labelpars a vector, a list of vectors or a list of lists of vectors indicating the labels
#' integers (starting at 1) of the parameters to be optimized and fixed. These integers correspond to
#' the position in \code{initpars} and \code{idpars}. The order of the labels in
#' the vector/list is first the metacommunity parameters (theta, and phi (for protracted speciation) or alpha
#' (for density-dependence or abundance-dependent speciation)), then the dispersal parameters (I).
#' See the example and the vignette for more explanation.
#' @param model the chosen combination of metacommunity model and local community model
#' as a vector, e.g. c('pm','dl') for a model with point mutation in the metacommunity and
#' dispersal limitation.
#' The choices for the metacommunity model are: 'pm' (point mutation), 'rf' (random fission),
#' 'pr' (protracted speciation), 'dd' (density-dependence).
#' The choices for the local community model are: 'dl' (dispersal limitation), 'dd' (density-dependence).
#' @param mult When set to 'single' (the default), the loglikelihood for a single sample and single guild
#' is computed. When set to 'mg', the loglikelihood for multiple guilds is computed.
#' When set to 'ms' the loglikelihood for multiple samples from the same metacommunity is computed.
#' @details Not all combinations of metacommunity model and local community model have been implemented yet.
#' because this requires checking for numerical stability of the integration. The currently available model combinations are, for a single sample, c('pm','dl'), c('pm','rf'), c('dd','dl'),
#' c('pr','dl'), c('pm','dd'), and for multiple samples, c('pm','dl').
#' @param tol a vector containing three numbers for the relative tolerance in the parameters, the relative tolerance in the function, and the absolute tolerance in the parameters.
#' @param maxiter sets the maximum number of iterations
#' @param optimmethod sets the optimization method to be used, either subplex (default) or an alternative implementation of simplex.
#' @param num_cycles the number of cycles of opimization. If set at Inf, it will do as many cycles as needed to meet the tolerance set for the target function.
#' @keywords model species-abundance-distribution
#' @references Haegeman, B. & R.S. Etienne (2017). A general sampling formula for community structure data. Methods in Ecology & Evolution 8: 1506-1519. doi: 10.1111/2041-210X.12807
#' @examples
#' utils::data(datasets);
#' utils::data(fitresults);
#' result <- SADISA_ML(
#'    abund = datasets$dset1.abunvec[[1]],
#'    initpars = fitresults$fit1a.parsopt[[1]],
#'    idpars = c(1,1),
#'    labelpars = c(1,2),
#'    model = c('pm','dl'),
#'    tol = c(1E-1, 1E-1, 1E-1)
#'    );
#' # Note that tolerances should be set much lower than 1E-1 to get the best results.
#' @export
#'
SADISA_ML <- function(
   abund,
   initpars,
   idpars,
   labelpars,
   model = c('pm','dl'),
   mult = 'single',
   tol = c(1E-6, 1E-6, 1E-6),
   maxiter = min(1000 * round((1.25)^sum(idpars)),100000),
   optimmethod = 'subplex',
   num_cycles = 1
   )
{
   if(!(is.nonnegativewholenumber(maxiter) && maxiter > 0))
   {
      stop('The number of iterations should be positive.')
   }
   if(mult != 'single' && !is.list(labelpars))
   {
      stop('The labels of the parameters should be in a list when there are multiple samples or guilds.')
   }
   dim_abund <- NULL
   dim_labelpars <- NULL
   abundff <- abund
   labelparsff <- labelpars
   for(i in 1:4)
   {
      dim_abund[i] <- length(abundff)
      dim_labelpars[i] <- length(labelparsff)
      abundff <- abundff[[1]]
      labelparsff <- labelparsff[[1]]
   }
   rm(abundff,labelparsff)
   #print(dim_abund)
   #print(dim_labelpars)
   depth_abund <- min(which(dim_abund == 1))
   depth_labelpars <- min(which(dim_labelpars == 1))
   if(depth_abund != depth_labelpars)
   {
      stop('The abundance data and the parameter labels do not match.')
   }
   initpars <- DDD::transform_pars(initpars);
   if(!is.null(which(idpars == 1)))
   {
      trparsopt <- initpars[which(idpars == 1)];
   } else
   {
      trparsopt <- NULL;
      stop('You are not optimizing anything.')
   }
   if(!is.null(which(idpars == 0)))
   {
      trparsfix <- initpars[which(idpars == 0)];
   } else
   {
      trparsfix <- NULL;
   }
   initloglik <- SADISA_loglik_choosepar(trparsopt = trparsopt,trparsfix = trparsfix,idpars = idpars,labelpars = labelpars,abund = abund,model = model,mult = mult);
   cat("The loglikelihood for the initial parameter values is ",initloglik,".\n",sep = '');
   utils::flush.console();
   if(initloglik == -Inf)
   {
      cat("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n");
      out <- NA;
      return(out);
   } else {
      optimpars <- c(tol,maxiter);
      out <- DDD::optimizer(optimmethod = optimmethod,optimpars = optimpars,num_cycles = num_cycles,fun = SADISA_loglik_choosepar,trparsopt = trparsopt,trparsfix = trparsfix,idpars = idpars,labelpars = labelpars,abund = abund,model = model,mult = mult);
   }
   if(out$conv != 0)
   {
      cat("Optimization has not converged. Try again with different initial values.\n");
      out <- NA;
      return(out);
   }
   MLtrpars <- as.numeric(unlist(out$par));
   MLpars <- DDD::untransform_pars(MLtrpars);
   parsfix <- DDD::untransform_pars(trparsfix);
   ML <- as.numeric(unlist(out$fvalues));
   pars <- make_list_of_pars(MLpars,parsfix,idpars,labelpars)
   out <- list(pars = pars, loglik = ML, conv = unlist(out$conv));
   cat('\nParameters after likelihood maximization:\n');
   print(pars);
   cat('\nMaximum loglikelihood:\n',ML,'\n\n');
   return(out);
}

SADISA_loglik_choosepar <- function(trparsopt,trparsfix,idpars,labelpars,abund,model,mult)
{
   if(!is.list(abund))
   {
      abund <- list(abund);
   }
   trpars1 <- make_list_of_pars(trparsopt,trparsfix,idpars,labelpars)
   unlisttrpars1 <- unlist(trpars1);
   if(max(unlisttrpars1) > 1 | min(unlisttrpars1) < -(model[1] == 'dd'))
   {
      loglik = -Inf;
      return(loglik);
   } else {
      pars <- untransform_list_of_pars(trpars1)
      loglik <- SADISA_loglik(abund = abund, pars = pars, model = model, mult = mult);
      if(is.nan(loglik) | is.na(loglik))
      {
         cat("There are parameter values which cause numerical problems.\n")
         loglik <- -Inf;
         return(loglik);
      }
   }
   return(loglik);
}

make_list_of_pars <- function(trparsopt,trparsfix,idpars,labelpars)
{
   trpars1 <- labelpars;
   j <- 1;
   k <- 1;
   if(!is.list(trpars1))
   {
      for(i in 1:length(idpars))
      {
         lb <- which(labelpars == i)
         if(idpars[i] == 0)
         {
            trpars1[lb] <- trparsfix[j];
            j <- j + 1;
         } else
         {
            trpars1[lb] <- trparsopt[k];
            k <- k + 1;
         }
      }
   } else if(!is.list(trpars1[[1]]))
   {
      for(i in 1:length(idpars))
      {
         if(idpars[i] == 0)
         {
            for(sg in 1:length(trpars1))
            {
               lb <- which(labelpars[[sg]] == i)
               trpars1[[sg]][lb] <- trparsfix[j];
            }
            j <- j + 1;
         } else
         {
            for(sg in 1:length(trpars1))
            {
               lb <- which(labelpars[[sg]] == i)
               trpars1[[sg]][lb] <- trparsopt[k];
            }
            k <- k + 1;
         }
      }
   } else
   {
      dim1 <- length(trpars1);
      dim2 <- length(trpars1[[1]]);
      for(i in 1:length(idpars))
      {
         if(idpars[i] == 0)
         {
            for(s in 1:dim1)
            {
               for(g in 1:dim2)
               {
                  lb <- which(labelpars[[s]][[g]] == i)
                  trpars1[[s]][[g]][lb] <- trparsfix[j];
               }
            }
            j <- j + 1;
         } else
         {
            for(s in 1:dim1)
            {
               for(g in 1:dim2)
               {
                  lb <- which(labelpars[[s]][[g]] == i)
                  trpars1[[s]][[g]][lb] <- trparsopt[k];
               }
            }
            k <- k + 1;
         }
      }
   }
   return(trpars1)
}

untransform_list_of_pars <- function(trpars1)
{
   pars <- trpars1;
   if(!is.list(trpars1))
   {
      pars <- DDD::untransform_pars(trpars1);
   } else if(!is.list(trpars1[[1]]))
   {
      for(i in 1:length(trpars1))
      {
         pars[[i]] <- DDD::untransform_pars(trpars1[[i]]);
      }
   } else
   {
      dim1 <- length(trpars1);
      dim2 <- length(trpars1[[1]]);
      for(i in 1:dim1)
      {
         for(j in 1:dim2)
         {
            pars[[i]][[j]] <- DDD::untransform_pars(trpars1[[i]][[j]]);
         }
      }
   }
   return(pars)
}
