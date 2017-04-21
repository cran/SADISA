#' @title Tests SADISA for data sets included in the paper by Haegeman & Etienne
#' @description Tests SADISA for data sets included in the paper by Haegeman & Etienne
#' @keywords model species-abundance-distribution
#' @references Haegeman, B. & Etienne, R.S. (2016). A general sampling formula for community abundance data. Methods in Ecology & Evolution. In review.
#' @export

SADISA_test <- function()
{
   datasets = NULL; rm(datasets);
   fitresults = NULL; rm(fitresults);
   utils::data('datasets', package = 'SADISA');
   utils::data('fitresults', package = 'SADISA');
   cat('\n\nTesting PM+DL model - unconditional (Table 1):\n')
   for(i in 1:6)
   {
      nn <- datasets$dset1.abunvec[[i]];
      po <- fitresults$fit1a.parsopt[[i]];
      result <- SADISA_loglik(abund = nn,pars = po,model = c('pm','dl'));
      cat('\nThe difference is:',result - fitresults$fit1a.llikopt[[i]],'  ');
      testthat::expect_equal(result,fitresults$fit1a.llikopt[[i]]);
   }
   cat('\n\nTesting PM+DL model - conditional (Table 1):\n')
   for(i in 1:6)
   {
      nn <- datasets$dset1.abunvec[[i]];
      po <- fitresults$fit1b.parsopt[[i]];
      result <- SADISA_loglik(abund = nn,pars = po,model = c('pmc','dl'));
      cat('\nThe difference is:',result - fitresults$fit1b.llikopt[[i]],'  ');
      testthat::expect_equal(result,fitresults$fit1b.llikopt[[i]]);
   }

   cat('\n\nTesting RF+DL model (Table 2):\n')
   for(i in 1:6)
   {
      nn <- datasets$dset1.abunvec[[i]];
      po <- fitresults$fit2.parsopt[[i]];
      result <- SADISA_loglik(abund = nn,pars = po,model = c('rf','dl'));
      cat('\nThe difference is:',result - fitresults$fit2.llikopt[[i]],'  ');
      testthat::expect_equal(result,fitresults$fit2.llikopt[[i]]);
   }

   cat('\n\nTesting MDD+DL model (Table 3):\n')
   for(i in 1:6)
   {
      nn <- datasets$dset1.abunvec[[i]];
      po <- fitresults$fit3.parsopt[[i]];
      result <- SADISA_loglik(abund = nn,pars = po,model = c('dd','dl'))
      cat('\nThe difference is:',result - fitresults$fit3.llikopt[[i]],'  ');
      testthat::expect_equal(result,fitresults$fit3.llikopt[[i]]);
   }

   cat('\n\nTesting multiple-samples model (Table 4):\n');
   for(i in 1:11)
   {
      nn <- datasets$dset2.abunvec[[i]];
      po <- fitresults$fit4.parsopt[[i]];
      result <- SADISA_loglik(abund = nn,pars = po,model = c('pm','dl'), mult = 'ms');
      cat('\nThe difference is:',result - fitresults$fit4.llikopt[[i]],'  ');
      testthat::expect_equal(result,fitresults$fit4.llikopt[[i]]);
   }

   cat('\n\nTesting multiple-guilds model (Table 5):\n')
   for(i in 1:6)
   {
      nn <- datasets$dset3.abunvec[[i]];
      po <- fitresults$fit5.parsopt[[i]];
      result <- SADISA_loglik(abund = nn,pars = po,model = c('pm','dl'));
      cat('\nThe difference is:',result - fitresults$fit5.llikopt[[i]],'  ');
      testthat::expect_equal(result,fitresults$fit5.llikopt[[i]]);
   }

   cat('\n\nTesting PR+DL model (Table S1):\n')
   for(i in 1:6)
   {
      nn <- datasets$dset1.abunvec[[i]];
      po <- fitresults$fit6.parsopt[[i]];
      result <- SADISA_loglik(abund = nn,pars = po,model = c('pr','dl'));
      cat('\nThe difference is:',result - fitresults$fit6.llikopt[[i]],'  ');
      testthat::expect_equal(result,fitresults$fit6.llikopt[[i]]);
   }

   cat('\n\nTesting PM+LDD model (Table S2):\n')
   for(i in 1:6)
   {
      nn <- datasets$dset1.abunvec[[i]];
      po <- fitresults$fit7.parsopt[[i]];
      result <- SADISA_loglik(abund = nn,pars = po,model = c('pm','dd'));
      cat('\nThe difference is:',result - fitresults$fit7.llikopt[[i]],'  ');
      testthat:: expect_equal(result,fitresults$fit7.llikopt[[i]]);
   }

   cat('\n\nTesting large data sets (Table S3):\n');
   # large data set
   for(i in 1:6)
   {
      nn <- datasets$dset4a.abunvec[[i]];
      po <- fitresults$fit8a.parsopt[[i]];
      result <- SADISA_loglik(abund = nn,pars = po,model = c('pm','dl'));
      cat('\nThe difference is:',result - fitresults$fit8a.llikopt[[i]],'  ');
      testthat::expect_equal(result,fitresults$fit8a.llikopt[[i]]);
   }
   # small data set
   for(i in 1:6)
   {
      nn <- datasets$dset4b.abunvec[[i]];
      po <- fitresults$fit8b.parsopt[[i]];
      result <- SADISA_loglik(abund = nn,pars = po,model = c('pm','dl'));
      cat('\nThe difference is:',result - fitresults$fit8b.llikopt[[i]],'  ');
      testthat::expect_equal(result,fitresults$fit8b.llikopt[[i]]);
   }
}
