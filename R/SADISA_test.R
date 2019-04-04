#' @title Tests SADISA for data sets included in the paper by Haegeman & Etienne
#' @description Tests SADISA for data sets included in the paper by Haegeman & Etienne
#' @param tol tolerance of the test
#' @keywords model species-abundance-distribution
#' @references Haegeman, B. & R.S. Etienne (2017). A general sampling formula for community structure data. Methods in Ecology & Evolution. In press.
#' @export

SADISA_test <- function(tol = 1E-3)
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
      #cat('\nThe difference is:',result - fitresults$fit1a.llikopt[[i]],'  ');
      testthat::expect_equal(result,fitresults$fit1a.llikopt[[i]],tolerance = tol);
   }
   cat('\n\nTesting PM+DL model - conditional (Table 1):\n')
   for(i in 1:6)
   {
      nn <- datasets$dset1.abunvec[[i]];
      po <- fitresults$fit1b.parsopt[[i]];
      result <- SADISA_loglik(abund = nn,pars = po,model = c('pmc','dl'));
      #cat('\nThe difference is:',result - fitresults$fit1b.llikopt[[i]],'  ');
      testthat::expect_equal(result,fitresults$fit1b.llikopt[[i]],tolerance = tol);
   }

   cat('\n\nTesting RF+DL model (Table 2):\n')
   for(i in 1:6)
   {
      nn <- datasets$dset1.abunvec[[i]];
      po <- fitresults$fit2.parsopt[[i]];
      result <- SADISA_loglik(abund = nn,pars = po,model = c('rf','dl'));
      #cat('\nThe difference is:',result - fitresults$fit2.llikopt[[i]],'  ');
      testthat::expect_equal(result,fitresults$fit2.llikopt[[i]],tolerance = tol);
   }

   cat('\n\nTesting MDD+DL model (Table 3):\n')
   for(i in 1:6)
   {
      nn <- datasets$dset1.abunvec[[i]];
      po <- fitresults$fit3.parsopt[[i]];
      result <- SADISA_loglik(abund = nn,pars = po,model = c('dd','dl'))
      #cat('\nThe difference is:',result - fitresults$fit3.llikopt[[i]],'  ');
      testthat::expect_equal(result,fitresults$fit3.llikopt[[i]],tolerance = tol);
   }

   cat('\n\nTesting multiple-samples model (Table 4):\n');
   for(i in 1:11)
   {
      nn <- datasets$dset2.abunvec[[i]];
      po <- fitresults$fit4.parsopt[[i]];
      result <- SADISA_loglik(abund = nn,pars = po,model = c('pm','dl'), mult = 'ms');
      #cat('\nThe difference is:',result - fitresults$fit4.llikopt[[i]],'  ');
      testthat::expect_equal(result,fitresults$fit4.llikopt[[i]],tolerance = tol);
   }

   cat('\n\nTesting multiple-guilds model (Table 5):\n')
   for(i in 1:6)
   {
      nn <- datasets$dset3.abunvec[[i]];
      po <- fitresults$fit5.parsopt[[i]];
      result <- SADISA_loglik(abund = nn,pars = po,model = c('pm','dl'));
      #cat('\nThe difference is:',result - fitresults$fit5.llikopt[[i]],'  ');
      testthat::expect_equal(result,fitresults$fit5.llikopt[[i]],tolerance = tol);
   }

   cat('\n\nTesting PR+DL model (Table S1):\n')
   for(i in 1:6)
   {
      nn <- datasets$dset1.abunvec[[i]];
      po <- fitresults$fit6.parsopt[[i]];
      result <- SADISA_loglik(abund = nn,pars = po,model = c('pr','dl'));
      #cat('\nThe difference is:',result - fitresults$fit6.llikopt[[i]],'  ');
      testthat::expect_equal(result,fitresults$fit6.llikopt[[i]],tolerance = tol);
   }

   cat('\n\nTesting PM+LDD model (Table S2):\n')
   for(i in 1:6)
   {
      nn <- datasets$dset1.abunvec[[i]];
      po <- fitresults$fit7.parsopt[[i]];
      result <- SADISA_loglik(abund = nn,pars = po,model = c('pm','dd'));
      #cat('\nThe difference is:',result - fitresults$fit7.llikopt[[i]],'  ');
      testthat::expect_equal(result,fitresults$fit7.llikopt[[i]],tolerance = tol);
   }

   cat('\n\nTesting large data sets (Table S3):\n');
   # large data set
   for(i in 1:6)
   {
      nn <- datasets$dset4a.abunvec[[i]];
      po <- fitresults$fit8a.parsopt[[i]];
      result <- SADISA_loglik(abund = nn,pars = po,model = c('pm','dl'));
      #cat('\nThe difference is:',result - fitresults$fit8a.llikopt[[i]],'  ');
      testthat::expect_equal(result,fitresults$fit8a.llikopt[[i]],tolerance = tol);
   }
   # small data set
   for(i in 1:6)
   {
      nn <- datasets$dset4b.abunvec[[i]];
      po <- fitresults$fit8b.parsopt[[i]];
      result <- SADISA_loglik(abund = nn,pars = po,model = c('pm','dl'));
      #cat('\nThe difference is:',result - fitresults$fit8b.llikopt[[i]],'  ');
      testthat::expect_equal(result,fitresults$fit8b.llikopt[[i]],tolerance = tol);
   }

   cat('\n\nTesting multiple-samples + multiple-guilds model (new):\n')
   nn1 <- datasets$dset2.abunvec[[1]];
   nn2 <- datasets$dset2.abunvec[[2]];
   po1 <- fitresults$fit4.parsopt[[1]];
   po2 <- fitresults$fit4.parsopt[[2]];
   nn <- list();
   po <- list();
   for(i in 1:3)
   {
      nn[[i]] <- list(nn1[[i]],nn2[[i]]);
      po2[[i]][1] <- po1[[i]][1];
      po[[i]] <- list(po1[[i]],po2[[i]]);
   }
   result <- SADISA_loglik(abund = nn,pars = po,model = c('pm','dl'), mult = 'both');
   result2 <- SADISA_loglik(abund = nn1,pars = po1,model = c('pm','dl'), mult = 'ms');
   result2 <- result2 + SADISA_loglik(abund = nn2,pars = po2,model = c('pm','dl'), mult = 'ms');
   #cat('\nThe difference is:',result - result2,'  ');
   testthat::expect_equal(result,result2,tolerance = tol);

   cat('\n\nTesting multiple-samples model for protracted speciation (new):\n');
   for(i in 1:11)
   {
      nn <- datasets$dset2.abunvec[[i]];
      po <- fitresults$fit4.parsopt[[i]];
      for(j in 1:3) {po[[j]] <- c(po[[j]][1],1e8,po[[j]][-1])}
      result <- SADISA_loglik(abund = nn,pars = po,model = c('pr','dl'), mult = 'ms');
      #cat('\nThe difference is:',result - fitresults$fit4.llikopt[[i]],'  ');
      testthat::expect_equal(result,fitresults$fit4.llikopt[[i]],tolerance = tol);
   }

   cat('\n\nTesting multiple-samples model for protracted speciation (new):\n');
   ph <- c(10,100,1000)
   llik <- rbind(
      c(-2.7436e+03,-1.6509e+03,-1.2287e+03),
      c(-2.0773e+03,-1.1438e+03,-7.9608e+02),
      c(-2.0670e+03,-1.1337e+03,-7.8593e+02),
      c(-2.1189e+03,-1.1563e+03,-7.9584e+02),
      c(-2.1544e+03,-1.1734e+03,-8.0533e+02),
      c(-2.1587e+03,-1.1772e+03,-8.0737e+02),
      c(-2.1064e+03,-1.1455e+03,-7.8261e+02),
      c(-2.0874e+03,-1.1397e+03,-7.8053e+02),
      c(-2.0451e+03,-1.1162e+03,-7.6576e+02),
      c(-2.0801e+03,-1.1334e+03,-7.7544e+02),
      c(-2.1001e+03,-1.1558e+03,-8.0473e+02))
   for(i in 1:11)
   {
      for(k in 1:3)
      {
         nn <- datasets$dset2.abunvec[[i]];
         po <- list();
         for(j in 1:3) {po[[j]] <- c(100,ph[k],50)}
         result <- SADISA_loglik(abund = nn,pars = po,model = c('pr','dl'), mult = 'ms');
         testthat::expect_equal(result,llik[i,k],tolerance = tol);
      }
   }

   cat('\n\nTesting multiple-samples model for density-dependence (new):\n');
   for(i in 1:11)
   {
      nn <- datasets$dset2.abunvec[[i]];
      po <- fitresults$fit4.parsopt[[i]];
      for(j in 1:3) {po[[j]] <- c(po[[j]][1],1e-6,po[[j]][-1])}
      result <- SADISA_loglik(abund = nn,pars = po,model = c('dd','dl'), mult = 'ms');
      #cat('\nThe difference is:',result - fitresults$fit4.llikopt[[i]],'  ');
      testthat::expect_equal(result,fitresults$fit4.llikopt[[i]],tolerance = tol);
   }

   cat('\n\nTesting multiple-samples model for density-dependence (new):\n');
   al <- c(-0.5,0,0.5)
   llik <- rbind(
      c(-1.3065e+03,-1.1669e+03,-1.1384e+03),
      c(-8.6008e+02,-7.4643e+02,-7.1899e+02),
      c(-8.5003e+02,-7.3598e+02,-7.0714e+02),
      c(-8.6279e+02,-7.4400e+02,-7.1155e+02),
      c(-8.7401e+02,-7.5251e+02,-7.1812e+02),
      c(-8.7646e+02,-7.5329e+02,-7.1616e+02),
      c(-8.5023e+02,-7.2892e+02,-6.9198e+02),
      c(-8.4722e+02,-7.2693e+02,-6.9118e+02),
      c(-8.3058e+02,-7.1367e+02,-6.7994e+02),
      c(-8.4190e+02,-7.2173e+02,-6.8463e+02),
      c(-8.6966e+02,-7.5544e+02,-7.2795e+02))
   for(i in 1:11)
   {
      for(k in 1:3)
      {
         nn <- datasets$dset2.abunvec[[i]];
         po <- list();
         for(j in 1:3) {po[[j]] <- c(100,al[k],50)}
         result <- SADISA_loglik(abund = nn,pars = po,model = c('dd','dl'), mult = 'ms');
         testthat::expect_equal(result,llik[i,k],tolerance = tol);
      }
   }
}
