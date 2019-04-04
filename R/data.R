#' @name datasets
#' @title Data sets of various tropical forest communities
#' @docType data
#' @description Various tree commnunity abundance data sets to test and illustrate the Independent Species approach.
#' \itemize{
#' \item dset1.abunvec contains a list of 6 samples of tree abundances from 6 tropical forest plots (BCI, Korup, Pasoh, Sinharaja, Yasuni, Lambir).
#' \item dset2.abunvec contains a list of 11 lists with one of 11 samples from BCI combined with samples from Cocoli and Sherman.
#' \item dset3.abunvec contains a list of 6 lists with 2 samples, each from one dispersal guild, for 6 tropical forest communities (BCI, Korup, Pasoh, Sinharaja, Yasuni, Lambir).
#' \item dset4a.abunvec contains a list of 6 samples from 6 censuses of BCI (1982, 1985, 1990, 1995, 200, 2005) with dbh > 1 cm.
#' \item dset4b.abunvec contains a list of 6 samples from 6 censuses of BCI (1982, 1985, 1990, 1995, 200, 2005) with dbh > 10 cm.
#' }
#' @usage data(datasets)
#' @author Rampal S. Etienne & Bart Haegeman
#' @source Condit et al. (2002). Beta-diversity in tropical forest trees. Science 295: 666-669. See also 11.	Janzen, T., B. Haegeman & R.S. Etienne (2015). A sampling formula for ecological communities with multiple dispersal syndromes. Journal of Theoretical Biology 387, 258-261.
#' @format A list of 5 data sets. See description for information on each of these data sets.
#' @keywords datasets
NULL

#' @name fitresults
#' @title Maximum likelihood estimates and corresponding likelihood values for various fits to various tropical forest communities
#' @docType data
#' @description Maximum likelihood estimates and corresponding likelihood values for various fits to various tropical forest communities, to test and illustrate the Independent Species approach.
#' \itemize{
#' \item fit1a.llikopt contains maximum likelihood values of fit of pm-dl model to dset1.abunvec
#' \item fit1a.parsopt contains maximum likelihood parameter estimates of fit of pm-dl model to dset1.abunvec
#' \item fit1b.llikopt contains maximum likelihood values of fit of pmc-dl model to dset1.abunvec
#' \item fit1b.parsopt contains maximum likelihood parameter estimates of fit of pmc-dl model to dset1.abunvec
#' \item fit2.llikopt contains maximum likelihood values of fit of rf-dl model to dset1.abunvec
#' \item fit2.parsopt contains maximum likelihood parameter estimates of fit of rf-dl model to dset1.abunvec
#' \item fit3.llikopt contains maximum likelihood values of fit of dd-dl model to dset1.abunvec
#' \item fit3.parsopt contains maximum likelihood parameter estimates of fit of dd-dl model to dset1.abunvec
#' \item fit4.llikopt contains maximum likelihood values of fit of pm-dl model to dset2.abunvec (multiple samples)
#' \item fit4.parsopt contains maximum likelihood parameter estimates of fit of pm-dl model to dset1.abunvec (multiple samples)
#' \item fit5.llikopt contains maximum likelihood values of fit of pm-dl model to dset3.abunvec (multiple guilds)
#' \item fit5.parsopt contains maximum likelihood parameter estimates of fit of pm-dl model to dset3.abunvec (multiple guilds)
#' \item fit6.llikopt contains maximum likelihood values of fit of pr-dl model to dset1.abunvec
#' \item fit6.parsopt contains maximum likelihood parameter estimates of fit of pr-dl model to dset1.abunvec
#' \item fit7.llikopt contains maximum likelihood values of fit of pm-dd model to dset1.abunvec
#' \item fit7.parsopt contains maximum likelihood parameter estimates of fit of pm-dd model to dset1.abunvec
#' \item fit8a.llikopt contains maximum likelihood values of fit of pm-dd model to dset4a.abunvec
#' \item fit8a.parsopt contains maximum likelihood parameter estimates of fit of pm-dd model to dset4a.abunvec
#' \item fit8b.llikopt contains maximum likelihood values of fit of pm-dd model to dset4b.abunvec
#' \item fit8b.parsopt contains maximum likelihood parameter estimates of fit of pm-dd model to dset4b.abunvec
#' }
#' @usage data(fitresults)
#' @author Rampal S. Etienne & Bart Haegeman
#' @source Condit et al. (2002). Beta-diversity in tropical forest trees. Science 295: 666-669.
#' @format A list of 20 lists, each containing either likelihood values or the corresponding parameter estimates. See description.
#' @keywords datasets
NULL
