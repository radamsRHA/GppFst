#' read.posterior: function to read a posterior distribution of MCMC samples
#'
#' This function returns list of each step from a MCMC sampling distribution
#' @param posterior.file path to the posterior samples file (header == T)
#' @param format specify whether the file is "csv" (comma seperared) or "tab" (tab delimited)
#' @param burnin specify a proportion of samples to remove as burnin [0,1]
#' @keywords population genetics, Fst, PPS, Bayesian
#' @export
#' @examples
#' # Read MCMC sample from tab-delimited file
#' MCMC.samples <- read.posterior(posterior.file = '~/Desktop/GppFST_Tutorial/atrox_snap_gamma.log', format = "tab", burnin = .95)
#' [1] "Removing the first 381 posterior samples as burn-in"
#' [1] "Leaving 20 steps for Fst simulations/analyses"
#' 
#' MCMC.samples[1:2,1:8]  # Look at the first two MCMC steps for the first 8 columns
#' Sample posterior likelihood     prior    theta0     theta1     theta2 TreeHeightLogger
#' 	382 381000 -214844.6  -214827.1 -17.56937  0.2060823 0.05711522  0.06334359      0.005438821
#' 	383 382000 -214846.3  -214827.2 -19.05971  0.2108340 0.06029055  0.06172205      0.005967549

#####################################################################
################ Read posterior file (MCMC samples) #################
#####################################################################


read.posterior <- function(posterior.file, format, burnin){
	
	if (format == "csv"){
		posterior.samples <- read.csv(file = posterior.file, header = T) # Read as CSV
	}
	if (format == "tab"){
		posterior.samples <- read.table(file = posterior.file, header = T) # Read as tbd
	}
	
	total.samples <- length(posterior.samples[,1]) # Get total number of MCMC samples [length columnd[,1]]
	remove.burn = burnin*total.samples # Compute Burnin percentage
	burn.step <- round(x = remove.burn, digits = 0) # Compute number of steps to remove
	print(gsub(pattern = "XXX", replacement = burn.step, x = "Removing the first XXX posterior samples as burn-in")) # Print readout
	print(gsub(pattern = "XXX", replacement = total.samples-burn.step, x = "Leaving XXX steps for Fst/Dxy PPS")) # Print readout
	
	Posterior.Samples <- posterior.samples[(burn.step+1):total.samples,] # Get only post burnin samples
	
	names.vec <- vector() # Inialize vector to check names
	for (parameter in c("pop0.theta", "pop1.theta", "pop12.theta", "TreeHeightLogger")){ # For each of the 4 coalescent parameters
		names.vec <- c(names.vec, length(grep(pattern = parameter, x = names(Posterior.Samples)))) # Check if in file
	}
	
	if (sum(names.vec) < 4){
		cancel_statment = gsub(pattern = "XXX", replacement = posterior.file, x = "ERROR: CHECK parameter names in file: XXX. Should be set to c(pop0.theta, pop1.theta, pop12.theta, TreeHeightLogger)")
		stop(cancel_statment) # STOP FUNCTION if parameter names are incorrect!
	}
	
	return(Posterior.Samples) # Return only post burnin samples with new parameter names from param.legend
	
}

