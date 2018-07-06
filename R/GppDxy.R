#' GppDxy: function to simulate posterior predictive simulated distributions of dxy under the neutral coalescent model for two species
#'
#' This function returns a vector of dxy values simulated from a posterior distribution
#' @param posterior.samples List of posterior samples with columns for single column for each parameter and a row for each step of the MCMC chain [see read.posterior]
#' @param loci.per.step Number of genealogies to simulate for each joint parameter combination of the posterior MCMC chain
#' @param sample.vec.0 List of number of sampled individuals for population 0 for each locus in the empirical distribution 
#' @param sample.vec.1 List of number of sampled individuals for population 1 for each locus in the empirical distribution
#' @param sequence.length.vec List of the locus length for each locus in the empirical dataset; can be set to c(190,190) to simulate all loci with length of 190 nucleotides
#' @keywords population genetics, dxy, Bayesian, Posterior Predictive Simulation
#' @export
#' @examples
#' library(GppFst)
#' library(Geneland)
#' library(phybase)
#' 
#' experimental_params <- read.table(file = '~/Desktop/GppFST_Tutorial/ExperimentalParameters.txt', header = T) # Read tab-delimited file with experimental parameters
#' pop0.samples <- experimental_params$pop0.samples # Extract pop0 samples per empirical locus
#' pop1.samples <- experimental_params$pop1.samples # Extract pop1 samples per empirical locus
#' locus.lengths <- experimental_params$locus.length # Extract locus lengths per empirical locus
#'
#' MCMC.samples <- read.posterior(posterior.file = '~/Desktop/GppFST_Tutorial/atrox_snap_gamma2.log', format = "tab", burnin = .95)
#' 
#' Gppdxy.results <- GppDxy(posterior.samples = MCMC.samples, loci.per.step = 10, sample.vec.0 = pop0.samples, sample.vec.1 = pop1.samples, sequence.length.vec = locus.lengths)
#' mean(as.numeric(Gppdxy.results) # get mean
#' hist(as.numeric(Gppdxy.results) # plot distribution
#' 

##################################################################################
################ Genomic Posteror Predicive distributions of Dxy #################
##################################################################################

GppDxy <- function(posterior.samples, loci.per.step, sample.vec.0, sample.vec.1, sequence.length.vec){
	
	names.vec <- vector() # Inialize vector to check names
	for (parameter in c("pop0.theta", "pop1.theta", "pop12.theta", "TreeHeightLogger")){ # For each of the 4 coalescent parameters
		names.vec <- c(names.vec, length(grep(pattern = parameter, x = names(posterior.samples)))) # Check if in file
	}
	if (sum(names.vec) < 4){ # Missing parameters in posterior.samples
		print("CHECK THESE PARAMETER NAMES:")
		for (parameter in names(posterior.samples)){
			print(parameter) # Missing parameters in posterior.samples
		}
		print(gsub(pattern = "XXX", replacement = (4-sum(names.vec)), x = "The posterior.samples object is missing XXX parameters"))
		stop("ERROR: CHECK parameter names in read.posterior object. Should be set to c(pop0.theta, pop1.theta, pop12.theta, TreeHeightLogger)") # STOP FUNCTION if parameter names are incorrect!
	}
	
	
	num.MCMC.samples <- length(posterior.samples[,1]) # Get number of MCMC posterior samples
	GppDxy.vec <- character(loci.per.step*num.MCMC.samples) # Inialize list to collect dxy values from simulated loci
	locus.dxy.count = 0 # Iniatilize counter for dxy vector
	for (MCMC.step in 1:num.MCMC.samples){ # For each MCMC step in the posterior distribution
		print.step.1 <- gsub(pattern = "XXX", replacement = MCMC.step, x = "Simulating MCMC.step: XXX, Locus.per.step: YYY")
		MCMC.step.params <- posterior.samples[MCMC.step,] # Get parameter estimates for each step of the MCMC
		MCMC.step.params <- posterior.samples[MCMC.step,] # Get parameter estimates for each step of the MCMC
		theta0 <- MCMC.step.params[names(MCMC.step.params) == "pop0.theta"] # Get theta0 estimate for step
		theta1 <- MCMC.step.params[names(MCMC.step.params) == "pop1.theta"] # Get theta1 estimate for step
		theta2 <- MCMC.step.params[names(MCMC.step.params) == "pop12.theta"] # Get theta2 estimate for step
		TreeHeightLogger <- MCMC.step.params[names(MCMC.step.params) == "TreeHeightLogger"]  # Get TreeHeightLogger estimate for step
		
		pop0.step.vec <- sample(x = sample.vec.0, size = loci.per.step, replace = T) # Draw random number of diploid individuals from population 0 for each locus
		pop1.step.vec <- sample(x = sample.vec.1, size = loci.per.step, replace = T) # Draw random number of diploid individuals from population 1 for each locus 
		locus.length.vec <- sample(x = sequence.length.vec, loci.per.step, replace = T) # Draw random number sites to simulate for each locus
		for (locus in 1:loci.per.step){ # For each simulated genealogy per MCMC step
			print.step.2 <- gsub(pattern = "YYY", replacement = locus, x = print.step.1)
			print(print.step.2)
			pop0.sample.count <- pop0.step.vec[locus]*2 # Get number of population 0 haploid samples
			pop1.sample.count <- pop1.step.vec[locus]*2 # Get number of population 1 haploid samples
			total.sample.vec <- c(pop0.sample.count, pop1.sample.count) # Concatenate sample counts for population 0 and population 1
			num.sites <- locus.length.vec[locus] # Get number of sites to simulate in alignment
			
			species.model <- sptree<-"(pop0:TreeHeightLogger#theta0, pop1:TreeHeightLogger#theta1)#theta2;" # Iniatilze species model with blank parameters
			species.model <- gsub(pattern = "TreeHeightLogger", replacement = TreeHeightLogger, x = species.model) # Set TreeHeight in model [Tau]  for MCMC step estimate
			species.model <- gsub(pattern = "theta0", replacement = theta0, x = species.model) # Set theta0 for MCMC step estimate
			species.model <- gsub(pattern = "theta1", replacement = theta1, x = species.model) # Set theta1 for MCMC step estimate
			species.model <- gsub(pattern = "theta2", replacement = theta2, x = species.model) # Set theta2 for MCMC step estimate
			Simulated.Locus.Aln <- simSeqfromSp2(sptree = species.model, spname = c('pop0','pop1'), ntaxasp = total.sample.vec, ngene = 1, outfile = '', format = "nexus", seqlength = num.sites) # Simulate nucleotid alignment for step and tree
			pop0.aln <- Simulated.Locus.Aln[1:total.sample.vec[1],] # Get alignment for only population 0
			pop1.aln <- Simulated.Locus.Aln[(total.sample.vec[1]+1):sum(total.sample.vec),] # Get alignment for only population 1
			
			locus.dxy.vec <- character(pop0.sample.count*pop1.sample.count) # Init dxy vector for locus
			count.compare = 0
			pop0.hap.vec <- rep(NA, pop0.sample.count)
			pop1.hap.vec <- rep(NA, pop1.sample.count)
			
			for (pop0.sample in 1:pop0.sample.count){
				pop0.0.sample.sequence <- pop0.aln[pop0.sample,] # Get X sample
				pop0.0.sample.sequence <- paste( unlist(pop0.0.sample.sequence), collapse='')
				pop0.hap.vec[pop0.sample] <- pop0.0.sample.sequence
			}
			for (pop1.sample in 1:pop1.sample.count){
				pop1.sample.sequence <- pop1.aln[pop1.sample,] # Get X sample
				pop1.sample.sequence <- paste( unlist(pop1.sample.sequence), collapse='')
				pop1.hap.vec[pop1.sample] <- pop1.sample.sequence
			
			}
			
			for (pop0.sample in 1:pop0.sample.count){
				pop0.sample.sequence <- pop0.aln[pop0.sample,] # Get X sample
				pop0.2.sample.sequence <- paste( unlist(pop0.sample.sequence), collapse='')
				di <- table(pop0.hap.vec)[names(table(pop0.hap.vec))==pop0.2.sample.sequence]/pop0.sample.count
				
				for (pop1.sample in 1:pop1.sample.count){
					count.compare = count.compare +1 # Up count
					pop1.sample.sequence <- pop1.aln[pop1.sample,] # Get X sample
					pop1.2.sample.sequence <- paste( unlist(pop1.sample.sequence), collapse='')
					dj <- table(pop1.hap.vec)[names(table(pop1.hap.vec))==pop1.2.sample.sequence]/pop1.sample.count
					
					#x = length(pop0.sample.sequence[pop0.sample.sequence != pop1.sample.sequence]) # Compare two sequences, count number of differences
					#diff <- x/num.sites # Divide the number of different nucleotides by the length
					
					diff = (length(pop0.sample.sequence[pop0.sample.sequence != pop1.sample.sequence]))/num.sites # Compare two sequences, count number of differences
					
					locus.dxy.vec[count.compare] <- di*dj*diff # Append difference count

				}
			}
			locus.dxy <- sum(as.numeric(locus.dxy.vec))/(pop0.sample.count*pop1.sample.count) # Get locus final locus Dxy
			locus.dxy.count = locus.dxy.count + 1 # Up Count
			GppDxy.vec[locus.dxy.count] <- locus.dxy	# Append Dxy for this locus
		}
	}	
	return(GppDxy.vec) # Return Dxy counts for PPS
}
