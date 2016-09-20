chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
#' GppFst: function to simulate posterior predictive distributions of Fst under the neutral coalescent model for two species
#'
#' This function returns a vector of Fst values simulated from a posterior distribution
#' @param posterior.samples List of posterior samples with columns for single column for each parameter and a row for each step of the MCMC chain [see read.posterior]
#' @param SNP.per.locus Total number of fst values taken. This is arbitrary, and can be set to the locus length to get all simulated SNPs
#' @param loci.per.step Number of genealogies to simulate for each joint parameter combination of the posterior MCMC chain
#' @param sample.vec.0 List of number of sampled individuals for population 0 for each locus in the empirical distribution 
#' @param sample.vec.1 List of number of sampled individuals for population 1 for each locus in the empirical distribution
#' @param sequence.length.vec List of the number of sites for each locus in the empirical dataset (can be set to c(190,190) to simulate all loci with length of 190 nucleotides)
#' @keywords population genetics, Fst, Bayesian, PPS
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
#' Gppfst.results <- GppFst(posterior.samples = MCMC.samples, loci.per.step = 10, SNP.per.locus = 1, sample.vec.0 = pop0.samples, sample.vec.1 = pop1.samples, sequence.length.vec = locus.lengths)
#' mean(as.numeric(Gppfst.results) # get mean
#' hist(as.numeric(Gppfst.results) # plot distribution
#' 
#' Compare empirical Fst vs. simulated Fst
#' emp.fst <- experimental_params$WEIR_AND_COCKERHAM_FST
#' 
#' par(mfrow=c(2,1))
#' 
#' hist(as.numeric(Gppfst.results) # plot distribution
#' hist(as.numeric(emp.fst) # plot distribution

##################################################################################
################ Genomic Posteror Predicive distributions of Fst #################
##################################################################################

GppFst <- function(posterior.samples, loci.per.step, SNP.per.locus, sample.vec.0, sample.vec.1, sequence.length.vec){
	
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
	limit = SNP.per.locus-1
	fst.vec <- character(loci.per.step*num.MCMC.samples*SNP.per.locus) # Inialize list to collect fst values from loci
	fst.count = 0
	for (MCMC.step in 1:num.MCMC.samples){ # For each MCMC step in the posterior distribution
		print.step.1 <- gsub(pattern = "XXX", replacement = MCMC.step, x = "Simulating MCMC.step: XXX, Locus.per.step: YYY")
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
			
			poly.site.count = 0
			for (site in 1:num.sites){ # For each site in simulated locus
				if (poly.site.count <= limit){
					total.site.aln <- Simulated.Locus.Aln[,site] # Get site alignment across both populations
					
					if (length(unique(total.site.aln)) > 1){ # If site is polymorphic
						poly.site.count = poly.site.count+1
						pop0.site.aln <- pop0.aln[,site] # Get population 0 only alignment
						pop1.site.aln <- pop1.aln[,site] # Get population 1 only alignment
						
						site.genotype.aln <- matrix(nrow = sum(total.sample.vec)/2, ncol = 2) # Iniatlize a matrix with rows for each individual and columns for each haplotype
						
						#### if pop0.samples > 2 ####
						if (pop0.sample.count > 2){ # If more than 2 samples
							pop0.handle <- chunk2(x = 1:pop0.sample.count, n = (pop0.sample.count/2)) # Iniatilze list that is chunked into genotype items
							for (diploid.sample in 1:length(pop0.handle)){ # For each haplotype in individual
								site.genotype.aln[diploid.sample,] <- c(pop0.site.aln[pop0.handle[[diploid.sample]][1]], pop0.site.aln[pop0.handle[[diploid.sample]][2]]) # Append pop0 data
							}
						}
						
							#### if pop0.samples == 2 ####
						if (pop0.sample.count == 2){ # If only a single sample
							pop0.handle <- 1:2
							site.genotype.aln[1,] <- c(	pop0.site.aln[1], 	pop0.site.aln[2]) # Append pop0 site data
						}
						
						#### if pop1.samples > 2 ####
						if (pop1.sample.count > 2){ # If multiple pop1 samples
							pop1.handle <- chunk2(x = 1:pop1.sample.count, n = (pop1.sample.count/2)) # Get east sample handle
							
							for (diploid.sample in 1:length(pop1.handle)){ # for diploid sample in pop1
								pop1.samples <- (pop1.handle[diploid.sample]) # Get samples positions from pop1.handle
								diploid.aln <- c(pop1.site.aln[pop1.samples[[1]][1]], pop1.site.aln[pop1.samples[[1]][2]])  # Get diploid.aln for pop1
								site.genotype.aln[(diploid.sample+pop0.sample.count/2),] <- diploid.aln # Append pop1 site data for individual
							}
						}	
						#### if pop1.samples == 2	####
						if (pop1.sample.count	== 2){ # If only a single sample
							pop1.handle = 1:2
							site.genotype.aln[pop0.sample.count/2+1,] <- c(pop1.site.aln[1], pop1.site.aln[2]) # Append pop1 site.alignment
						}
						FST.calc <- Fstat(genotypes = site.genotype.aln, npop = 2, pop.mbrship = c(rep(1, pop0.sample.count/2), rep(2,pop1.sample.count/2))) # Get FST for data
						fst.data <- FST.calc$Fst[1,2] # FST value
						if (length(fst.data) > 0){
							fst.count = fst.count + 1
							fst.vec[fst.count] <- fst.data	
						}	
					} 
				}	else {break}
			}
		}	
	}	
	return(fst.vec)
}