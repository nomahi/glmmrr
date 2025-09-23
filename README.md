
# glmmrr package


## Logistic mixed-effects model analysis with pseudo-observations for estimating risk ratios in clustered binary data analysis

Logistic mixed-effects models are widely used as multilevel statistical models for analyzing clustered binary outcome data (e.g., longitudinal studies, cluster-randomized trials, and multi-center clinical trials). However, the resulting odds ratio can only be interpreted as an approximation of the risk ratio for rare events, and cannot be directly used as an effect measure. To address this limitation, modified Poisson regression and its extension to GEE methodology have been widely applied in recent clinical and epidemiological research. Nonetheless, these estimating-equation–based semiparametric methods cannot be readily extended to mixed-effects models. Noma and Gosho (2025, https://doi.org/10.1002/sim.70280) proposed a new approach that yields consistent risk ratio estimators within a multilevel framework by incorporating pseudo-observations into logistic mixed-effects models. A key advantage of this method is that it can be implemented with standard GLMM software simply by modifying the dataset. This package provides computational tools to perform risk ratio estimation in a multilevel modeling framework using this novel approach.



## Installation

Please download "glmmrr_2.1-1.tar.gz" and install it by R menu: "packages" -> "Install package(s) from local files...".

Download: [please click this link](https://github.com/nomahi/glmmrr/raw/main/glmmrr_2.1-1.tar.gz)

Manual: [please click this link](https://github.com/nomahi/glmmrr/raw/main/glmmrr_2.1-1.pdf)



```r

#  R example code for implementing the modified logistic mixed-effects model analysis

# The "glmmrr" package
# GitHub webpage: https://github.com/nomahi/glmmrr/

###

# Download the R package file from the following URL:
# https://github.com/nomahi/glmmrr/raw/main/glmmrr_2.1-1.tar.gz

# Then, install the package (tar.gz format) by R menu: "packages" -> "Install package(s) from local files...".

###

library("glmmrr")				# load the "glmmrr" package
library("glmmML")				# load the "glmmML" package

##

# Analyzing the example dataset "resp" using the modified logistic mixed-effects model analysis

data(resp)

gmm1 <- glmmML(y ~ treat + baseline + center + sex + age, data=resp, family =binomial, cluster=id, method="ghq")		# Ordinary logistic mixed model analysis

resp.t <- adpdt(y, resp)		# Adding pseudo-observations to the original dataset
	
gmm2 <- glmmML(y ~ treat + baseline + center + sex + age, data=resp.t, family =binomial, cluster=id, method="ghq")		# Modified logistic mixed model analysis with pseudo-observations

scoef(gmm1, eform=TRUE)			# Odds-ratio estimates
scoef(gmm2, eform=TRUE)			# Risk-ratio estimates; 95%CIs and P-values are incorrect (based on the naive model variances)

###

library("doSNOW")					# load the "doSNOW" package
library("doParallel")				# load the "doParallel" package

cl <- makeSOCKcluster(max(detectCores()-1,1))
registerDoSNOW(cl)

B <- 1000		# number of resampling

opts <- list(progress = function(x) print(paste0(x,"th bootstrap is completed.")))

R1 <- foreach(b = 1:B, .combine = rbind, .options.snow = opts) %dopar% {

	library("glmmML")
	library("glmmrr")
	
	resp.b <- cboot(id, resp)		# Creating cluster-bootstrap resamling dataset
	resp.t <- adpdt(y, resp.b) 		# Adding the pseudo-observations

	gmm.b <- glmmML(y ~ treat + baseline + center + sex + age, data=resp.t, family =binomial, cluster=cluster.b, method="ghq")		# Modified logistic mixed model analysis with pseudo-observations
	gmm.b$coefficients

}

stopCluster(cl)

scoef(gmm2, eform=TRUE)			# Risk-ratio estimates; 95%CIs and P-values are incorrect (based on the naive model variances)
sboot(R1,eform=TRUE)			# 95%CIs and P-values of risk-ratios by bootstrap




##

# Analyzing the example dataset "mch" using the modified logistic mixed-effects model analysis

data(mch)

mch2 <- mch[mch$ses==2,]		# Specify the 2nd quintile subgroup

gmm1 <- glmmML(y ~ x, data=mch2, family =binomial, cluster=SOUM, method="ghq")		# Ordinary logistic mixed model analysis

mch.t <- adpdt(y, mch2)		# Adding pseudo-observations to the original dataset
	
gmm2 <- glmmML(y ~ x, data=mch.t, family =binomial, cluster=SOUM, method="ghq")		# Modified logistic mixed model analysis with pseudo-observations

scoef(gmm1, eform=TRUE)			# Odds-ratio estimates
scoef(gmm2, eform=TRUE)			# Risk-ratio estimates; 95%CIs and P-values are incorrect (based on the naive model variances)

###

library("doSNOW")					# load the "doSNOW" package
library("doParallel")				# load the "doParallel" package

cl <- makeSOCKcluster(max(detectCores()-1,1))
registerDoSNOW(cl)

B <- 1000		# number of resampling

opts <- list(progress = function(x) print(paste0(x,"th bootstrap is completed.")))

R1 <- foreach(b = 1:B, .combine = rbind, .options.snow = opts) %dopar% {

	library("glmmML")
	library("glmmrr")
	
	mch.b <- cboot(SOUM, mch2)		# Creating cluster-bootstrap resamling dataset
	mch.t <- adpdt(y, mch.b)  		# Adding the pseudo-observations

	gmm.b <- glmmML(y ~ x, data=mch.t, family =binomial, cluster=cluster.b, method="ghq")		# Modified logistic mixed model analysis with pseudo-observations
	gmm.b$coefficients

}

stopCluster(cl)

scoef(gmm2, eform=TRUE)			# Risk-ratio estimates; 95%CIs and P-values are incorrect (based on the naive model variances)
sboot(R1,eform=TRUE)			# 95%CIs and P-values of risk-ratios by bootstrap


