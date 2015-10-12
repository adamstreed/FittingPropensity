#FUNCTIONS FOR MEASURING FITTING PROPENSITY USING MODEL GUIDED SAMPLING
#Primary outputs:
#	1) Heat map of propensity to fit simulated model for region of parameter space
#		Biased walk: mean direction and variance
#		Unbiased walk: variance
#		Stasis: distribution mean and variance
#	2) Histogram of propensities for one model at a particular point in parameter space		


#Remaining Questions:
#How densely should each parameter range be sampled?
#Does it matter how many time steps we use in the simulations?
#How much variation is there in how a time series gets classified within the relevant parameter ranges?
#What value of within-population variance is appropriate?
#Does it matter that the GRW is initialized to have y-intercept 0 at time 0? How is fit calculated?

sampleUnifTimeSeries <- function(sampleSize, ns, min, max)
{
	
}

propensityHeatMap <- function(modelType, gridSize1, gridSize2, sampleSize, ns, min1, max1, min2, max2, vp, nn)
{
	heatMap <- array(rep(0,gridSize1*gridSize2),dim=c(gridSize1,gridSize2))
	param1Range <- seq(min1, max1, (max1-min1) / (gridSize1 - 1))
	param2Range <- seq(min2, max2, (max2-min2) / (gridSize2 - 1))
	
	for(i in 1:gridSize1)		
		for(j in 1:gridSize2)
		{
			result <- propensityMapPoint(modelType,sampleSize,ns,param1Range[i],param2Range[2],vp,nn)
			heatMap[i,j] <- result[modelType]
		}

	heatMap
}
#DEBUG: if modelType == 2, then param1 = 0

propensityMap <- function(modelType, gridSize1, gridSize2, sampleSize, ns, min1, max1, min2, max2, vp, nn)
{
	map <- array(rep(0,gridSize1*gridSize2),dim=c(gridSize1,gridSize2,3))
	
	if(min1 == max1)
		param1Range <- min1
	else
		param1Range <- seq(min1, max1, (max1-min1) / (gridSize1 - 1))
	if(min2 == max2)
		param2Range <- min2
	else
		param2Range <- seq(min2, max2, (max2-min2) / (gridSize2 - 1))
	
	for(i in 1:gridSize1)		
		for(j in 1:gridSize2)
			map[i,j,] <- propensityMapPoint(modelType,sampleSize,ns,param1Range[i],param2Range[2],vp,nn)

	map
}
#DEBUG: if modelType == 2, then param1 = 0; gridSize1 == 1 iff min1 == max1; gridSize2 == 1 iff min2 == max2

propensityMapPoint <- function(modelType, sampleSize, ns, param1, param2, vp, nn)
{
	if(modelType == 1 || modelType == 2) #GRW and URW
		data <- sampleGRWPropensity(sampleSize, ns, param1, param2, vp, nn)
	else if (modelType == 3)
		data <- sampleStasisPropensity(sampleSize, ns, param1, param2, vp, nn)
		
	rowMeans(data)
}
#DEBUG: if modelType == 2, then param1 = 0

#----Stasis-----

#INPUT: number of time steps (ns), distribution mean (theta), distribution variance (omega), within-population trait variance (vp)
#OUTPUT: Returns an array of size (3,1), vector of sample size of individuals for each time step (nn)
propensityStasis <- function(ns, theta, omega, vp, nn)
{
	x <- sim.Stasis(ns,theta,omega,vp,nn)
	xfit <- fit3models(x, silent = TRUE)
	xfit$modelFits$Akaike.wt	#selects akaike weights in order GRW,URW,Stasis
}
#DEBUG: ns > 0, omega > 0, vp > 0, nn = vector of length ns with entries > 0

#INPUT: number of sample fits to draw (sampleSize), number of time steps (ns), distribution mean (theta), distribution variance (omega), within-population trait variance (vp), vector of sample size of individuals for each time step (nn)
#OUTPUT: array of length (3,sampleSize) with Akaike weights of three models to simulated stasis data
sampleStasisPropensity <- function(sampleSize, ns, theta, omega, vp, nn)
{
	simData <- array(0,dim = c(3,sampleSize))
	for(i in 1:sampleSize)
		simData[,i]<- propensityStasis(ns,theta,omega,vp,nn)	#selects akaike weights in order GRW,URW,Stasis
		
	simData
}
#DEBUG: sampleSize > 0, ns > 0, omega > 0, vp > 0, nn = vector of length ns with entries > 0

#----GRW-----

#INPUT: number of time steps (ns), distribution mean (ms), distribution variance (vs), within-population trait variance (vp)
#OUTPUT: Returns an array of size (3,1), vector of sample size of individuals for each time step (nn)
propensityGRW <- function(ns, ms, vs, vp, nn)
{
	x <- sim.GRW(ns,ms,vs,vp,nn)
	xfit <- fit3models(x, silent = TRUE)
	xfit$modelFits$Akaike.wt	#selects akaike weights in order GRW,URW,Stasis
}
#DEBUG: ns > 0, ms > 0, vs > 0, vp > 0, nn = vector of length ns with entries > 0

#INPUT: number of sample fits to draw (sampleSize), number of time steps (ns), distribution mean (ms), distribution variance (vs), within-population trait variance (vp), vector of sample size of individuals for each time step (nn)
#OUTPUT: array of length (3,sampleSize) with Akaike weights of three models to simulated stasis data
sampleGRWPropensity <- function(sampleSize, ns, ms, vs, vp, nn)
{
	simData <- array(0,dim = c(3,sampleSize))
	for(i in 1:sampleSize)
		simData[,i]<- propensityGRW(ns,ms,vs,vp,nn)	#selects akaike weights in order GRW,URW,Stasis
		
	simData
}
#DEBUG: ns > 0, ms > 0, vs > 0, vp > 0, nn = vector of length ns with entries > 0


distAkaikeWts <- function(simData)
{
	
}

calcFreq <- function()
{
		partitionSize = 10 #num intervals in parameter range
	ns = 30 #num steps in simulated time series
	numReps = 10 #num time series simulated for each parameter combination
	
	msRange = seq(msMin,msMax,(msMax-msMin)/partitionSize)
	vsRange = seq(vsMin,vsMax,(vsMax-vsMin)/partitionSize)
	
	simData = array(0,dim=c(partitionSize*numReps,ns))
	
	count = 1
	for(i = 1 to partitionSize)
	{
		for(j = 1 to numReps)
		{
}

###############################################################################
#AICc distribution estimation
#input: vector of paleoTS objects

#UNIFORM SAMPLING
#input: num steps in time series, min and max of trait range
#Questions
#What value of within-population variance should be attached to the simulated data?
sampleRange <- function(ns, ymin, ymax)
{
	numReps = 10
	x <- rep(runif(ns,ymin,ymax),numReps)
}
#DEBUG: is distribution uniform over the correct range?
#DEBUG: ns > 0 & integer, ymax > ymin

############################################################################
#MODEL GUIDED SAMPLING
# input: ranges for parameters for the three models:
# Generalized Random Walk, Unbiased Random Walk, Stasis

#INPUT: min step bias, max step bias, min variance, max variance, sampling variance within time step
#OUTPUT: matrix with simulated time series, dimensions = (partitionSize*numReps, ns)
#Reduces to unbiased random walk if step bias is set to zero
sampleGRW <- function(msMin, msMax, vsMin, vsMax, vp)
{
	partitionSize = 10 #num intervals in parameter range
	ns = 30 #num steps in simulated time series
	numReps = 10 #num time series simulated for each parameter combination
	
	msRange = seq(msMin,msMax,(msMax-msMin)/partitionSize)
	vsRange = seq(vsMin,vsMax,(vsMax-vsMin)/partitionSize)
	
	simData = array(0,dim=c(partitionSize*numReps,ns))
	
	count = 1
	for(i = 1 to partitionSize)
	{
		for(j = 1 to numReps)
		{
			x <- sim.GRW(ns,msRange[i],vsRange[i],vp)
			simData[count,] <- x$mm
			count <- count + 1
		}
	}
}
#DEBUG: msMax >= msMin, vsMax >= vsMin >= 0, vp > 0

