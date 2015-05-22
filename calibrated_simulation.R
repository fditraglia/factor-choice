library("mvtnorm")

# ML Estimation of SUR model
# (In this case, equivalent to equation-by-equation OLS)
SURidentical <- function(Factors, Portfolios, inData){
  n_Obs <- nrow(inData)
  n_Portfolios <- length(Portfolios)
  n_Factors <- length(Factors)
  FactorLoadings <- matrix(NA_real_, n_Portfolios, n_Factors + 1)
  colnames(FactorLoadings) <- c("alpha", Factors)
  rownames(FactorLoadings) <- Portfolios
  FactorResiduals <- matrix(NA_real_, n_Obs, n_Portfolios)
  colnames(FactorResiduals) <- Portfolios
  
  for(i in 1:n_Portfolios){
    reg <- lm(reformulate(Factors, Portfolios[i]), inData)
    FactorLoadings[i,] <- coefficients(reg)
    FactorResiduals[,i] <- residuals(reg)
  }
  return(list(Gamma = FactorLoadings, Sigma = cov(FactorResiduals)))
}

# Calibrated simulation
SURsim <- function(Factors, Portfolios, inData, nu = Inf){
  n_Obs <- nrow(inData)
  SURresults <- SURidentical(Factors, Portfolios, inData)
  FactorLoadings <- SURresults$Gamma
  FactorData <- cbind(rep(1, n_Obs), as.matrix(inData[,Factors]))
  errorCov <- SURresults$Sigma
  if(nu == Inf){
    sim_errors <- rmvt(n_Obs, sigma = errorCov) 
  }else{
    sim_errors <- rmvt(n_Obs, sigma = errorCov, df = nu) 
  }
  sim_returns <- FactorData %*% t(FactorLoadings) + sim_errors
  return(list(x = FactorData[,-1], y = sim_returns))
}

value <- read.csv("data_value.csv")
equal <- read.csv("data_equal.csv")
FF3 <- c("Mkt.RF", "SMB", "HML")
Size10 <- c("Lo10", paste0("Dec", 2:9), "Hi10")
Industry10 <- c("NoDur", "Durbl", "Manuf", "Enrgy", 
                "HiTec", "Telcm", "Shops", "Hlth", 
                "Utils", "Other")

valueSize <- SURidentical(FF3, Size10, value)
equalSize <- SURidentical(FF3, Size10, equal)
valueIndustry <- SURidentical(FF3, Industry10, value)
equalIndustry <- SURidentical(FF3, Industry10, equal)

set.seed(8372)

# Test Run
alphaTrue <- valueSize$Gamma[,1]
testSim <- function(df){
  simValueSize <- SURsim(FF3, Size10, value, nu = df)
  simData <- with(simValueSize, cbind(x, y))
  return(SURidentical(FF3, Size10, data.frame(simData))$Gamma[,1])
}
testy <- replicate(100, testSim(df = 2))
rowMeans(testy) - alphaTrue
