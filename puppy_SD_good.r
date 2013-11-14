##Transmission model for puppy project

require(deSolve)

puppy.model <- function (t, x, params) {
  #state variables
  S1 <- x[1]
  S2 <- x[2]
  S3 <- x[3]
  E1 <- x[4]
  E2 <- x[5]
  E3 <- x[6]
  I1 <- x[7]
  I2 <- x[8]
  I3 <- x[9]
  V1 <- x[10]
  V2 <- x[11]
  V3 <- x[12]
  Nd <- x[13]
  Sw <- x[14]
  Ew <- x[15]
  Iw <- x[16]
  #parameters
  f3 <- params["f3"]
  g1 <- params["g1"]
  g2 <- params["g2"]
  mu1 <- params["mu1"]
  mu2 <- params["mu2"]
  muw <- params["muw"]
  sigma <- params["sigma"]
  alpha <- params["alpha"]
  k11 <- params["k11"] 
  k12 <- params["k12"]
  k21 <- params["k21"]
  k22 <- params["k22"]
  #equations
  dS1dt <- f3*(S3+E3+V3) - g1*S1 - mu1*S1 - k11*(I1+I2+I3)*S1 - k12*Iw*S1
  dS2dt <- g1*S1 - g2*S2 - mu2*S2 - k11*(I1+I2+I3)*S2 - k12*Iw*S2
  dS3dt <- g2*S2 - mu2*S3 - k11*(I1+I2+I3)*S3 - k12*Iw*S3
  dE1dt <- -g1*E1 - mu1*E1 + k11*(I1+I2+I3)*S1 + k12*Iw*S1 - sigma*E1
  dE2dt <- g1*E1 - g2*E2 - mu2*E2 + k11*(I1+I2+I3)*S2 + k12*Iw*S2 - sigma*E2
  dE3dt <- g2*E2 - mu2*E3 + k11*(I1+I2+I3)*S3 + k12*Iw*S3 - sigma*E3
  dI1dt <- sigma*E1 - mu1*I1 - alpha*I1
  dI2dt <- sigma*E2 - mu2*I2 - alpha*I2
  dI3dt <- sigma*E3 - mu2*I3 - alpha*I3
  dV1dt <- -g1*V1 - mu1*V1
  dV2dt <- g1*V1 - g2*V2 - mu2*V2
  dV3dt <- g2*V2 - mu2*V3
  dNddt <- f3*(S3+E3+V3) - mu1*(S1+E1+I1+V1) - mu2*(S2+E2+I2+V2+S3+E3+I3+V3) - alpha*(I1+I2+I3)
  dSwdt <- muw*(Sw+Ew+Iw) - muw*Sw + alpha*Iw - k21*(I1+I2+I3)*Sw - k22*Iw*Sw
  dEwdt <- -muw*Ew + k21*(I1+I2+I3)*Sw + k22*Iw*Sw - sigma*Ew
  dIwdt <- sigma*Ew - muw*Iw - alpha*Iw
  #results
  puppy.output <- c(dS1dt, dS2dt, dS3dt, dE1dt, dE2dt, dE3dt, dI1dt, dI2dt, 
    dI3dt, dV1dt, dV2dt, dV3dt, dNddt, dSwdt, dEwdt, dIwdt)
  #list it!
  list(puppy.output)

}

times <- seq(0, 365, by = 1)
params <- c(
  f3 = 0.0049,
  g1 = 0.0067,
  g2 = 0.01, 
  mu1 = 0.01,
  mu2 = 0.0016,
  muw = 0.0016,
  sigma = 1/22.3, 
  alpha = 1/3.1,
  k11 = 1.09/(3.1*9.38),			#dog to dog transmission
  k12 = 0.09/(3.1*9.38),			#wildlife to dog transmission
  k21 = 0.95/(3.1*3),			#dog to wildlife transmission
  k22 = 0.23/(3.1*3)			  #wildlife to wildlife transmission 
)

avacc <- 0
pvacc <- 0

xstart <- c(
  S1 = 1.9*(1 - pvacc),
  S2 = 1.1*(1 - avacc),
  S3 = 6.5*(1 - avacc),
  E1 = 0,
  E2 = 0,
  E3 = 0,
  I1 = 0,
  I2 = 0,
  I3 = 0.0001,
  V1 = 1.9*pvacc,
  V2 = 1.1*avacc,
  V3 = 6.5*avacc,
  Nd = 9.5,
  Sw = 3,
  Ew = 0,
  Iw = 0
)

SD.out <- as.data.frame(
  ode(func = puppy.model, y = xstart, times = times, parms = params)
)

for (i in 1:9) {
  attach(SD.out[length(SD.out$S1),])
  newstart <- c(
    S1 = (S1+V1)*(1-pvacc),
    S2 = (S2+V2)*(1-avacc),
    S3 = (S3+V3)*(1-avacc),
    E1 = E1,
    E2 = E2,
    E3 = E3,
    I1 = I1,
    I2 = I2,
    I3 = I3,
    V1 = (S1+V1)*pvacc,
    V2 = (S2+V2)*avacc,
    V3 = (S3+V3)*avacc, #work on this - how do previously vacc pups get accounted for?
    Nd = Nd,
    Sw = Sw,
    Ew = Ew,
    Iw = Iw
  )
  detach(SD.out[length(SD.out$S1),])  
  SD.temp <- as.data.frame(
    ode(func = puppy.model, y = newstart, times = times, parms = params)
  )
  SD.temp$time <- SD.temp$time + i*365
  SD.out <- rbind(SD.out, SD.temp[-1,])
}