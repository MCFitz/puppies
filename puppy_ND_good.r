##Transmission model for puppy project

require(deSolve)

puppy.model <- function (t, x, params) {
  #state variables
  S1 <- x[1]
  S2 <- x[2]
  S3 <- x[3]
  S4 <- x[4]
  E1 <- x[5]
  E2 <- x[6]
  E3 <- x[7]
  E4 <- x[8]
  I1 <- x[9]
  I2 <- x[10]
  I3 <- x[11]
  I4 <- x[12]
  V1 <- x[13]
  V2 <- x[14]
  V3 <- x[15]
  V4 <- x[16]
  Sw <- x[17]
  Ew <- x[18]
  Iw <- x[19]
  #parameters
  f3 <- params["f3"]
  f4 <- params["f4"]
  g1 <- params["g1"]
  g2 <- params["g2"]
  g3 <- params["g3"]
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
  dS1dt <- f3*(S3+E3+V3) + f4*(S4+E4+V4) - g1*S1 - mu1*S1 - k11*(I1+I2+I3+I4)*S1 - k12*Iw*S1
  dS2dt <- g1*S1 - g2*S2 - mu2*S2 - k11*(I1+I2+I3+I4)*S2 - k12*Iw*S2
  dS3dt <- g2*S2 - g3*S3 - mu2*S3 - k11*(I1+I2+I3+I4)*S3 - k12*Iw*S3
  dS4dt <- g3*S3 - mu2*S4 - k11*(I1+I2+I3+I4)*S4 - k12*Iw*S4
  dE1dt <- -g1*E1 - mu1*E1 + k11*(I1+I2+I3+I4)*S1 + k12*Iw*S1 - sigma*E1
  dE2dt <- g1*E1 - g2*E2 - mu2*E2 + k11*(I1+I2+I3+I4)*S2 + k12*Iw*S2 - sigma*E2
  dE3dt <- g2*E2 - g3*E3 - mu2*E3 + k11*(I1+I2+I3+I4)*S3 + k12*Iw*S3 - sigma*E3
  dE4dt <- g3*E3 - mu2*E4 + k11*(I1+I2+I3+I4)*S4 + k12*Iw*S4 - sigma*E4
  dI1dt <- sigma*E1 - mu1*I1 - alpha*I1
  dI2dt <- sigma*E2 - mu2*I2 - alpha*I2
  dI3dt <- sigma*E3 - mu2*I3 - alpha*I3
  dI4dt <- sigma*E4 - mu2*I4 - alpha*I4
  dV1dt <- -g1*V1 - mu1*V1
  dV2dt <- g1*V1 - g2*V2 - mu2*V2
  dV3dt <- g2*V2 - g3*V3 - mu2*V3
  dV4dt <- g3*V3 - mu2*V4
  dSwdt <- muw*(Sw+Ew+Iw) - muw*Sw + alpha*Iw - k21*(I1+I2+I3+I4)*Sw - k22*Iw*Sw
  dEwdt <- -muw*Ew + k21*(I1+I2+I3+I4)*Sw + k22*Iw*Sw - sigma*Ew
  dIwdt <- sigma*Ew - muw*Iw - alpha*Iw
  #results
  puppy.output <- c(dS1dt, dS2dt, dS3dt, dS4dt, dE1dt, dE2dt, dE3dt, dE4dt,
    dI1dt, dI2dt, dI3dt, dI4dt, dV1dt, dV2dt, dV3dt, dV4dt, dSwdt, dEwdt, 
    dIwdt)
  #list it!
  list(puppy.output)

}

times <- seq(0, 1, by = 1/48)
params <- c(
  f3 = 0.45,
  f4 = 1.72,
  g1 = 2.79,
  g2 = 3.79,
  g3 = 0.79, 
  mu1 = 2.76,
  mu2 = .45,
  muw = .45,
  sigma = 365/22.3, 
  alpha = 365/3.1,
  k11 = 1.16*365/(3.1*1.36),			#dog to dog transmission
  k12 = 0.49*365/(3.1*1.36),			#wildlife to dog transmission
  k21 = 0.13*365/(3.1*4.5),			#dog to wildlife transmission
  k22 = 0.39*365/(3.1*4.5)			  #wildlife to wildlife transmission 
)

avacc <- .3
pvacc <- .0

xstart <- c(
  S1 = 0.24*(1 - pvacc),
  S2 = 0.15*(1 - avacc),
  S3 = 0.44*(1 - avacc),
  S4 = 0.67*(1 - avacc),
  E1 = 0,
  E2 = 0,
  E3 = 0,
  E4 = 0,
  I1 = 0,
  I2 = 0,
  I3 = 0,
  I4 = 0.001,
  V1 = 0.24*pvacc,
  V2 = 0.15*avacc,
  V3 = 0.44*avacc,
  V4 = 0.67*avacc,
  Sw = 4.5,
  Ew = 0,
  Iw = 0
)

ND.out <- as.data.frame(
  ode(func = puppy.model, y = xstart, times = times, parms = params)
)

for (i in 1:9) {
  attach(ND.out)
  newstart <- c(
    S1 = S1[length(S1)],
    S2 = S2[length(S1)],
    S3 = S3[length(S1)],
    S4 = S4[length(S1)],
    E1 = E1[length(S1)],
    E2 = E2[length(S1)],
    E3 = E3[length(S1)],
    E4 = E4[length(S1)],
    I1 = I1[length(S1)],
    I2 = I2[length(S1)],
    I3 = I3[length(S1)],
    I4 = I4[length(S1)],
    V1 = V1[length(S1)],
    V2 = V2[length(S1)],
    V3 = V3[length(S1)],
    V4 = V4[length(S1)],
    Sw = Sw[length(S1)],
    Ew = Ew[length(S1)],
    Iw = Iw[length(S1)]
  )
  detach(ND.out) 
  ND.temp <- as.data.frame(
    ode(func = puppy.model, y = newstart, times = times, parms = params)
  )
  ND.temp$time <- ND.temp$time + i
  ND.out <- rbind(ND.out, ND.temp[-1,])
}  