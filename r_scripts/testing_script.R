#############################################################
## nlammers@berkeley.edu - April 2020
#############################################################

##############################
# 0 - Load librairies
##############################
library(deGradInfer)
library(deSolve)
library(data.table)

##############################
# 1 - Test example code chunk from manual
##############################
dataTest <- LV_example_dataset$data
timeTest <- LV_example_dataset$time
noiseTest <- LV_example_dataset$noise

LV_func = function(t, X, params) {
  dxdt = cbind(
    X[,1]*(params[1] - params[2]*X[,2]),
    - X[,2]*(params[3] - params[4]*X[,1])
  )
  return(dxdt)
}
# Example run only; to achieve convergence the number of iterations and
# chains must be increased.
param.result = agm(data=dataTest,time=timeTest,noise.sd=0.31,ode.system=LV_func,
                   numberOfParameters=4,temperMismatchParameter=TRUE,
                   chainNum=4, maxIterations=150,originalSignalOnlyPositive=TRUE,
                   logPrior="Gamma",defaultTemperingScheme="LB10")
print(param.result$posterior.mean)


##############################
# 2 - Appears to work. Next let's try to just solve an ODE we define numerically
##############################
parameters <- c(k1 = 1.2,
                k2 = 2.5)

state <- c(A = 1,
           S = 200)

Auto.Cat<-function(t, state, parameters) {
    with(as.list(c(state, parameters)),{
      # rate of change
        dA <- A*(10-A)*k1
        dS <- -A*S*k2
           # return the rate of change
           list(c(dA, dS))
         }) # end with(as.list ...
      }

times <- seq(0, 1, by = 0.01)
out <- ode(y = state, times = times, func = Auto.Cat, parms = parameters)
tail(out)
outDf <- as.data.frame(out)
g<-ggplot(outDf, aes(x=time, y=200-S)) + geom_line()+ geom_point()
#g1 <- g + coord_cartesian(xlim=c(0,0.1), ylim=c(0, 1000000))  # zooms in

# Add Title and Labels
#g + labs(title="Simple NCR Reaction",  y="cleaved reporter", x="time (seconds)")

# or

g + ggtitle("Simple NCR Reaction", subtitle="cleaved reporter") + xlab("time") + ylab("fluorescence")
plot(g)


##############################
# 3 - Shockingly straight forward. 
#     Now, let's see if I can infer back the parameters from my simulated system
##############################

ncrTest <- data.matrix(outDf[,2:3])

Auto.Cat.Fun = function(t, X, params) {
  dxdt = cbind(
    X[,1]*(10-X[,1])*params[1],
    - X[,1]*X[,2]*params[2]
  )
  return(dxdt)
}
# Example run only; to achieve convergence the number of iterations and
# chains must be increased.
param.result2 = agm(data=ncrTest,time=times,  noise.sd=0.31,  ode.system=Auto.Cat.Fun,
                   numberOfParameters=2,  temperMismatchParameter=TRUE,
                   chainNum=4, maxIterations=150,  originalSignalOnlyPositive=TRUE,
                   logPrior="Gamma",defaultTemperingScheme="LB10")
print(param.result2$posterior.mean)
