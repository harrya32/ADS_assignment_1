#######Q1

#1
set.seed(1)
TC1 = scale(rep(c(rep(1, 15), rep(0,15)), 8))

TC2 = c(rep(0, 20), rep(c(rep(1, 20), rep(0, 25)), 5))
TC2 = scale(TC2[1:240])

TC3 = scale(rep(c(rep(1, 25), rep(0, 35)), 4))

TC4 = scale(rep(c(rep(1, 15), rep(0, 25)), 6))

TC5 = scale(rep(c(rep(1, 20), rep(0, 20)), 6))

TC6 = scale(rep(c(rep(1, 25), rep(0, 15)), 6))

TC = cbind(TC1, TC2, TC3, TC4, TC5, TC6)

par(mfrow = c(2,3))
plot(TC1, type = 'l', xlab = NA)
plot(TC2, type = 'l', xlab = NA)
plot(TC3, type = 'l', xlab = NA)
plot(TC4, type = 'l', xlab = NA)
plot(TC5, type = 'l', xlab = NA)
plot(TC6, type = 'l', xlab = NA)

#2
par(mfrow = c(1,1))
par(mar=c(5,6,4,4)+.1)
plot(cor(TC), xlab = "TC", ylab = "TC", main = "TC correlation")

#3
SM1 = t(matrix(c(rep(0, 21), rep(c(0,rep(1,5),rep(0,15)), 5), rep(0,315)), 21, 21))

SM2 = t(matrix(c(rep(0, 21), rep(c(rep(0,14),rep(1,5),rep(0,2)), 5), rep(0,315)), 21, 21))

SM3 = t(matrix(c(rep(0, 21*7), rep(c(0,rep(1,5),rep(0,15)), 6), rep(0,168)), 21, 21))

SM4 = t(matrix(c(rep(0, 21*7), rep(c(rep(0,14),rep(1,5),rep(0,2)), 6), rep(0,168)), 21, 21))

SM5 = t(matrix(c(rep(0, 21*14), rep(c(0,rep(1,5),rep(0,15)), 5), rep(0,42)), 21, 21))

SM6 = t(matrix(c(rep(0, 21*14), rep(c(rep(0,14),rep(1,5),rep(0,2)), 5), rep(0,42)), 21, 21))

library('plot.matrix')
par(mfrow = c(2,3))
par(mar=c(4,4,3,4)+.1)
plot(SM1, border  = NA, xlab = NA, ylab = NA)
plot(SM2, border  = NA, xlab = NA, ylab = NA)
plot(SM3, border  = NA, xlab = NA, ylab = NA)
plot(SM4, border  = NA, xlab = NA, ylab = NA)
plot(SM5, border  = NA, xlab = NA, ylab = NA)
plot(SM6, border  = NA, xlab = NA, ylab = NA)

library(pracma)
SM1 = Reshape(SM1, 1, 441)
SM2 = Reshape(SM2, 1, 441)
SM3 = Reshape(SM3, 1, 441)
SM4 = Reshape(SM4, 1, 441)
SM5 = Reshape(SM5, 1, 441)
SM6 = Reshape(SM6, 1, 441)

SM = rbind(SM1, SM2, SM3, SM4, SM5, SM6)

par(mfrow = c(1,1))
plot(cor(t(SM)), xlab = "SM", ylab = "SM", main = "SM correlation")

#4
TC_noise = matrix(rnorm(1440, mean = 0, sd = sqrt(0.25)), 240,6)
SM_noise = matrix(rnorm(2646, mean = 0, sd = sqrt(0.015)), 6,441)

plot(cor(TC_noise), xlab = "TC", ylab = "TC", main = "TC noise correlation", key = NULL)
plot(cor(t(SM_noise)), xlab = "SM", ylab = "SM", main = "SM noise correlation")

#using sturges binning method
hist(TC_noise, breaks = "Sturges")
lines(c(-0.98, -0.98), c(0, 100), col = 2, lwd = 3)
lines(c(0.98, 0.98), c(0, 100), col = 2, lwd = 3)

hist(SM_noise, breaks = "Sturges")
lines(c(-0.24, -0.24), c(0, 200), col = 2, lwd = 3)
lines(c(0.24, 0.24), c(0, 200), col = 2, lwd = 3)

prod = TC_noise %*% SM_noise #?? what to do w this

#avg correlation between non diagonals
(sum(cor(prod)) - 441)/(441*441 - 441)
#5
X = (TC + TC_noise) %*% (SM + SM_noise)

s = X[, sample(441, 100)]
plot(s[,1], type = 'l', ylim = c(-3,3), xlab = NA, ylab = NA, 
     main = "100 randomly selected time-series")
for(i in 2:100){lines(s[,i], col = i)}

vars = rep(0, 441)
for(i in 1:441){
  vars[i] = var(X[,i])
}
plot(vars, xlab = NA, ylab = 'Variance', col = 4)

X = scale(X)


########Q2

#1
D = TC
A_lsr = abs(solve(t(D)%*%D)%*%t(D)%*%X)
D_lsr = X%*%t(A_lsr)

library('plot.matrix')
library(pracma)

par(mfrow = c(3,2))

plot(Reshape(A_lsr[1,],21,21), border  = NA, xlab = NA, ylab = NA, main = 'Retrieved SM1')
plot(D_lsr[,1], type = 'l', xlab = NA, ylab = NA, main = 'Retrieved TC1')

plot(Reshape(A_lsr[2,],21,21), border  = NA, xlab = NA, ylab = NA, main = 'Retrieved SM2')
plot(D_lsr[,2], type = 'l', xlab = NA, ylab = NA, main = 'Retrieved TC2')

plot(Reshape(A_lsr[3,],21,21), border  = NA, xlab = NA, ylab = NA, main = 'Retrieved SM3')
plot(D_lsr[,3], type = 'l', xlab = NA, ylab = NA, main = 'Retrieved TC3')

plot(Reshape(A_lsr[4,],21,21), border  = NA, xlab = NA, ylab = NA, main = 'Retrieved SM4')
plot(D_lsr[,4], type = 'l', xlab = NA, ylab = NA, main = 'Retrieved TC4')

plot(Reshape(A_lsr[5,],21,21), border  = NA, xlab = NA, ylab = NA, main = 'Retrieved SM5')
plot(D_lsr[,5], type = 'l', xlab = NA, ylab = NA, main = 'Retrieved TC5')

plot(Reshape(A_lsr[6,],21,21), border  = NA, xlab = NA, ylab = NA, main = 'Retrieved SM6')
plot(D_lsr[,6], type = 'l', xlab = NA, ylab = NA, main = 'Retrieved TC6')

par(mfrow = c(1,1))
plot(D_lsr[,3], X[,30], col = 4, xlab = "3rd retrieved TC", ylab = '30th column of X')

plot(D_lsr[,4], X[,30])

#2
lambda = 0.9
lambda_v = 6 * lambda
A_rr = abs(solve(t(D)%*%D + lambda_v * diag(6)) %*% t(D)%*%X) #take abs
D_rr = X%*%t(A_rr)

c_tlsr = diag(cor(TC, D_lsr))
c_trr = diag(cor(TC, D_rr))

sum(c_tlsr)
sum(c_trr)

#plot for lambda = 1000
lambda = 1000
lambda_v = 6 * lambda
A_rr_1000 = abs(solve(t(D)%*%D + lambda_v * diag(6)) %*% t(D)%*%X)
D_rr_1000 = X%*%t(A_rr_1000)

par(mfrow = c(1,2))
plot(Reshape(A_rr_1000[1,],21,21), border  = NA, xlab = NA, ylab = NA,
     main = "Ridge Regression with lambda = 1000")

plot(Reshape(A_lsr[1,],21,21), border  = NA, xlab = NA, ylab = NA,
     main = "Least Squares Regression")

#3
rhos = seq(0,1,0.05)
nsrcs = 6
N = 240
m <- matrix(rep(0, 21), 1,21)
all_MSE = matrix(list(m), 1,10)

for(v in 1:10){
  #generate noise
  TC_noise = matrix(rnorm(1440, mean = 0, sd = sqrt(0.25)), 240,6)
  SM_noise = matrix(rnorm(2646, mean = 0, sd = sqrt(0.015)), 6,441)
  X = (TC + TC_noise) %*% (SM + SM_noise)
  
  #find MSEs for each rho
  MSEs = matrix(0,1,21)
  for(j in 1:21){
    rho = rhos[j]
    
    #find Alr
    step <- 1/(norm(TC %*% t(TC)) * 1.1)
    thr <- rho*N*step
    Ao <- matrix(0, nsrcs, 1)
    A <- matrix(0, nsrcs, 1)
    Alr <- matrix(0, nsrcs, 441)
    
    for (k in 1:441) {
      A <- Ao+step*(t(TC) %*% (X[,k]-(TC%*%Ao)))
      A <- (1/(1+thr)) * (sign(A)*pmax(replicate(nsrcs, 0), abs(A)-thr))
      
      for (i in 1:10) {
        Ao <- A
        A <- Ao+step * (t(TC)%*%(X[,k]-(TC%*%Ao)))
        A <- (1/(1+thr)) * (sign(A)*pmax(replicate(nsrcs, 0), abs(A)-thr))
      }
      Alr[,k] <- A
    }
    
    Alr = abs(Alr)
    Dlr = X%*%t(Alr)
    MSE = sum(sum((X - Dlr%*%Alr)^2)) / (N*441)
    MSEs[j] = MSE
  }
  
  all_MSE[[1,v]] = MSEs
}

mean_MSEs = rep(0,21)
for(i in 1:21){
  MSE = 0
  for(j in 1:10){
      MSE = MSE + all_MSE[[1,j]][i]
  }
  mean_MSEs[i] = MSE / 10
}

par(mfrow = c(1,1))
plot(rhos, mean_MSEs, xlab = 'Rho', ylab = 'MSE')

plot(rhos[11:17], mean_MSEs[11:17], xlab = 'Rho', ylab = 'MSE')

min_ind = which.min(mean_MSEs)
min_rho = rhos[min_ind]


#4
step <- 1/(norm(TC %*% t(TC)) * 1.1)
thr <- min_rho*N*step
Ao <- matrix(0, nsrcs, 1)
A <- matrix(0, nsrcs, 1)
Alr <- matrix(0, nsrcs, 441)

for (k in 1:441) {
  A <- Ao+step*(t(TC) %*% (X[,k]-(TC%*%Ao)))
  A <- (1/(1+thr)) * (sign(A)*pmax(replicate(nsrcs, 0), abs(A)-thr))
  
  for (i in 1:10) {
    Ao <- A
    A <- Ao+step * (t(TC)%*%(X[,k]-(TC%*%Ao)))
    A <- (1/(1+thr)) * (sign(A)*pmax(replicate(nsrcs, 0), abs(A)-thr))
  }
  Alr[,k] <- A
}
Alr = abs(Alr)
Dlr = X%*%t(Alr)

c_srr = diag(cor(t(SM), t(A_rr)))
c_tlr = diag(cor(TC, Dlr))
c_slr = diag(cor(t(SM), t(Alr)))

sum(c_trr)
sum(c_tlr) 

sum(c_srr)
sum(c_slr)

par(mfrow = c(3,4))
for(i in 1:3){
  plot(D_rr[,i], type = 'l', xlab = NA, ylab =NA, main = paste(c('Ridge Regression Retrieved TC', i), collapse = ''))
  plot(Dlr[,i], type = 'l', xlab = NA, ylab =NA,  main = paste(c('Lasso Regression Retrieved TC', i), collapse = ''))
  plot(Reshape(A_rr[i,],21,21), border  = NA, xlab = NA, ylab =NA,  main = paste(c('Ridge Regression Retrieved SM', i), collapse = ''))
  plot(Reshape(Alr[i,],21,21), border  = NA, xlab = NA, ylab =NA,  main = paste(c('Lasso Regression Retrieved SM', i), collapse = ''))
}

for(i in 4:6){
  plot(D_rr[,i], type = 'l', xlab = NA, ylab =NA, main = paste(c('Ridge Regression Retrieved TC', i), collapse = ''))
  plot(Dlr[,i], type = 'l', xlab = NA, ylab =NA,  main = paste(c('Lasso Regression Retrieved TC', i), collapse = ''))
  plot(Reshape(A_rr[i,],21,21), border  = NA, xlab = NA, ylab =NA,  main = paste(c('Ridge Regression Retrieved SM', i), collapse = ''))
  plot(Reshape(Alr[i,],21,21), border  = NA, xlab = NA, ylab =NA,  main = paste(c('Lasso Regression Retrieved SM', i), collapse = ''))
}

#5

PCs = svd(TC)

#plot eigenvalues
par(mfrow = c(1,1))

library(ggplot2)
qplot(c(1:6), PCs$d) + 
  geom_line() +
  xlab('Principle Component') +
  ylab('Eigenvalue')

#plot regressors and source TCs
par(mfrow = c(6,2))
for(i in 1:6){
  plot(PCs$u[,i], type = 'l', xlab = NA, ylab = paste(c('Z', i), collapse = ''))
  plot(TC[,i], type = 'l', xlab = NA, ylab = paste(c('TC', i), collapse = ''))
}

Z = PCs$u
rho = 0.001
step <- 1/(norm(Z %*% t(Z)) * 1.1)
thr <- rho*N*step
Ao <- matrix(0, nsrcs, 1)
A <- matrix(0, nsrcs, 1)
A_pc <- matrix(0, nsrcs, 441)

for (k in 1:441) {
  A <- Ao+step*(t(Z) %*% (X[,k]-(Z%*%Ao)))
  A <- (1/(1+thr)) * (sign(A)*pmax(replicate(nsrcs, 0), abs(A)-thr))
  
  for (i in 1:10) {
    Ao <- A
    A <- Ao+step * (t(Z)%*%(X[,k]-(Z%*%Ao)))
    A <- (1/(1+thr)) * (sign(A)*pmax(replicate(nsrcs, 0), abs(A)-thr))
  }
  A_pc[,k] <- A
}
A_pc = abs(A_pc)
D_pc = X%*%t(A_pc)

par(mfrow = c(6,2))
for(i in 1:6){
  plot(D_pc[,i], type = 'l', xlab = NA, ylab = NA, 
       main = paste(c('Principal Component Regression Retrieved TC', i), collapse = ''))
  plot(Reshape(A_pc[i,],21,21), border  = NA, xlab = NA, ylab = NA, 
       main = paste(c('Principal Component Regression Retrieved SM', i), collapse = ''))
}

c_spc = diag(cor(t(SM), t(A_pc)))
c_tpc = diag(cor(TC, D_pc))
sum(c_spc)
sum(c_tpc)
