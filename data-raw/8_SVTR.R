
M <- 200            #   number of rows
N <- 500            #   number of columns
Ns <- 50            #   number of Monte Carlo samples
Nl <- 50            #   number of values for lambda
lambda_max <- 50            #   maximum value for lambda
SNR <- c(0.5, 1, 2, 4)#        signal-to-noise ratio

set.seed(1)
X1          <-   matrix(rnorm(M * N, 0, 1), nrow = M, ncol = N)
X1_svd     <-   svd(X1) # USV

# rank = 100
idx <- round(0.5 * M, 0)
S <- X1_svd$d
S[(idx+1):M] <- 0
X2 <- X1_svd$u %*% diag(S) %*% t(X1_svd$v)

# rank = 10
idx         <-   round(0.05*M, 0)
S <- X1_svd$d
S[(idx+1):M] <- rep(0, M - idx)
X3 <- X1_svd$u %*% diag(S) %*% t(X1_svd$v)

# \sigma_i = \sqrt{M}/(1 + exp{(i - M/2)/20}), i = 1, ..., 200
S <- sqrt(M)/(1 + exp(((1:M) - 0.5*M)/20))
X4 <- X1_svd$u %*% diag(S) %*% t(X1_svd$v)

X1 <- X1/sqrt(sum(X1^2))
X2 <- X2/sqrt(sum(X2^2))
X3 <- X3/sqrt(sum(X3^2))
X4 <- X4/sqrt(sum(X4^2))

X0  <-  list(X1, X2, X3, X4) 
lambda <- array(NA, dim = c(4, 4, Nl))   # thresholds
tau_w <- matrix(0, 4, 4) #  noise standard-deviation
MCR <- array(0, dim = c(Nl, 4, 4)) #  Monte Carlo estimated risk
MCS <- array(0, dim = c(Nl, 4, 4, Ns)) # Monte Carlo samples
SURE <- array(0, dim = c(Nl, 4, 4)) # SURE estimate

for (Ik in 1) { # loop over matrices 1:4
  cat("Matrix", Ik, "...\n")
  X <- X0[[Ik]]
  for (In in 1) { # loop over SNR 1:4
    tau_w[Ik, In] <- 1/(SNR[In] * sqrt(M*N)) # noise standard-deviation
    lambda[Ik, In, ]  <- seq(0, lambda_max * tau_w[Ik, In], length.out = Nl) # thresholds
    
     for (Il in 1) { # loop over lambda 1:Nl
      for (Is in 1) { # Monte Carlo sample 1:Ns
        Y <- X + tau_w[Ik, In] * matrix(rnorm(M * N, 0, 1), nrow = M, ncol = N)
        Y_svd <- svd(Y)
        Sy <- Y_svd$d - lambda[Ik, In, Il]
        Sy[Sy < 0] <- 0
        SVT_Y <- Y_svd$u %*% diag(Sy) %*% t(Y_svd$v)

        MCS[Il, In, Ik, Is] <- sum((SVT_Y - X)^2)
        MCR[Il, In, Ik] <- MCR[Il, In, Ik] + MCS[Il, In, Ik, Is]
      }
      MCR[Il, In, Ik] = MCR[Il, In, Ik]/Ns
       Y <- X + tau_w[Ik, In] * matrix(rnorm(M * N, 0, 1), nrow = M, ncol = N)
       SURE[Il, In, Ik] <- sure_svt(lambda[Ik, In, Il], tau_w[Ik, In], Y)
    }
  }
}


for (In in 1) {
  for (Ik in 1:4) {
    if (Ik == 1) plot(tau_w[Ik, In] * lambda[Ik, In, ], SURE[, In, Ik])
    else points(tau_w[Ik, In] * lambda[Ik, In, ], SURE[, In, Ik])
  }
}

# synthesize some fake examples to understand the behavior

truth <- matrix(c(rep(1:10, 10)), nrow = 10)
set.seed(100)
noises <- matrix(rnorm(100, mean = 0, sd = 0.01), nrow = 10)
added_noise <- truth + noises
MSE <- sum((added_noise - truth)^2)
added_noise_svd <- rsvd::rsvd(added_noise, k = 10, q = 2)
tau <- 0.01
lambda <- seq(0, 10 * tau, length.out = 50)
SURE <- purrr::map_dbl(lambda, ~sure_svt(.x, tau, added_noise, s = added_noise_svd$d))
plot(lambda, SURE)

res <- added_noise_svd$u %*% diag(scDNAmClock::soft_threshold(added_noise_svd$d, lambda[which.min(SURE)])) %*% t(added_noise_svd$v)
sum((res - truth)^2) # 0.004

# when tau is set to the correct value, the noise reduction is optimal