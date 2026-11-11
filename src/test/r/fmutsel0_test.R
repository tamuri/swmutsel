# EXAMPLE MODEL ----------------------------------------------------------------

# The eigendecomposition of the Q matrix, where Q is
# the FMutSel0 model of Yang & Nielsen (2008), with
# kappa = 1, omega = 1, pi = {0.25, 0.25, 0.25, 0.25} and all fitness = { 0 }

# Eigendecomposition of Q where Q = U L V
# Columns of U are the right eigenvectors of rate matrix Q
U <- matrix(scan("U.dat"), byrow = T, nrow = 64, ncol = 64)
# Eigenvalues
L <- scan("L.dat")
# Rows of V = U^(-1) are the left eigenvectors of Q
V <- matrix(scan("V.dat"), byrow = T, nrow = 64, ncol = 64)

# The stationary frequencies of the model
Freqs <- scan("Freqs.dat")

# EXAMPLE DATA -----------------------------------------------------------------

# The conditional likelihood at tip p and length of p to its parent (r)
p <- scan("pCond.dat")
p_len <- scan("pLen.dat")

# The conditional likelihood at tip q and length of q to its parent (r)
q <- scan("qCond.dat")
q_len <- scan("qLen.dat")

# CONDITIONAL LIKELIHOOD AT PARENT ---------------------------------------------

calc_cond_like <- function(p, p_len, q, q_len, U, L, V) {
  r <- (U %*% diag(exp(L * p_len)) %*% V %*% p) * # Tip p
       (U %*% diag(exp(L * q_len)) %*% V %*% q)   # Tip q
  r[,1]
}

get_log_like <- function(freqs, cond) {
  log(sum(freqs * cond))
}

# CHECK RESULT -----------------------------------------------------------------
results <- scan("result.dat")
r <- calc_cond_like(p, p_len, q, q_len, U, L, V)
tol = 1e-12
sum((r - results) <= tol) == 64

# likelihood
get_log_like(Freqs, r)

# A STRESS TEST ----------------------------------------------------------------
branches <- 50 # make this a large number
props <- runif(branches * 2, min = 0, max = 1)
len_props <- matrix(props, ncol = 2)
conditionals <- apply(len_props, MARGIN = 1, FUN = function(x) {
  calc_cond_like(p = p, p_len = x[1], q = q, q_len = x[2], U = U, L = L, V = V)
})

log_likes <- apply(conditionals, MARGIN = 2, FUN = function(x) { get_log_like(Freqs, x) })
sum(len_props[which.max(log_likes),])
max(log_likes)

# MAXIMUM LIKELIHOOD EXAMPLE ---------------------------------------------------
opt <- optim(par = c(p_len, q_len), 
      control = list(fnscale = -1), 
      method="L-BFGS-B", 
      fn = function(x) {
        get_log_like(Freqs, calc_cond_like(p, x[1], q, x[2], U, L, V))
      }, 
      lower = c(0,0), 
      upper = c(1,1))
sum(opt$par)
opt$value

