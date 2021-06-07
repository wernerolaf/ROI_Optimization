library(ROI)
library(slam)
library(ROI.plugin.qpoases)
library(lars)

data('diabetes')
x <- diabetes$x # 10 variables, standarized
y <- diabetes$y # output

############################ LASSO #############################

#----------- with glmnet -------------#

library(glmnet)

lambda <- 1
mo <- glmnet(x, y, family = "gaussian", alpha = 1, lambda = lambda, 
             intercept = FALSE, standardize = FALSE)
glmnet_beta <- setNames(as.vector(coef(mo)), rownames(coef(mo)))
round(glmnet_beta, 4) # I guess the result is equivalent to the 3rd row in object$beta but the cofficients are different

#------------ with ROI -----------------#


# I guess we need to make solver using:

# ROI_plugin_register_solver_method(signatures, solver, method, plugin = solver)


# Proba z glmnet #

glmnet_solve_OP <- function(x, control = list()) {
  
  control$canonicalize_status <- FALSE
  mo <- glmnet(x = constraints(x)[1]$L[1:sum(constraints(x)[[3]]!=0),1:(sum(constraints(x)[[3]]==0)/2)], y=constraints(op)[[3]][1:sum(constraints(op)[[3]]!=0)], family = "gaussian", alpha = 1, lambda = x$objective[2]$L[x$objective[2]$L$ncol], 
               intercept = FALSE, standardize = FALSE)
  glmnet_beta <- setNames(as.vector(coef(mo)), rownames(coef(mo)))
  ROI_plugin_canonicalize_solution(solution = glmnet_beta,
                                   optimum = mo$dev.ratio, 
                                   status = 0L, 
                                   solver = "glmnet3")
}

glmnet_signature <- ROI_plugin_make_signature(objective = "Q",
                                            constraints = "L", 
                                            types = c("C", "I", "B", "CI", "CB", "IB", "CIB"),
                                            bounds = c("X", "V"), maximum = c(TRUE, FALSE))

ROI_plugin_register_solver_method(glmnet_signature, "glmnet3", glmnet_solve_OP)


solver <- "glmnet3"
ROI_plugin_add_status_code_to_db(solver, 0L, "OK", "Solution is optimal.", 0L)
ROI_plugin_add_status_code_to_db(solver, 1L, "Error", "Solution is undefined.", 1L)

dbind <- function(...) {
  .dbind <- function(x, y) {
    A <- simple_triplet_zero_matrix(NROW(x), NCOL(y))
    B <- simple_triplet_zero_matrix(NROW(y), NCOL(x))
    rbind(cbind(x, A), cbind(B, y))
  }
  Reduce(.dbind, list(...))
}

glmnet_lasso <- function(x, y, lambda) {
  stzm <- simple_triplet_zero_matrix
  stdm <- simple_triplet_diag_matrix
  m <- NROW(x); n <- NCOL(x)
  Q0 <- dbind(stzm(n), stdm(1, m), stzm(n))
  a0 <- c(b = double(n), g = double(m), t = lambda * rep(1, n))
  op <- OP(objective = Q_objective(Q = Q0, L = a0))
  ## y - X %*% beta = gamma  <=>  X %*% beta + gamma = y
  A1 <- cbind(x, stdm(1, m), stzm(m, n))
  LC1 <- L_constraint(A1, eq(m), y)
  ##  -t <= beta  <=>  0 <= beta + t
  A2 <- cbind(stdm(1, n), stzm(n, m), stdm(1, n))
  LC2 <- L_constraint(A2, geq(n), double(n))
  ##   beta <= t  <=>  beta - t <= 0
  A3 <- cbind(stdm(1, n), stzm(n, m), stdm(-1, n))
  LC3 <- L_constraint(A3, leq(n), double(n))
  constraints(op) <- rbind(LC1, LC2, LC3)
  bounds(op) <- V_bound(ld = -Inf, nobj = ncol(Q0))
  op
}
op <- glmnet_lasso(x, y, 1)

(qp0 <- ROI_solve(op, "glmnet3"))

mo <- glmnet(x = constraints(op)[1]$L[1:sum(constraints(op)[[3]]!=0),1:(sum(constraints(op)[[3]]==0)/2)], y=constraints(op)[[3]][1:sum(constraints(op)[[3]]!=0)], family = "gaussian", alpha = 1, lambda = op$objective[2]$L[op$objective[2]$L$ncol], 
             intercept = FALSE, standardize = FALSE)
