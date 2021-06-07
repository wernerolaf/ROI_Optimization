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


dbind <- function(...) {
  .dbind <- function(x, y) {
    A <- simple_triplet_zero_matrix(NROW(x), NCOL(y))
    B <- simple_triplet_zero_matrix(NROW(y), NCOL(x))
    rbind(cbind(x, A), cbind(B, y))
  }
  Reduce(.dbind, list(...))
}

qp_lasso <- function(x, y, lambda) {
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


op <- qp_lasso(x, y, lambda * NROW(x))
(qp1 <- ROI_solve(op, "qpoases"))
head(solution(qp1), 10) #the outputs similar to glmnet_beta

# I guess we need to make solver using:

# ROI_plugin_register_solver_method(signatures, solver, method, plugin = solver)


# from https://epub.wu.ac.at/5858/1/ROI_StatReport.pdf:

#Extending ROI with a new solver method can be split into three parts: 
#First, a function to be called by ROI has to be written. 
#Second, the function plus information about the function are added into the ROI solver registry. 
#Third, a mapping from the solver specific arguments and the status codes to their ROI counterpart has to be provided.


# 1.)
#In order to establish a connection between the OP and the solvers provided via plug-ins, 
#both are equipped with a signature. The signature captures all the information necessary to determine
#which solver is applicable to a given problem.

OP_signature(op)

#New signatures are created by the function ROI_plugin_make_signature(). 
#The following shows how to create the signature for the glpk solver, 
#where the objective and the constraints have to be linear. 
glpk_signature <- ROI_plugin_make_signature(objective = "L",
                                            constraints = "L", 
                                            types = c("C", "I", "B", "CI", "CB", "IB", "CIB"),
                                            bounds = c("X", "V"), maximum = c(TRUE, FALSE))

#Furthermore, this signature indicates that, the variable types are allowed to be:
#binary ("B"), integer ("I"), continuous ("C") or any combinations of them. 
#The bounds have to be variable bounds ("V") or no bounds at all encoded by "X". 
#The last argument maximum specifies that GLPK can find both maxima and minima.

# 2.) Writing a new solver method

#Any function supposed to add a solver to ROI has to take the arguments x and control, 
#where x is of class ‘OP’ and control a list containing the additional arguments. 
#Furthermore, the solution has to be canonicalized before it is returned. 
#The following shows the code from ROI.plugin.glpk for solving linear problems.

#Rglpk_solve_LP to funkcja z pakietu Rglpk i te obiekty z listy glpk to argumenty tej funkcji

glpk_solve_OP <- function(x, control = list()) {
    
  control$canonicalize_status <- FALSE
  glpk <- list(Rglpk_solve_LP, obj = terms(objective(x))[["L"]],
               mat = constraints(x)$L, dir = constraints(x)$dir,
               rhs = constraints(x)$rhs, bounds = bounds(x),
               types = types(x), max = maximum(x), control = control)
  mode(glpk) <- "call"
  if ( isTRUE(control$dry_run) )
    return(glpk)
  
  out <- eval(glpk)
  ROI_plugin_canonicalize_solution(solution = out$solution,
                                   optimum = out$optimum, 
                                   status = out$status, 
                                   solver = "glpk", 
                                   message = out)
  
    
}

ROI_plugin_register_solver_method(glpk_signature, "glpk", glpk_solve_OP)


# Proba z glmnet #
glmnet_signature = ROI_plugin_make_signature(objective = 'Q', 
                                           constraints = 'L', 
                                           bounds = 'X', 
                                           types = 'C',
                                           maximum = FALSE)

glmnet_solve_OP <- function(x, control = list()) {
  
  control$canonicalize_status <- FALSE
  glmnet <- list(glmnet::glmnet, obj = terms(objective(x))[["L"]],
               mat = constraints(x)$L, dir = constraints(x)$dir,
               rhs = constraints(x)$rhs, bounds = bounds(x),
               types = types(x), max = maximum(x), control = control)
  mode(glmnet) <- "call"
  if ( isTRUE(control$dry_run) )
    return(glmnet)
  
  out <- eval(glmnet)
  ROI_plugin_canonicalize_solution(solution = out$solution,
                                   optimum = out$optimum, 
                                   status = out$status, 
                                   solver = "glmnet", 
                                   message = out)
  
  
}

ROI_plugin_register_solver_method(glpk_signature, "glpk", glpk_solve_OP)

