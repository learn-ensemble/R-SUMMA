# ================================================================
# ================================================================
# Estimate eigenvector and eigenvalue of rank one matrix
# ================================================================
# ================================================================
# Input:
# m method x m method covariance matrix
# Returns:
# data frame of eigenvalue(s) (value(s)) and eigenvector(s) (vector(s)):
# if full_out is FALSE, the approximate rank one eigenvalue 'value' and eigenvector 'vector'
# if full_out is TRUE, the entire spectrum of the approximate rank one matrix


find_eig <- function(covMatrix, tol=1e-6,niter_max=10000)
{
  eig <- eigen(covMatrix, symmetric=TRUE)
  eigen_all=vector()
  l <- 0
  # Iterate procedure until the lead eigenvalue does not, effectively, change.
  iter=0
  while (abs(l - eig$values[1]) > tol & (iter<niter_max))
  {
    l <- eig$values[1]
    # estimate the rank one matrix
    r <- get_rank_one_matrix(covMatrix, sqrt(l[1])*eig$vectors[, 1])
    eig <- eigen(r, symmetric=TRUE)
    eigen_all=c(eigen_all,eig$values[1])
    iter=iter+1
  }
  # if (full_out == FALSE)
  # {
  #   # we assume that the majority of methods have a positive eigenvector value
  #   # apply this assumption to the solution of the Rank 1 matrix
  #   sign_assumption <- sign(length(eig$values)/2 - sum(eig$vectors[, 1] < 0))
  #   if (sign_assumption < 0)
  #   {
  #     print('Eigenvector elements are a majority negative, multiply by negative 1')
  #   }
  #   output <- data.frame('value'=eig$values[1],
  #                        'vector'=sign_assumption * eig$vectors[, 1],
  #                        'method_id'=rownames(covMatrix))
  # }
  # else
  # {
  #   output <- list(eig=eig,eigen_all=eigen_all)
  # }

  output <- list(eig=eig,eigen_all=eigen_all)
  return(output)
}

# ================================================================
# ================================================================
# estimates the rank one matrix, used in the 'find_eig' function

get_rank_one_matrix <- function(c, r)
{
  return(c - diag(diag(c)) + diag(r^2))
}

# ================================================================
# ================================================================
# Compute maximum likelihood score and predictions
# ================================================================
# ================================================================
# scores
