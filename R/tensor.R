


#' Construct an unsupervised ensemble
#'
#'
#' This is the main summa function, which estimates a ranking of a given number
#' of methods in an unsupervised setting. It also constructs an unsupervised
#' ensemble of the given methods.
#' If type is "binary" we implement
#' the SML method by Parisi et. al. for more detail see [1,2].
#' @param type A character specifying the nature of the analysis and can be
#' either "binary" or "rank".
#' @param predictions A matrix of size samples by methods. If
#' type="binary", this should consist of binary values. By convention
#' samples belonging to the negative class should be
#' denoted by -1 and to the positive denoted by 1. If type="rank" each
#' entry corresponds to a confidence score
#' (not necessarily ranked by the user) with samples having higher scores
#' are more likely to belong to the positive class
#' @return An object of class summa
#' @examples
#' load(system.file("extdata", "example_rank.Rdata", package = "UnsupervisedEnsemble"))
#' summa=summa(rank_matrix,"rank")
#' @references [1] F. Parisi et. al. "Ranking and combining multiple predictors without
#' labeled data". doi:10.1073, Proceedings of the National Academy of Sciences, 2014.
#' @references [2] A. Jaffe et. al. "Estimating the accuracies of multiple classifiers
#' without labeled data. In AISTATS (Vol. 2,p. 4).
#' @export
#'
summa <- function(predictions,type,use_tensor=TRUE) {

  ##### This is the main summa function
  ##### You only need to input the prediction matrix
  ##### and the type is a character of binary or rank
  ##### If binary please let positive class be 1
  ##### Negative class be
  ##### If it is of type rank then the patients with higher scores
  ##### are most probable to be of positive class.

  predictions=as.matrix(predictions)

  if(type=="binary")
  {
    nclasses=length(unique(as.vector(predictions)))
    if(nclasses != 2)
    {
      print("The Predictions shoul be binary valued!")
      return(0)
    }
    covar_matrix=cov(predictions)
    tensorr_list=calculate_3d_covar(predictions)
    npos_tensor=tensorr_list[[2]]
    nneg_tensor=tensorr_list[[3]]
    covar_tensor=tensorr_list[[1]]
  }


  else
  {
    #convert to ranks
    predictions=apply(predictions,2,convert_to_rank)
    covar_matrix=cov(predictions)
    tensorr_list=calculate_3d_covar(predictions)
    npos_tensor=tensorr_list[[2]]
    nneg_tensor=tensorr_list[[3]]
    covar_tensor=tensorr_list[[1]]
  }


  covar_list=find_eig(covar_matrix)





  eigenvector=covar_list$eig$vectors[,1]
  eigen_2d_all=covar_list$eigen_all
  niterations_2d=length(eigen_2d_all)
  eigen_2d=eigen_2d_all[niterations_2d]
  error_covariance=100*abs(eigen_2d_all[niterations_2d]-eigen_2d_all[niterations_2d-1])/abs(eigen_2d_all[niterations_2d])
  if(error_covariance<0.01)
  {
    error_covariance="< 0.01"
  }
  else
  {
    error_covariance=round(error_covariance,2)
  }

 # plot(1:niterations_2d,eigen_2d_all,type="b",pch=16,
  #     xlab="Iterations",ylab="Eigenvalue Covariance Matrix",
  #      main=paste0("Convergence error (",error_covariance," % )"))


  npositive_eigen=length(which(eigenvector>=0))
  nnegative_eigen=length(which(eigenvector<0))

  if( nnegative_eigen > npositive_eigen)
  {
    eigenvector=-eigenvector
  }

  eigenvector=(eigenvector/sqrt(sum(eigenvector^2)))

  #' @export
  summa=new("summa",predictions=predictions,type=type)
  summa@covariance_matrix=covar_matrix
  summa@covariance_tensor=covar_tensor
  #summa@weights=eigenvector/sqrt(sum(eigenvector^2))
  summa@weights=eigenvector/sqrt(sum(eigenvector^2))
  summa@nmethods=dim(predictions)[2]
  summa@nsamples=dim(predictions)[1]

  if(summa@type=="binary")
  {
    summa@estimated_rank=as.vector(predictions %*% summa@weights)
    summa@estimated_label=as.vector(predictions %*% summa@weights)
    summa@estimated_label[summa@estimated_label>=0]=1
    summa@estimated_label[summa@estimated_label<0]=-1
  }

  if(summa@type=="rank")
  {
    summa@estimated_rank=as.vector(((0.5*(summa@nsamples+1)-summa@predictions) %*% summa@weights))
    # npositive=floor(summa@nsamples * summa@estimated_prevelance)
    # idx_summa=order(summa@estimated_rank,decreasing=FALSE)
    summa@estimated_label=as.vector((0.5*(summa@nsamples+1)-summa@predictions) %*% summa@weights)
    summa@estimated_label[summa@estimated_label>=0]=1
    summa@estimated_label[summa@estimated_label<0]=-1
    #summa@estimated_label[idx_summa[1:npositive]]=1
    #summa@estimated_label[idx_summa[(npositive+1):summa@nsamples]]=-1
  }


  ##### tensor



  if(use_tensor==TRUE)
  {
  tensor_list=find_eig_tensor(covar_tensor)
  eigen_3d_all=tensor_list$eigen
  niterations_3d=length(eigen_3d_all)
  eigen_3d=eigen_3d_all[niterations_3d]


  error_tensor=100*abs(eigen_3d_all[niterations_3d]-eigen_3d_all[niterations_3d-1])/abs(eigen_3d_all[niterations_3d])

  if(error_tensor<0.01)
  {
  error_tensor="< 0.01"
  }
  else
  {
  error_tensor=round(error_tensor,2)
  }




  #plot(3:niterations_3d,eigen_3d_all[3:niterations_3d],type="b",pch=16,
  #     xlab="Iterations",ylab="Eigenvalue Covariance Tensor",
  #     main=paste0("Convergence error (",error_tensor," % )"))

  ### Check the number of positive eigenvector elements
  ### and convert



  ### Calculate the norm of the eigenvector
  ### and then use it to calculate the eigenvector
  ###

  norm_eigenvector=calculate_norm(eigen_3d,eigen_2d)

  #### This gives 3 outputs
  ####
  root1=norm_eigenvector[2]
  root2=norm_eigenvector[3]
  norm_eigenvector=norm_eigenvector[1]

  eigenvector=(eigenvector/sqrt(sum(eigenvector^2)))*norm_eigenvector


  # if(summa@type=="binary")
  # {
  # eigenvector=eigenvector #0.5
  # }

  #summa@eigenvalue_matrix=eigen_2d
  #summa@eigenvalue_tensor=eigen_3d





  if(summa@type=="binary")
  {
    summa@estimated_performance=0.25*as.matrix(eigenvector)+0.5
    summa@estimated_rank=as.vector(predictions %*% summa@weights)
    summa@estimated_label=as.vector(predictions %*% summa@weights)
    summa@estimated_label[summa@estimated_label>=0]=1
    summa@estimated_label[summa@estimated_label<0]=-1
    summa@estimated_prevelance=(1-eigen_3d/(root1*root2*(norm_eigenvector)^3))/2
    if(npos_tensor < nneg_tensor)
    {
      summa@estimated_prevelance=1-summa@estimated_prevelance
    }

  }

  if(summa@type=="rank")
  {
    summa@estimated_performance=as.matrix((eigenvector/summa@nsamples+0.5))
    summa@estimated_rank=as.vector(((0.5*summa@nsamples)-summa@predictions) %*% summa@weights)
    npositive=floor(summa@nsamples * summa@estimated_prevelance)
    idx_summa=order(summa@estimated_rank,decreasing=FALSE)
    summa@estimated_label=as.vector(((0.5*summa@nsamples)-summa@predictions) %*% summa@weights)
    summa@estimated_label[summa@estimated_label>=0]=1
    summa@estimated_label[summa@estimated_label<0]=-1
    summa@estimated_prevelance=(1-eigen_3d/(root1*root2*(norm_eigenvector)^3))/2
    if(npos_tensor > nneg_tensor)
    {
      summa@estimated_prevelance=1-summa@estimated_prevelance
    }
    #summa@estimated_label[idx_summa[1:npositive]]=1
    #summa@estimated_label[idx_summa[(npositive+1):summa@nsamples]]=-1
  }

}

  return(summa)
}
## This function converts scores to ranks if they have already not been
## done.
convert_to_rank <- function(prediction) {
  #idx_order=order(prediction,decreasing= TRUE)
  #prediction[idx_order]=1:length(idx_order)
  prediction=rank(-prediction,ties.method = "random")
  return(prediction)
}

#' An S4 class to represent the summa output.
#' @slot predictions matrix of predictions
#' @slot covariance_matrix Square matrix of size learners containing
#' covariance of learners
#' @slot covariance_tensor Cube tensor of size learners containing
#' the covariance of three learners
#' @slot weights Contains weights for each method
#' @slot nsamples Number representing number of samples
#' @slot nmethods Number representing number of learners
#' @slot woc A numeric vector of samples constructed by weighting each
#' learner equally
#' @slot summa A numeric vector of samples constructed by weighting each
#' learner by their estimated performance
#' @slot type A user defined character vector describing the data under analysis
#' @slot actual_performance A numeric vector of learnears representing actual
#' performance of the learners, for binary data this is balanced accuracy
#' for ranked data this is AUC
#' @slot estimated_performance A numeric vector of learners representing
#' the summa estimated performance of learners
#' @slot estimated_prevelance A number corresponding to estimated prevelance
#' @slot sampe_rank A numeric vector ranking each sample
#' by decreased confidence
#' belonging to positive class using the summa ensemble
#' @export

setClass("summa",
         slots=c(
           predictions="array",  # a single string
           covariance_matrix="array",
           covariance_tensor="array",
           weights="numeric",
           nsamples="numeric",
           nmethods="numeric",
           estimated_label="numeric",
           woc="numeric",
           summa="numeric",
           type="character",        # an integer vector of length N
           actual_performance="array",
           estimated_performance="array",
           estimated_prevelance="numeric",
           estimated_rank="numeric",
           summa_combined_performance="numeric",
           woc_combined_performance="array"
         )
)


setClass("summa",
         slots=c(
           predictions="array",  # a single string
           covariance_matrix="array",
           covariance_tensor="array",
           weights="numeric",
           nsamples="numeric",
           nmethods="numeric",
           estimated_label="numeric",
           woc="numeric",
           summa="numeric",
           type="character",        # an integer vector of length N
           actual_performance="array",
           estimated_performance="array",
           estimated_prevelance="numeric",
           estimated_rank="numeric",
           summa_combined_performance="numeric",
           woc_combined_performance="array"
         )
)



get_scores <- function(X, weights)
{
  weighted_vote <- apply(X, 1, function(x, w) 2*(x - 0.5)*w, w=weights)
  return(apply(weighted_vote, 2, sum))
}

# ================================================================
# ================================================================
# prediction

mle_class <- function(X, weights)
{
  sc <- get_scores(X, weights)
  return(0.5*(sign(sc) + 1))
}



