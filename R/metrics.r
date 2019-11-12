
#' Calculate actual performance of methods
#'
#' This function calculates actual performance of a given set of learners using
#' a user-defined gold standard. If problem type is binary this calculates
#' balanced accuracy, if it is of ranks then the performance measure is AUC.
#' @param summa A summa object created by the summa() function
#' @param true_val A numeric vector of samples corresponding to true labels with 1
#' representing positive class and -1 representing negative class.
#' @return An object of class summa
#' @examples
#' load(system.file("extdata", "example_rank.Rdata", package = "UnsupervisedEnsemble"))
#' summa=summa(rank_matrix,"rank")
#' summa=calculate_performance(summa,actual_labels,"rank")
#' @export


calculate_performance<-function(summa,true_val,nrandom=50)
{
if(summa@type == "binary")
{
predictions=summa@predictions
performance_all=apply(predictions,2,calculate_ba,true_val=true_val)
summa@actual_performance=as.matrix(performance_all)
summa@woc=calculate_ba(apply(predictions,1,mean),true_val)
summa@summa=calculate_ba(predictions %*% summa@weights,true_val)
summa@summa_combined_performance=get_combination_scores(summa,true_val)
summa@woc_combined_performance=get_random_woc_performance(summa,true_val,nrandom)
}
else if(summa@type == "rank")
{
predictions=summa@predictions
performance_all=apply(predictions,2,calculate_auc,true_val=true_val)
summa@actual_performance=as.matrix(performance_all)
summa@woc=calculate_auc(apply(predictions,1,mean),true_val)
summa@summa=calculate_auc(as.vector(predictions %*% summa@weights),true_val)
summa@summa_combined_performance=get_combination_scores(summa,true_val)
summa@woc_combined_performance=get_random_woc_performance(summa,true_val,nrandom)
}

return(summa)
}


get_combination_scores <- function(summa, true_val)
{
type=summa@type
nmethods=summa@nmethods
summa_combined_performance=as.vector(matrix(0,nmethods,1))
idx_1=which.max(summa@estimated_performance)
summa_combined_performance[1]=summa@actual_performance[idx_1]
idx=order(abs(summa@estimated_performance),decreasing=TRUE)
predictions=summa@predictions
if(type=="binary")
{

for(i in 2:length(idx))
{
summa_combined_performance[i]=calculate_ba(as.vector(predictions[,idx[1:i]] %*% summa@weights[idx[1:i]]),true_val)
#woc_combined_performance[i]=calculate_ba(apply(predictions[,idx[1:i]],1,mean) ,true_val)

}
}

if(type=="rank")
{
  for(i in 2:length(idx))
  {
summa_combined_performance[i]=calculate_auc(as.vector(predictions[,idx[1:i]] %*% summa@weights[idx[1:i]]),true_val)
#woc_combined_performance[i]=calculate_auc(apply(predictions[,idx[1:i]],1,mean) ,true_val)
}
}

return(summa_combined_performance)
}

get_random_woc_performance <- function(summa,true_val,nrandom=50)
{
  type=summa@type
  nmethods=summa@nmethods
  predictions=summa@predictions
  woc_combined_performance=matrix(0,nmethods*nrandom,2)
  colnames(woc_combined_performance)=c("n_methods","performance")
  kk=1

  for(i in 1:nmethods)
  {
  for (j in 1:nrandom)
  {
  woc_combined_performance[kk,1]=i
  idx=sample(1:nmethods,i,replace=FALSE)

  if(type=="binary")
  {
  woc_combined_performance[kk,2]=calculate_ba(apply(as.matrix(predictions[,idx]),1,mean) ,true_val)
  }

  if(type=="rank")
  {
  woc_combined_performance[kk,2]=calculate_auc(apply(as.matrix(predictions[,idx]),1,mean) ,true_val)
  }

  kk=kk+1
  }

  }

  return(woc_combined_performance)
}


calculate_auc <- function(predict, true_val)
{
roc_obj <- pROC::roc(true_val, predict,direction=">")
auc=roc_obj$auc[1]
return(auc)
}

# compute balanced accuracy
calculate_ba <- function(predict, true_val)
{
  pp <- which(predict >= 0)
  p  <- which(true_val >= 0)
  pn <- which(predict < 0)
  n   <- which(true_val < 0)
  tpr=length(intersect(pp,p))/length(p)
  tnr=length(intersect(pn,n))/length(n)
  ba=0.5*(tpr + tnr)
  return(ba)
}

get_SE <- function(x) {
  v <- c( mean(x) - sd(x)/sqrt(length(x)), mean(x), mean(x) + sd(x)/sqrt(length(x)))
  names(v) <- c("ymin",  "middle",  "ymax")
  v
}

# ===============================================================
# ===============================================================
# rates
# ===============================================================
# ===============================================================

#####  Calculate AUC
#####
#####

