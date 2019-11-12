#' Create synthetic data
#'
#' This function creates a synthetic data set as well as the actual labels
#' @param nsamples A number representing the number of samples
#' @param nclassifiers A number representing the number of classifiers
#' @param prevelance A number representing the percentage of positive samples
#' @param type A character either "binary" or "rank"
#' @return A list where the first element is the prediction matrix and the
#' second element is actual labels
#' @examples
#' data_binary=create_predictions(100,20,0.4,"binary")
#' data_rank=create_predictions(100,20,0.7,"rank")
#' @export
create_predictions <- function(nsamples, nclassifiers,prevelance,type)
{

nb1=floor(0.5*nclassifiers)
nb2=floor(0.4*nclassifiers)
nb3=nclassifiers-nb1-nb2;
balanced_accuracy_classifiers=rep(0,nclassifiers);
balanced_accuracy_classifiers[1:nb1]=0.2*runif(nb1)+0.4;
balanced_accuracy_classifiers[(nb1+1):(nb1+nb2)]=0.2*runif(nb2,1)+0.45
balanced_accuracy_classifiers[(nb1+nb2+1):(nclassifiers)]=0.1*runif(nb3,1)+0.6
sigma_classifiers=rep(0,nclassifiers)
k_classifiers=rep(0,nclassifiers)
mu_p_classifiers=rep(0,nclassifiers)
mu_n_classifiers=rep(0,nclassifiers)
sigma_d_classifiers=rep(0,nclassifiers)


PD=prevelance

PN=1-PD

npoints=20
npoints1=50
sigma_d=seq(0.01,10,length.out=npoints1)
mu_p=seq(1.9,3,length.out=npoints)
mu_n=seq(1,2.1,length.out=npoints)
k_thres=seq(1,3,length.out=npoints1)

baccuracy_all=array(0,c(npoints,npoints,npoints1,npoints1))

for(j in 1:npoints)
{
  for(k in 1:npoints)
  {
    for(l in 1:npoints1)
    {
      for(ll in 1:npoints1)
      {
     #   P_TP=(1-pnorm(k_thres[ll],mu_p[j],sigma_d[l]));
     #P_TN=(pnorm(k_thres[ll],mu_n[k],sigma_d[l]));
      baccuracy_all[j,k,l,ll]=pnorm((mu_p[j]-mu_n[k])/(sqrt(2)*sigma_d[l]))
      }
    }
  }
}

        for(i in 1:nclassifiers)
        {
          dummy=balanced_accuracy_classifiers[i]
      A=baccuracy_all-dummy
      minx=min(abs(A))
      indices=which(abs(A)==minx,arr.ind=TRUE)
      mu_p_classifiers[i]=mu_p[indices[1,1]]
      mu_n_classifiers[i]=mu_n[indices[1,2]]
      sigma_d_classifiers[i]=sigma_d[indices[1,3]]
      k_classifiers[i]=0.5*(mu_p_classifiers[i]+mu_n_classifiers[i])#k_thres[indices[1,4]]
        }

### We create classifiers now lets create
##  patients



nd=floor(nsamples*PD)
nn=nsamples-nd


patient_risk_matrix=matrix(0,nsamples,nclassifiers)

patient_prediction_matrix=matrix(0,nsamples,nclassifiers)


for (i in 1:nd)
{
for (j in 1:nclassifiers)
{
patient_risk_matrix[i,j]=sigma_d_classifiers[j]*rnorm(1)+mu_p_classifiers[j]
}
}



for (i in (nd+1):nsamples)
{
for (j in 1:nclassifiers)
{
patient_risk_matrix[i,j]=sigma_d_classifiers[j]*rnorm(1)+mu_n_classifiers[j]
}
}

actual_labels=c(rep(1,nd),rep(-1,nn))
patient_prediction_matrix=patient_risk_matrix


if(type=="binary")
{
  for (j in 1:nclassifiers)
  {
    dummy=patient_risk_matrix[,j]
    idx=which(dummy<k_classifiers[j])
    patient_risk_matrix[idx,j]=-1
    patient_risk_matrix[-idx,j]=1
  }
}

data=list(predictions=patient_risk_matrix,
          actual_labels=actual_labels)



}
