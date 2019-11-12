
calculate_3d_covar <-function(predictions)
{
  nmethods=dim(predictions)[2]
  npos=0
  nneg=0
  covar_tensor=array(0,c(nmethods,nmethods,nmethods))
  for(i in 1:nmethods)
  {
    nmethod1=predictions[,i]
    nmethod1=nmethod1-mean(nmethod1)
    for(j in 1:nmethods)
    {
      nmethod2=predictions[,j]
      nmethod2=nmethod2-mean(nmethod2)
      for(k in 1:nmethods)
      {


        nmethod3=predictions[,k]
        nmethod3=nmethod3-mean(nmethod3)
        ## This is the covariance of the three vectors
        d_cov=(1/length(nmethod1))*sum(nmethod1*nmethod2*nmethod3)

        ## this is 6 way symetric so no need to recalculate
        ## 6 times
        covar_tensor[i,j,k]=d_cov
        covar_tensor[i,k,j]=d_cov
        covar_tensor[j,i,k]=d_cov
        covar_tensor[k,i,j]=d_cov
        covar_tensor[k,j,i]=d_cov
        covar_tensor[j,k,i]=d_cov

        if(length(union(i,union(j,k)))==3)
        {
          if (d_cov>=0)
          {
            npos=npos+1
          }
          else
          {
            nneg=nneg+1
          }
        }

      }
    }
  }
  return(list(covar_tensor,npos,nneg))

}


# ================================================================
# ================================================================
# Estimate eigenvector and eigenvalue of rank one matrix
# using tensor recursion
# ================================================================
# ================================================================


find_eig_tensor <- function(covTensor, tol=50, full_out=FALSE)
{

  #tolerance could be an integer or tolerance



  t=rTensor::as.tensor(covTensor)

  if(tol<1)
  {

    eig=c(0,0)
    i=2

    error_now=0

    while( (error_now>tol) || (i<4))
    {
      i=i+1
      output=tensor_svd(t)
      t=output$t
      eig[i]=output$eig
      error_now=abs(eig[i]-eig[i-1])/max(abs(eig[i-1]),1)

    }
  }

  if(tol>1)
  {
    eig=c(0,0)
    for(i in 3:(tol+2))
    {
      output=tensor_svd(t)
      t=output$t
      eig[i]=output$eig
      error_now=abs(eig[i]-eig[i-1])/max(abs(eig[i-1]),1)
      if(error_now<0.00001)
      {
        break
      }

    }
  }

  return(list(t=t,eigen=eig))
}

####
####
tensor_svd <- function(covTensor)
{


  t=rTensor::hosvd(covTensor)

  nmethods=dim(covTensor)[3]

  ### Only maintain the 1D

  xx=t$Z[1,1,1]
  eig_1=norm(t$Z[,,1]@data,type="F")

  dzero=rTensor::as.tensor(array(0,c(nmethods,nmethods,nmethods)))
  dzero[1,1,1]=xx

  t$Z=dzero


  ####


  ## get the 1 d estimate
  t=rTensor::ttm(rTensor::ttm(rTensor::ttm(t$Z,t$U[[1]],1),t$U[[2]],2),t$U[[3]],3)
  #change the diagonals

  for(i in 1:nmethods)
  {
    for(j in 1:nmethods)
    {
      for(k in 1:nmethods)
      {
        if(length(union(i,union(j,k)))<3)
          covTensor[i,j,k]=t[i,j,k]
      }
    }
  }

  output=list(t=covTensor,eig=eig_1)
  return(output)
}


##Predictions here is the matrix where rows are patients columns
## are methods

###########
###########

calculate_norm <- function(tensor_eig,cov_eig)
{
  alpha=tensor_eig^2/cov_eig^3
  p=c(1,-4-alpha,4+alpha)
  roots=polyroot(p)
  ### Roots of polynomial p are rho and 1-rho where
  ### rho is the prevalence to get nomr delta we need
  ### one more operataion
  root1=Re(roots[1])
  root2=Re(roots[2])

  norm_delta=cov_eig*(alpha+4)#/(root1*root2)
  norm_delta=sqrt(norm_delta)
  return(c(norm_delta,root1,root2))

}
