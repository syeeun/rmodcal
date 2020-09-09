
#' Variance estimation of adjusted-weighted estimates of the expected number of cases
#' 
#' Internal use only.
#' 
#' @keywords ModelCal
#' @export 
var.Eadj = function(risk, auxinfo, ind.ph2, event, incl.prob, gcal, omega = NULL){
  
  n = length(ind.ph2)
  w0 = 1/incl.prob
  noncase_mar = (event != 1)
  
  tmp_risk = tmp_gcal = rep(0, n)
  tmp_risk[ind.ph2 == 1] = risk
  tmp_gcal[ind.ph2 == 1] = gcal
  risk = tmp_risk
  gcal = tmp_gcal
  
  aux = cbind(1, auxinfo)
  
  Eqd21 = colSums(risk*ind.ph2*w0*gcal*aux)
  Eqd11 = t(aux) %*% diag(ind.ph2*w0*gcal) %*% aux
  # matrix(c(sum(ind.ph2*w0*gcal), sum(ind.ph2*w0*gcal*auxinfo),
  #                  sum(ind.ph2*w0*gcal*auxinfo), sum(ind.ph2*w0*gcal*auxinfo^2)), 2, 2)
  # 
  bhat = Eqd21 %*% solve(Eqd11) #2x1
  
  bX = tcrossprod(aux, bhat)
  rstar = gcal*(risk-bX)
  Vwrstar = ind.ph2*w0*rstar
  
  exp_var = sum((1-incl.prob)*ind.ph2*(w0*rstar)^2)
  if(!is.null(omega)){
    Vwrstar_ncc = Vwrstar[noncase_mar == 1 & ind.ph2 == 1]
    exp_var = as.numeric(t(Vwrstar_ncc) %*% omega %*% Vwrstar_ncc)
  }
  
  var_exp = n/(n-1)*sum(ind.ph2*w0*(rstar-mean(Vwrstar))^2) + n/(n-1)*sum((bX-mean(bX))*ind.ph2*w0*(rstar-mean(Vwrstar))) + n*var(bX)
  varE = exp_var + var_exp
  
  return(as.numeric(varE))
}

#' Standard error of the estimated observed-to-expected ratio
#' 
#' Internal use only.
#' 
#' @keywords ModelCal
#' @export 
ratiosd = function(O, E, vE, covOE){
  sqrt(O/E^2 + O^2/E^4*vE - 2*O/E^3*covOE)
}

#' Covariance between nested case-control samples
#' 
#' Internal use only. Source from modelbasedVar() in 'multipleNCC' package.
#' 
#' @keywords ModelCal
#' @export 
wcov.ncc = function(survtime, samplestat, m, psample, left.time = 0, match.var = 0, match.int = 0)  {
  
  n.cohort = length(survtime)
  status = rep(0,n.cohort)
  status[samplestat > 1] = 1
  o = order(survtime)
  samplestat = samplestat[o]
  status = status[o]
  survtime = survtime[o]
  psample = psample[o]
  o2 = order(o[samplestat==1 & status == 0])
  
  endpoint = length(unique(status))-1-1*any(is.na(status))
  
  samplestat = samplestat!=0
  renkont = sum(samplestat!=0 & status==0)
  failuretimes = survtime[which(status != 0)]
  tidcase<-survtime[status==endpoint]
  leftcase = left.time[status==endpoint]
  tidkont<-survtime[status==0 & samplestat!=0]
  leftkont = left.time[status==0 & samplestat!=0]
  
  #With matching without left truncation
  
  if(sum(left.time) == 0 & length(match.var) > 1)  {
    
    match.var = as.matrix(match.var)
    fail.match = as.matrix(match.var[which(status != 0 ),])
    kont.match = as.matrix(match.var[which(status == 0 & samplestat!=0),])
    nfail = rep(0,length(failuretimes))
    
    if(dim(match.var)[2] == 1)  {
      
      for(k in 1:length(failuretimes))  {
        nfail[k] = length(survtime[which(survtime >= failuretimes[k]&
                                           match.var[,1] >= fail.match[k,1]+match.int[1] & match.var[,1] <= fail.match[k,1]+match.int[2])])
      }
      
      H = 1-2*m/(nfail-1)+m*(m-1)/((nfail-1)*(nfail-2))
      H = H/(1-m/(nfail-1))^2
      
      Rho<-matrix(rep(0,renkont^2),ncol=renkont)
      for(i in 1:(renkont-1)) {
        for(j in ((i+1):renkont))  {
          ind = survtime[i] > failuretimes & survtime[j] > failuretimes &
            match.var[i,1] >= fail.match[,1]+match.int[1] & match.var[i,1] <= fail.match[,1]+match.int[2] &
            match.var[j,1] >= fail.match[,1]+match.int[1] & match.var[j,1] <= fail.match[,1]+match.int[2]
          Rho[i,j] = prod(H[ind])-1
        }
      }
    }
    
    
    if(dim(match.var)[2] == 2)  {
      
      for(k in 1:length(failuretimes))  {
        nfail[k] = length(survtime[which(survtime >= failuretimes[k] &
                                           match.var[,1] >= fail.match[k,1]+match.int[1] & match.var[,1] <= fail.match[k,1]+match.int[2] &
                                           match.var[,2] >= fail.match[k,2]+match.int[3] & match.var[,2] <= fail.match[k,2]+match.int[4])])
      }
      
      H = 1-2*m/(nfail-1)+m*(m-1)/((nfail-1)*(nfail-2))
      H = H/(1-m/(nfail-1))^2
      
      Rho<-matrix(rep(0,renkont^2),ncol=renkont)
      for(i in 1:(renkont-1)) {
        for(j in ((i+1):renkont))  {
          ind = survtime[i] > failuretimes & survtime[j] > failuretimes &
            match.var[i,1] >= fail.match[,1]+match.int[1] & match.var[i,1] <= fail.match[,1]+match.int[2] &
            match.var[j,1] >= fail.match[,1]+match.int[1] & match.var[j,1] <= fail.match[,1]+match.int[2] &
            match.var[i,2] >= fail.match[,2]+match.int[3] & match.var[i,2] <= fail.match[,2]+match.int[4] &
            match.var[j,2] >= fail.match[,2]+match.int[3] & match.var[j,2] <= fail.match[,2]+match.int[4]
          Rho[i,j] = prod(H[ind])-1
        }
      }
    }
    
    
    if(dim(match.var)[2] == 3)  {
      
      for(k in 1:length(failuretimes))  {
        nfail[k] = length(survtime[which(survtime >= failuretimes[k] &
                                           match.var[,1] >= fail.match[k,1]+match.int[1] & match.var[,1] <= fail.match[k,1]+match.int[2] &
                                           match.var[,2] >= fail.match[k,2]+match.int[3] & match.var[,2] <= fail.match[k,2]+match.int[4] &
                                           match.var[,3] >= fail.match[k,3]+match.int[5] & match.var[,3] <= fail.match[k,3]+match.int[6])])
      }
      
      H = 1-2*m/(nfail-1)+m*(m-1)/((nfail-1)*(nfail-2))
      H = H/(1-m/(nfail-1))^2
      
      Rho<-matrix(rep(0,renkont^2),ncol=renkont)
      for(i in 1:(renkont-1)) {
        for(j in ((i+1):renkont))  {
          ind = survtime[i] > failuretimes & survtime[j] > failuretimes &
            match.var[i,1] >= fail.match[,1]+match.int[1] & match.var[i,1] <= fail.match[,1]+match.int[2] &
            match.var[j,1] >= fail.match[,1]+match.int[1] & match.var[j,1] <= fail.match[,1]+match.int[2] &
            match.var[i,2] >= fail.match[,2]+match.int[3] & match.var[i,2] <= fail.match[,2]+match.int[4] &
            match.var[j,2] >= fail.match[,2]+match.int[3] & match.var[j,2] <= fail.match[,2]+match.int[4] &
            match.var[i,3] >= fail.match[,3]+match.int[5] & match.var[i,3] <= fail.match[,3]+match.int[6] &
            match.var[j,3] >= fail.match[,3]+match.int[5] & match.var[j,3] <= fail.match[,3]+match.int[6]
          Rho[i,j] = prod(H[ind])-1
        }
      }
    }
    
    
    if(dim(match.var)[2] > 3)  {
      
      tt = function(rad)  {
        sum(rad==1)==length(rad)
      }
      
      nfail = rep(0,length(failuretimes))
      ncont3 = function(M,T)  {
        mat = matrix(ncol = length(M)+1,nrow=dim(match.var)[1],0)
        mat[which(survtime >= T),1] = 1
        for(i in 1:length(M))  {
          mat[which(match.var[,i] >= M[i]+match.int[(2*i-1)] & match.var[,i] <= M[i]+match.int[(2*i)]),(i+1)] = 1
        }
        sum(apply(mat,1,tt))
      }
      
      for(k in 1:length(failuretimes))  {
        nfail[k] =  ncont3(fail.match[k,],failuretimes[k])
      }
      
      H = 1-2*m/(nfail-1)+m*(m-1)/((nfail-1)*(nfail-2))
      H = H/(1-m/(nfail-1))^2
      
      indeks3 = function(M1,T1,M2,T2)  {
        mat = matrix(ncol = length(M1)*2+2,nrow=dim(fail.match)[1],0)
        mat[which(tidcase <T1),1] = 1
        mat[which(tidcase <T2),2] = 1
        for(i in 1:length(M1))  {
          mat[which(M1[i] >= fail.match[,i]+match.int[(2*i-1)] & M1[i] <= fail.match[,i]+match.int[(2*i)]),(i+2)] = 1
          mat[which(M2[i] >= fail.match[,i]+match.int[(2*i-1)] & M2[i] <= fail.match[,i]+match.int[(2*i)]),(i+dim(match.var[2])+2)] = 1
        }
        apply(mat,1,tt)
      }
      
      
      Rho<-matrix(rep(0,renkont^2),ncol=renkont)
      for(i in 1:(renkont-1)) {
        for(j in ((i+1):renkont))  {
          ind = indeks3(kont.match[i,],tidkont[i],kont.match[j,],tidkont[j])
          Rho[i,j] = prod(H[ind])-1
        }
      }
    }
  }
  
  
  #With matching and left truncation
  if(sum(left.time) > 0 & length(match.var) > 1)  {
    
    match.var = as.matrix(match.var)
    fail.match = as.matrix(match.var[which(status != 0 ),])
    kont.match = as.matrix(match.var[which(status == 0 & samplestat!=0),])
    nfail = rep(0,length(failuretimes))
    if(dim(match.var)[2] == 1)  {
      
      for(k in 1:length(failuretimes))  {
        nfail[k] = length(survtime[which(survtime >= failuretimes[k] & left.time < failuretimes[k] &
                                           match.var[,1] >= fail.match[k,1]+match.int[1] & match.var[,1] <= fail.match[k,1]+match.int[2])])
      }
      
      H = 1-2*m/(nfail-1)+m*(m-1)/((nfail-1)*(nfail-2))
      H = H/(1-m/(nfail-1))^2
      
      Rho<-matrix(rep(0,renkont^2),ncol=renkont)
      for(i in 1:(renkont-1)) {
        for(j in ((i+1):renkont))  {
          ind = survtime[i] > failuretimes & survtime[j] > failuretimes & left.time[i] < failuretimes & left.time[j] < failuretimes &
            match.var[i,1] >= fail.match[,1]+match.int[1] & match.var[i,1] <= fail.match[,1]+match.int[2] &
            match.var[j,1] >= fail.match[,1]+match.int[1] & match.var[j,1] <= fail.match[,1]+match.int[2]
          Rho[i,j] = prod(H[ind])-1
        }
      }
    }
    
    
    if(dim(match.var)[2] == 2)  {
      
      for(k in 1:length(failuretimes))  {
        nfail[k] = length(survtime[which(survtime >= failuretimes[k] & left.time < failuretimes[k] &
                                           match.var[,1] >= fail.match[k,1]+match.int[1] & match.var[,1] <= fail.match[k,1]+match.int[2] &
                                           match.var[,2] >= fail.match[k,2]+match.int[3] & match.var[,2] <= fail.match[k,2]+match.int[4])])
      }
      
      H = 1-2*m/(nfail-1)+m*(m-1)/((nfail-1)*(nfail-2))
      H = H/(1-m/(nfail-1))^2
      
      Rho<-matrix(rep(0,renkont^2),ncol=renkont)
      for(i in 1:(renkont-1)) {
        for(j in ((i+1):renkont))  {
          ind = survtime[i] > failuretimes & survtime[j] > failuretimes & left.time[i] < failuretimes & left.time[j] < failuretimes &
            match.var[i,1] >= fail.match[,1]+match.int[1] & match.var[i,1] <= fail.match[,1]+match.int[2] &
            match.var[j,1] >= fail.match[,1]+match.int[1] & match.var[j,1] <= fail.match[,1]+match.int[2] &
            match.var[i,2] >= fail.match[,2]+match.int[3] & match.var[i,2] <= fail.match[,2]+match.int[4] &
            match.var[j,2] >= fail.match[,2]+match.int[3] & match.var[j,2] <= fail.match[,2]+match.int[4]
          Rho[i,j] = prod(H[ind])-1
        }
      }
    }
    
    
    if(dim(match.var)[2] == 3)  {
      
      for(k in 1:length(failuretimes))  {
        nfail[k] = length(survtime[which(survtime >= failuretimes[k] & left.time < failuretimes[k] &
                                           match.var[,1] >= fail.match[k,1]+match.int[1] & match.var[,1] <= fail.match[k,1]+match.int[2] &
                                           match.var[,2] >= fail.match[k,2]+match.int[3] & match.var[,2] <= fail.match[k,2]+match.int[4] &
                                           match.var[,3] >= fail.match[k,3]+match.int[5] & match.var[,3] <= fail.match[k,3]+match.int[6])])
      }
      
      H = 1-2*m/(nfail-1)+m*(m-1)/((nfail-1)*(nfail-2))
      H = H/(1-m/(nfail-1))^2
      
      Rho<-matrix(rep(0,renkont^2),ncol=renkont)
      for(i in 1:(renkont-1)) {
        for(j in ((i+1):renkont))  {
          ind = survtime[i] > failuretimes & survtime[j] > failuretimes & left.time[i] < failuretimes & left.time[j] < failuretimes &
            match.var[i,1] >= fail.match[,1]+match.int[1] & match.var[i,1] <= fail.match[,1]+match.int[2] &
            match.var[j,1] >= fail.match[,1]+match.int[1] & match.var[j,1] <= fail.match[,1]+match.int[2] &
            match.var[i,2] >= fail.match[,2]+match.int[3] & match.var[i,2] <= fail.match[,2]+match.int[4] &
            match.var[j,2] >= fail.match[,2]+match.int[3] & match.var[j,2] <= fail.match[,2]+match.int[4] &
            match.var[i,3] >= fail.match[,3]+match.int[5] & match.var[i,3] <= fail.match[,3]+match.int[6] &
            match.var[j,3] >= fail.match[,3]+match.int[5] & match.var[j,3] <= fail.match[,3]+match.int[6]
          Rho[i,j] = prod(H[ind])-1
        }
      }
    }
    
    if(dim(match.var)[2] > 3)  {
      tt = function(rad)  {
        sum(rad==1)==length(rad)
      }
      
      nfail = rep(0,length(failuretimes))
      ncont4 = function(M,T,V)  {
        mat = matrix(ncol = length(M)+2,nrow=dim(match.var)[1],0)
        mat[which(survtime >= T),1] = 1
        mat[which(left.time < T),2] = 1
        for(i in 1:length(M))  {
          mat[which(match.var[,i] >= M[i]+match.int[(2*i-1)] & match.var[,i] <= M[i]+match.int[(2*i)]),(i+2)] = 1
        }
        sum(apply(mat,1,tt))
      }
      for(k in 1:length(failuretimes))  {
        nfail[k] =  ncont4(fail.match[k,],failuretimes[k],survtime[status!=0][k])
      }
      
      H = 1-2*m/(nfail-1)+m*(m-1)/((nfail-1)*(nfail-2))
      H = H/(1-m/(nfail-1))^2
      
      indeks4 = function(M1,T1,V1,M2,T2,V2)  {
        mat = matrix(ncol = length(M1)*2+4,nrow=dim(fail.match)[1],0)
        mat[which(tidcase <T1),1] = 1
        mat[which(tidcase <T2),2] = 1
        mat[which(leftcase < T1),3] = 1
        mat[which(leftcase < T2),4] = 1
        for(k in 1:length(M1))  {
          mat[which(M1[k] >= fail.match[,k]+match.int[(2*k-1)] & M1[k] <= fail.match[,k]+match.int[(2*k)]),(k+4)] = 1
          mat[which(M2[k] >= fail.match[,k]+match.int[(2*k-1)] & M2[k] <= fail.match[,k]+match.int[(2*k)]),(k+dim(match.var[2])+4)] = 1
        }
        apply(mat,1,tt)
      }
      
      Rho<-matrix(rep(0,renkont^2),ncol=renkont)
      for(i in 1:(renkont-1)) {
        for(j in ((i+1):renkont))  {
          ind = indeks4(kont.match[i,],tidkont[i],left.time[i],kont.match[j,],tidkont[j],left.time[i])
          Rho[i,j] = prod(H[ind])-1
        }
      }
    }
  }
  
  #With left truncation without matching
  if(sum(left.time) > 0 & length(match.var) == 1)  {
    
    nfail = rep(0,length(failuretimes))
    for(i in 1:length(failuretimes))  {
      nfail[i] =  length(survtime[which(survtime >= failuretimes[i] & left.time < failuretimes[i])])
    }
    
    H<-1-2*m/(nfail-1)+m*(m-1)/((nfail-1)*(nfail-2))
    H<-H/(1-m/(nfail-1))^2
    
    rho<-cumprod(H)-1
    indeksT = order(tidkont)
    indeksL = order(leftkont)
    
    Rho<-matrix(rep(0,renkont^2),ncol=renkont)
    for(i in (1:renkont-1)) {
      for(j in ((i+1):renkont))  {
        
        ind = 1:length(failuretimes)
        ind = ind[which(failuretimes >= leftkont[i] & failuretimes >= leftkont[j] & failuretimes<tidkont[i] & failuretimes<tidkont[j])]
        
        Rho[i,j] = prod(H[ind])-1
        
      }
    }
  }
  
  
  
  #Without left truncation and without matching
  if(sum(left.time) == 0 & length(match.var) == 1)  {
    nfail = rep(0,length(failuretimes))
    for(i in 1:length(failuretimes))  {
      nfail[i] =  length(survtime[which(survtime >= failuretimes[i])])
    }
    
    H<-1-2*m/(nfail-1)+m*(m-1)/((nfail-1)*(nfail-2))
    H<-H/(1-m/(nfail-1))^2
    
    rho = rep(0,length(H))
    rho<-cumprod(H)-1
    
    Rho<-matrix(rep(0,renkont^2),ncol=renkont)
    
    for (i in 1:(renkont-1)){
      index<-sum(tidkont[i]>tidcase)
      Rho[i,(i+1):renkont] = rho[index]
    }
  }
  Rhony = Rho + t(Rho)
  p0 = psample[status==0 & samplestat==1]
  # Q = matrix((1-p0)/p0,ncol=1)%*%matrix((1-p0)/p0,nrow=1)
  # R = Rhony*Q+diag((1-p0))
  # R = R[o2, o2]
  
  Vij = Rhony * tcrossprod(1-p0) + diag((1-p0)*(p0)) # diag = var(Vi) = pi(1-pi); offdiag = cov(Vi, Vj) = pij-pi*pj
  Pij = Rhony * tcrossprod(1-p0) + tcrossprod(p0); diag(Pij) <- p0 # diag = E(Vi) = pi; offdiag = E(Vi*Vj) = pij
  R = Vij/Pij
  R = R[o2, o2]
  
  return(R)
}
