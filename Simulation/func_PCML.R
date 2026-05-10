##directly nlmnib does not work


PCML <- function(beta, data, alpha1) {
  ####some functions
  norm_vec <- function(x) sqrt(sum(x^2))
  
  inv.matrix=function(mat) { 
    MinEig=0.001
    inv=try(solve(mat),silent=T)
    if (is.character(inv)==TRUE){ 
      dec=eigen(mat, symmetric=T)
      iv=1/(dec$values*(dec$values>=MinEig)+MinEig*(dec$values<MinEig))
      inv=dec$vectors%*%diag(iv)%*%t(dec$vectors)
      print("inv func used")
    }
    return(inv)
  }
  
  
  NANthenINF <- function(x) {
    ifelse(is.nan(x),Inf,x)
  }
  
  U <- function(beta){
    as.vector(x_aug %*% beta - x_aug[, 1:(q + 1)] %*% alpha1) * z
  }
  
  U_temp_f <- function(beta,gamma){
    temp <- x_aug %*% beta
    U_temp <- as.vector(temp - z %*% alpha1)*z
    U_temp - matrix(gamma,nrow=n,ncol=length(gamma),byrow=TRUE)
  }
  
  Q_hat <- function(rho,U) {
    mean(log(1-as.vector(U%*%rho)))
  }
  
  Q1_hat <- function(lambda_all,beta,gamma,rho,U) {
    mean((y - x_aug %*% beta)^2) / 2 + 
      sum(lambda_all * abs(gamma) / abs(gamma_tilde)^w) + Q_hat(rho,U)
  }
  
  Q1_hat_mingamma <- function(lambda_all,beta,gamma,rho,U) {
    sum(lambda_all * abs(gamma) / abs(gamma_tilde)^w) + Q_hat(rho,U)
  }
  
  log.lik=function(beta,data){ 
    y <- as.numeric(data[1])
    data <- c(1, as.numeric(data[-1]))
    - (y - data %*% beta)^2 / 2
  }
  
  score.log.lik=function(beta,data){ 
    y <- as.numeric(data[1])
    data <- c(1, as.numeric(data[-1]))
    data * c(y - data %*% beta)
  }
  
  hessian.log.lik=function(beta,data){ 
    data <- c(1, as.numeric(data[-1]))
    - data %*% t(data)
  }
  
  
  moment=function(beta,data) {
    data <- c(1, as.numeric(data[-1]))
    data[1:(q + 1)] * c(data %*% beta - data[1:(q + 1)] %*% alpha1)
  }
  
  moment.jacobian=function(beta,data) { 
    data <- c(1, as.numeric(data[-1]))
    data[1:(q + 1)] %*% t(data)
  }
  
  
  innerloop=function(data,beta,gamma){ 
    MaxIters=100       ## max number of iterations
    MaxStepIters=10    ## max number to find step length 
    MaxTol=10^(-4)     ## tolerence for convergence
    
    N=dim(data)[1]        ## number of subjects
    
    e=NULL  ## moment condition
    
    for (i in 1:N){ 
      e=rbind(e,matrix(moment(beta,data[i,])-gamma,1,)) 
    }
    
    m=dim(e)[2]   ## number of moment conditions from external data
    
    
    lambda=matrix(0,m,1)  ## initial of lambda
    denom=matrix(1,N,1)   ## initial of the denominator of p_i
    ret.inner=0           ## error type
    lik=0                 ## initial of likelihood ratio function
    N1=1/N
    
    
    for (it in 1:MaxIters) { 
      temp.lambda=lambda   ## record the lambda value from last iteration for later comparison with updated lambda
      
      u=e/(denom%*%matrix(1,1,m))
      graddd=matrix(apply(u,2,sum),m,1)  ## updated gradient
      hessiannn=t(u)%*%u  ## updated hessian
      direc=inv.matrix(hessiannn)%*%graddd  ## direction of Newton's method
      
      
      step=1  ## initial step length
      for (si in 1:MaxStepIters){ 
        ind.lik=0   
        lambda1=lambda-direc*step  
        denom1=1-e%*%lambda1
        if (min(denom1)>N1)  { 
          lik1=-sum(log(denom1))  
          if (lik1<lik) { ind.lik=1
          break
          }
        }
        step=step/2
      }
      
      #### if cannot find appropriate step length, then lambda will not be updated##
      if (ind.lik==0) {
        ret.inner=2   
        break
      }
      
      lik=lik1
      lambda=lambda1
      denom=denom1
      
      test1=max(abs(lambda-temp.lambda)) 
      if (test1<MaxTol) {break}
      if (it==MaxIters) {ret.inner=1}  ## stop without convergence
      
    } 
    
    #if (ret.inner==1) {print("Maximum Inner Loop Iterations Reached before Convergence")}
    
    temp=0 
    
    for (i in 1:N) { 
      temp=temp+moment.jacobian(beta=beta,data=data[i,])/denom[i,1] 
    }
    
    M.b=t(temp)%*%lambda  ## 1st derivative wrt beta
    u=e/(denom%*%matrix(1,1,m))
    hessiannn=t(u)%*%u  
    M.bb=-t(temp)%*%inv.matrix(hessiannn)%*%temp    ## 2nd derivative wrt beta
    
    return(list(rho=lambda,lik=lik,M.b=M.b,M.bb=M.bb,e=e))
  }
  
  outerloop_gamma <- function(data,beta,gamma,seqnum,lambda_all) { ## Outer Loop for gamma_k
    gamma_temp <- gamma
    beta_temp <- beta
    
    lambda_n <- lambda_all[seqnum]
    gamma1_temp <- gamma_temp
    gamma1_temp[seqnum] <- 0
    gamma_tilde_sub <- gamma_tilde[seqnum]
    rho_temp <- as.vector(innerloop(data,beta_temp,gamma1_temp)$rho)
    U_temp <- U_temp_f(beta_temp,gamma1_temp)
    delta_1_gamma1 <- 0
    for (i in 1:n) {
      delta_1_gamma1 <- delta_1_gamma1 + rho_temp[seqnum]/as.numeric(1-U_temp[i,]%*%rho_temp)
    }
    delta_1_gamma1 <- delta_1_gamma1/n
    if (abs(delta_1_gamma1) < lambda_n/(abs(gamma_tilde_sub)^w)) {
      gamma_temp[seqnum] <- 0
      gamma1_est <- 0
      l <- 0
      #print(paste0("gamma_", seqnum, " estimated as 0"))
    } else {
      gamma1_est <- 1
      
      l <- 0
      tau <- 1
      maxit <- 300
      sdelta <- 1
      epsilon <- 1e-4
      
      if (gamma_temp[seqnum]==0) {
        gamma_temp[seqnum] <- gamma_tilde_sub
      }
      
      while ( (norm_vec(sdelta)>=epsilon) & (l<maxit) ) {
        l <- l+1
        
        U_temp <- U_temp_f(beta_temp,gamma_temp)
        rho_temp <- as.vector(innerloop(data,beta_temp,gamma_temp)$rho)
        if (sum(is.nan(U_temp%*%rho_temp))>0) {
          print(paste0("When updating gamma_", seqnum, " sum(is.nan(U_temp%*%rho_temp))>0"))
        } else if (max(U_temp%*%rho_temp)>=1) {
          print(paste0("When updating gamma_", seqnum, " max(U_temp%*%rho_temp)>=1"))
          break
        }
        Q1_hat_mingamma_temp <- Q1_hat_mingamma(lambda_all,beta_temp, gamma_temp, rho_temp, U_temp)
        if (is.nan(Q1_hat_mingamma_temp)) {
          print(paste0("When updating gamma_", seqnum, " Q1_hat_mingamma_temp NAN"))
        }
        
        delta_1_gamma1 <- sign(gamma_temp[seqnum])*lambda_n/(abs(gamma_tilde_sub)^w)
        a1 <- 0
        a2 <- matrix(0,nrow=q + 1,ncol= q + 1)
        u_gamma <- rep(0,q + 1)
        u_gamma[seqnum] <- -1
        u_gamma <- as.vector(u_gamma)
        for (i in 1:n) {
          delta_1_gamma1 <- delta_1_gamma1 + rho_temp[seqnum]/n/as.numeric(1-U_temp[i,]%*%rho_temp)
          a1 <- a1 + u_gamma/n/as.numeric(1-U_temp[i,]%*%rho_temp)
          a2 <- a2 + U_temp[i,]%*%t(U_temp[i,])/n/as.numeric(1-U_temp[i,]%*%rho_temp)^2
        }
        if (is.character(try(solve(a2),silent=T))==TRUE){
          print("a2 singular matrix")
          delta_1_gamma1gamma1 <- t(a1)%*%ginv(a2)%*%a1
        } else {
          delta_1_gamma1gamma1 <- t(a1)%*%solve(a2)%*%a1
        }
        delta_1_gamma1gamma1 <- t(a1)%*%solve(a2)%*%a1
        delta_2_gamma1 <- delta_1_gamma1/delta_1_gamma1gamma1
        if (is.infinite(delta_2_gamma1)) {
          break
        }
        
        tau <- 1
        sdelta <- rep(0,q + 1)
        sdelta[seqnum] <- tau*delta_2_gamma1
        
        MaxStepIters <- 5    ## max number to find step length 
        
        for (it in 1:MaxStepIters) {
          ind.lik <- 0   ## indicator of a successful step length
          U_temp_temp <- U_temp_f(beta_temp,c(gamma_temp-sdelta))
          rho_temp_temp <- as.vector(innerloop(data,beta_temp,gamma_temp-sdelta)$rho)
          if (sum(is.nan(U_temp_temp%*%rho_temp_temp))==0) {
            if (max(U_temp_temp%*%rho_temp_temp)<1) {
              Q1_hat_mingamma_temp1 <- Q1_hat_mingamma(lambda_all,beta_temp,c(gamma_temp-sdelta), rho_temp_temp, U_temp_temp)
              if (!is.nan(Q1_hat_mingamma_temp1)) {
                if (Q1_hat_mingamma_temp1 < NANthenINF(Q1_hat_mingamma_temp)) {
                  ind.lik <- 1
                  break
                }
              }
            }  
          }
          tau <- tau/2
          sdelta[seqnum] <- tau*delta_2_gamma1
        }
        #if (ind.lik==0) {print(paste0("Fail to find the step length for gamma_", seqnum))}
        
        ### even if no step length is found, still update gamma using the smallest step length
        gamma_temp <- gamma_temp-sdelta
      }
      #print(paste0("gamma_", seqnum, " updated"))
    }
    
    return (list(gamma=gamma_temp))
  }
  
  outerloop_beta=function(data,beta,gamma){ ## Outer Loop for beta 
    MaxIters=300
    MaxStepIters=5
    MaxTol=10^(-4)
    
    ret.outer=0 ## error type
    no.length=0  ## number of failures finding step length
    no.improve=0  ## number of iterations with no improvement of lik
    N=dim(data)[1]
    
    part.ext=innerloop(data,beta,gamma)
    lik.int=0
    for (i in 1:N){ 
      lik.int=lik.int+log.lik(beta,data[i,]) 
    }
    lik=lik.int+part.ext$lik
    
    for (it in 1:MaxIters) { 
      temp.beta=beta
      grad.int=0
      hessian.int=0
      for (i in 1:N) { 
        grad.int=grad.int+matrix(score.log.lik(beta=beta,data=data[i,]),ncol = 1)
        hessian.int=hessian.int+hessian.log.lik(beta=beta,data=data[i,])
      }
      graddd=grad.int+part.ext$M.b
      hessiannn=hessian.int+part.ext$M.bb
      
      direc=inv.matrix(hessiannn)%*%graddd  ## direction of Newton's method
      
      step=1                    ## initial step length
      for (si in 1:MaxStepIters){ 
        ind.lik=0
        beta1=beta-direc*step
        part.ext=innerloop(data,beta1,gamma)
        lik.int=0
        for (i in 1:N){ 
          lik.int=lik.int+log.lik(beta1,data[i,]) 
        }
        lik1=lik.int+part.ext$lik
        if (lik1>lik){ 
          ind.lik=1
          break 
        }
        step=step/2
      }
      
      beta=beta1   ## even if no step length is found, still update beta using the smallest step length
      
      test1=max(abs(beta-temp.beta))
      if (test1<MaxTol){ break }
      
      if (ind.lik==0){ 
        no.length=no.length+1  ## count number of failures to find length
        if (no.length==200) { 
          ## if happens alot, then stop iter
          ret.outer=2
          break
        }
      }
      
      ## if too many iterations with no improvement of lik, stop###
      if (abs(lik1-lik)<MaxTol) { 
        no.improve=no.improve+1
        if (no.improve==30){ 
          ret.outer=3
          break
        }
      }
      
      
      lik=lik1
      
      if (it==MaxIters){ ret.outer=1 }
    } 
    #if (ret.outer==1){print("Maximum Outer Loop Iterations Reached before Convergence")}
    #if (ret.outer==2){print("Number of Iterations with No Step Length Exceeds 200")}
    #if (ret.outer==3){print("Number of Iterations with No Increase of Likelihood Exceeds 30")}
    
    return (list(beta=beta))
  }
  
  argmin_betagamma <- function(lambda_all,beta_temp_int,gamma_temp_int) {
    beta_temp <- beta_temp_int
    gamma_temp <- gamma_temp_int
    
    
    t <- 0 ## record the iterated times for the whole procedure
    maxit_t <- 300
    
    epsilon <- 1e-4
    beta_diff <- 1
    gamma_diff <- 1
    gamma_old <- gamma_temp
    
    while ( ((norm_vec(beta_diff)>=epsilon) | (norm_vec(gamma_diff)>=epsilon) | sum((gamma_old==0)!=(gamma_temp==0))>0 ) & (t<maxit_t) ) {
      t <- t + 1
      beta_old <- beta_temp
      gamma_old <- gamma_temp
      
      ### Update gamma component-wise sequentially
      for (seqnum in 1:(q + 1)) {
        gamma_temp_temp <- outerloop_gamma(data,beta_old,gamma_temp,seqnum,lambda_all)$gamma
        gamma_temp <- gamma_temp_temp
      }
      ### Update beta
      beta_temp <- as.vector((outerloop_beta(data,beta_old,gamma_temp))$beta)
      
      beta_diff <- beta_old-beta_temp
      gamma_diff <- gamma_old-gamma_temp
    }
    
    rho_temp <- as.vector(innerloop(data,beta_temp,gamma_temp)$rho)
    list(beta=beta_temp,gamma=gamma_temp,rho=rho_temp,itnum=t)
  }
  
  s_beta <- function(beta){
    c(y - x_aug %*% beta) * x_aug
  }
  
  G_i_beta <- function(i,beta){
    temp <- x_aug[i, ] %*% beta
    x_aug[i, 1:(q + 1)] %*% t(x_aug[i, ])
  }
  
  Sigma_beta <- function(beta,gamma){
    U_temp <- U_temp_f(beta,gamma)
    U_temp <- cbind(U_temp[,gamma!=0],U_temp[,gamma==0])
    bigOmega <- t(U_temp)%*%U_temp/n
    bigG <- matrix(0,ncol=length(beta),nrow=length(gamma))
    for (i in 1:n) {
      bigG <- bigG + G_i_beta(i,beta)
    }
    bigG <- bigG/n
    bigG <- rbind(bigG[gamma!=0,],bigG[gamma==0,])
    if (length(which(gamma!=0))>0 & length(which(gamma==0))>0) {
      if (length(which(gamma!=0))==1) {
        bigG <- cbind(bigG,rbind(-1,matrix(0,ncol=length(which(gamma!=0)),nrow=length(which(gamma==0)))))
      } else {
        bigG <- cbind(bigG,rbind(diag(rep(-1,length(which(gamma!=0)))),matrix(0,ncol=length(which(gamma!=0)),nrow=length(which(gamma==0)))))
      }
      bigS <- (t(s_beta(beta))%*%s_beta(beta))/n
      bigS <- rbind(cbind(bigS,matrix(0,nrow=length(beta),ncol=length(which(gamma!=0)))),matrix(0,nrow=length(which(gamma!=0)),ncol=length(beta)+length(which(gamma!=0))))
    } else if (length(which(gamma!=0))==0) {
      bigS <- (t(s_beta(beta))%*%s_beta(beta))/n
    } else if (length(which(gamma==0))==0) {
      bigG <- cbind(bigG,diag(rep(-1,length(which(gamma!=0)))))
      bigS <- (t(s_beta(beta))%*%s_beta(beta))/n
      bigS <- rbind(cbind(bigS,matrix(0,nrow=length(beta),ncol=length(which(gamma!=0)))),matrix(0,nrow=length(which(gamma!=0)),ncol=length(beta)+length(which(gamma!=0))))
    }
    solve(bigS+t(bigG)%*%solve(bigOmega)%*%bigG)
  }
  
  cov_psi_phi_sqrt <- function(beta,gamma){
    U_temp <- U_temp_f(beta,gamma)
    U_temp <- cbind(U_temp[,gamma!=0],U_temp[,gamma==0])
    bigOmega <- t(U_temp)%*%U_temp/n
    mathcal_S <- t(s_beta(beta))%*%s_beta(beta)/n
    sqrtm(rbind(cbind(bigOmega,matrix(0,nrow=q + 1,ncol=p + 1)),cbind(matrix(0,nrow=p + 1,ncol=q + 1),mathcal_S)))
  }
  
  A_0 <- function(beta,gamma){
    U_temp <- U_temp_f(beta,gamma)
    U_temp <- cbind(U_temp[,gamma!=0],U_temp[,gamma==0])
    bigOmega <- t(U_temp)%*%U_temp/n
    bigG <- matrix(0,ncol=length(beta),nrow=length(gamma))
    for (i in 1:n) {
      bigG <- bigG + G_i_beta(i,beta)
    }
    bigG <- bigG/n
    bigG <- rbind(bigG[gamma!=0,],bigG[gamma==0,])
    if (length(which(gamma!=0))>0 & length(which(gamma==0))>0) {
      if (length(which(gamma!=0))==1) {
        bigG <- cbind(bigG,rbind(-1,matrix(0,ncol=length(which(gamma!=0)),nrow=length(which(gamma==0)))))
      } else {
        bigG <- cbind(bigG,rbind(diag(rep(-1,length(which(gamma!=0)))),matrix(0,ncol=length(which(gamma!=0)),nrow=length(which(gamma==0)))))
      }
      bigS <- (t(s_beta(beta))%*%s_beta(beta))/n
      bigS <- rbind(cbind(bigS,matrix(0,nrow=length(beta),ncol=length(which(gamma!=0)))),matrix(0,nrow=length(which(gamma!=0)),ncol=length(beta)+length(which(gamma!=0))))
    } else if (length(which(gamma!=0))==0) {
      bigS <- (t(s_beta(beta))%*%s_beta(beta))/n
    } else if (length(which(gamma==0))==0) {
      bigG <- cbind(bigG,diag(rep(-1,length(which(gamma!=0)))))
      bigS <- (t(s_beta(beta))%*%s_beta(beta))/n
      bigS <- rbind(cbind(bigS,matrix(0,nrow=length(beta),ncol=length(which(gamma!=0)))),matrix(0,nrow=length(which(gamma!=0)),ncol=length(beta)+length(which(gamma!=0))))
    }
    A_0_left <- diag(1,length(gamma_tilde))-bigG%*%solve(bigS+t(bigG)%*%solve(bigOmega)%*%bigG)%*%t(bigG)%*%solve(bigOmega)
    A_0_right <- bigG%*%solve(bigS+t(bigG)%*%solve(bigOmega)%*%bigG)%*%rbind(diag(1,length(beta)),matrix(0,nrow=nrow(bigS)-length(beta),ncol=length(beta)))
    cbind(A_0_left,A_0_right)
  }
  
  Pi_0 <- function(beta,gamma){
    U_temp <- U_temp_f(beta,gamma)
    U_temp <- cbind(U_temp[,gamma!=0],U_temp[,gamma==0])
    bigOmega <- t(U_temp)%*%U_temp/n
    solve(bigOmega)%*%A_0(beta,gamma)%*%cov_psi_phi_sqrt(beta,gamma)
  }
  
  #####begin calculation
  y <- as.numeric(data[, 1])
  x_aug <- cbind(1, as.matrix(data[, -1], nrow = n))

  z <- x_aug[, 1:(q + 1)]
  
  gamma_tilde <- colMeans(U(beta))
  
  w <- 2 ## Specify w in the adaptive Lasso penalty function
  lambda <- n^{-1/2 - w/4}
  
  ### First, get the preliminary PCML estimation results with tuning parameter \lambda_n=n^{-1/2-w/4}=1/n
  lambda_try <- rep(1/n,length(gamma_tilde))
  hat_rhobetagamma_lambda_try <- argmin_betagamma(lambda_all=lambda_try,beta,gamma_tilde)
  
  
  ### Then, calculate "C" as shown in the section of Tuning Parameter Selection
  
  Pi_est_temp <- Pi_0(hat_rhobetagamma_lambda_try$beta,hat_rhobetagamma_lambda_try$gamma)
  Pi_est <- Pi_est_temp
  for (k in 1:length(gamma_tilde)) {
    Pi_est[k,] <- Pi_est_temp[which(c(which(hat_rhobetagamma_lambda_try$gamma!=0),which(hat_rhobetagamma_lambda_try$gamma==0))==k),]
  }
  
  c_est <- c()
  for (k in 1:length(gamma_tilde)) {
    c_est[k] <- norm_vec(Pi_est[k,])
  }
  
  ### Now, with the selected tuning paramter, we compute the PCML estimator
  lambda_best <- 1/n*c_est
  hat_rhobetagamma_best <- argmin_betagamma(lambda_all=lambda_best,beta,gamma_tilde)
  hat_rhobetagamma_best$beta
}
