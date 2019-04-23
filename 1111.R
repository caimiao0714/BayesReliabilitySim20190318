
#####calculate s，n
set.seed(123)
alpha = matrix(c(0.001, 0.05, 0.0001, 0.08), ncol = 2, byrow = TRUE)
#alpha = matrix(c(0.005, 0.05, 0.0005, 0.08), ncol = 2, byrow = TRUE)
#alpha = matrix(c(0.008, 0.05, 0.0008, 0.08), ncol = 2, byrow = TRUE)

w = c(45, 55)
t = matrix(seq(10, 40, 10), ncol = 2, byrow = TRUE)
#t21=t[2,1]-t[1,2]+t[1,2]*(exp(alpha[1,1]*w[1]))/(exp(alpha[1,1]*w[2]))
#t22=t[2,2]-t[1,2]+t[1,2]*(exp(alpha[2,1]*w[1]))/(exp(alpha[2,1]*w[2]))
tn = matrix(c(10,20,30,40),ncol = 2, byrow = TRUE)
K = c(10, 100, 1000)
k = 3


lambda = function(r = 1:2, j = 1:2) {
  alpha_r0 = alpha[r, 1]
  alpha_r1 = alpha[r, 2]
  return(alpha_r0*exp(alpha_r1*w[j]))
}
p = function(i = 1:2, j = 1:2){
  lam1j = lambda(1, j) # failure due to factor 1
  lam2j = lambda(2, j) # failure due to factor 2
  p0 = exp(-(lam1j + lam2j)*tn[j, i])
  p1 = (lam1j/(lam1j + lam2j))*(1 - p0)
  p2 = (lam2j/(lam1j + lam2j))*(1 - p0)
  return(c(p0, p1, p2))
}
Tmatrxi1 = t(rmultinom(1000, K[k], prob = p(1, 1)))
Tmatrxi2 = t(rmultinom(1000, K[k], prob = p(1, 2)))
Tmatrxi3 = t(rmultinom(1000, K[k], prob = p(2, 1)))
Tmatrxi4 = t(rmultinom(1000, K[k], prob = p(2, 2)))

errorAll1 = 0
biasAll1 = 0
errorAll2 = 0
biasAll2 = 0
errorAll3 = 0
biasAll3 = 0
errorAll4 = 0
biasAll4 = 0

for (nn in 1:1000){ 
  s = matrix(c(mean(Tmatrxi1[nn,1]), mean(Tmatrxi3[nn,1]), mean(Tmatrxi2[nn,1]), mean(Tmatrxi4[nn,1])), ncol = 2, byrow = TRUE)
  n1 = matrix(c(mean(Tmatrxi1[nn,2]), mean(Tmatrxi3[nn,2]), mean(Tmatrxi2[nn,2]), mean(Tmatrxi4[nn,2])), ncol = 2, byrow = TRUE)
  n2 = matrix(c(mean(Tmatrxi1[nn,3]), mean(Tmatrxi3[nn,3]), mean(Tmatrxi2[nn,3]), mean(Tmatrxi4[nn,3])), ncol = 2, byrow = TRUE)
  
  
  ###### EM 
  alpha0 = matrix(c(0.01, 0.1, 0.001, 0.2), ncol = 2, byrow = TRUE)
  w = c(45, 55)
  dis=2
  ite=0
  alphaAll = matrix(c(0.01, 0.1, 0.001, 0.2), ncol = 4, byrow = TRUE)
  disAll = dis
  while(dis>0.0046){
    t = matrix(seq(10, 40, 10), ncol = 2, byrow = TRUE)
    t21=t[2,1]-t[1,2]+t[1,2]*(exp(alpha0[1,1]*w[1]))/(exp(alpha0[1,1]*w[2]))
    t22=t[2,2]-t[1,2]+t[1,2]*(exp(alpha0[2,1]*w[1]))/(exp(alpha0[2,1]*w[2]))
    tn = matrix(c(10,20,30,40),ncol = 2, byrow = TRUE)
    lambda = function(r = 1:2, j = 1:2) {
      alpha_r0 = alpha0[r, 1]
      alpha_r1 = alpha0[r, 2]
      return(alpha_r0*exp(alpha_r1*w[j]))
    }
    
    tji1 = function(i = 1:2, j = 1:2){
      if (j==1){
        return (tn[1,i])
      } else if (j==2){
        return (tn[2,i]-tn[1,2]+tn[1,2]*lambda(1,1)/lambda(1,2))
        #return (tn[2,i])
        
      }
    }
    tji2 = function(i = 1:2, j = 1:2){
      if (j==1){
        return (tn[1,i])
      } else if (j==2){
        return (tn[2,i]-tn[1,2]+tn[1,2]*lambda(2,1)/lambda(2,2))
        #return (tn[2,i])
        
      }
    }
    temp1 = function(i=1:2,j=1:2){
      return((1/(lambda(1,j)+lambda(2,j)))-((tji1(i,j)*exp(-(lambda(1,j)+lambda(2,j))*tji1(i,j)))/(1-exp(-(lambda(1,j)+lambda(2,j)))*tji1(i,j))))
    }
    temp2 = function(i=1:2,j=1:2){
      return((1/(lambda(1,j)+lambda(2,j)))-(tji2(i,j)*exp(-(lambda(1,j)+lambda(2,j))*tji2(i,j)))/(1-exp(-(lambda(1,j)+lambda(2,j)))*tji2(i,j)))
    }
    
    Et1j = function(j=1:2){
      return (s[1,j]*(tji1(1,j)+1/lambda(1,j))+n1[1,j]*(temp1(1,j))+n2[1,j]*((1/lambda(1,j))+temp1(1,j))+s[2,j]*(tji1(2,j)+1/lambda(2,j))+n1[2,j]*(temp1(2,j))+n2[2,j]*((1/lambda(2,j))+temp1(2,j)))
    }
    
    Et2j = function(j=1:2){
      return (s[1,j]*(tji2(1,j)+1/lambda(1,j))+n2[1,j]*(temp2(1,j))+n1[1,j]*((1/lambda(2,j))+temp2(1,j))+s[2,j]*(tji2(2,j)+1/lambda(2,j))+n2[2,j]*(temp2(2,j))+n1[2,j]*((1/lambda(2,j))+temp2(2,j)))
    }
    C1 = w[1]-((w[1]+w[2])/2)
    C2 = w[2]-((w[1]+w[2])/2)
    C = c(C1,C2)
    temp3 = function(r=1:2,j=1:2){
      a=0
      if (r==1){
        a = (C[j]*exp(alpha0[r,2]*w[j])*Et1j(j))
      }else{
        a = (C[j]*exp(alpha0[r,2]*w[j])*Et2j(j))
      }
      return (a)
    }
    
    alpha011 = alpha0[1,2]-(temp3(1,1)+temp3(1,2))/(temp3(1,1)*w[1]+temp3(1,2)*w[2])
    alpha021 = alpha0[2,2]-(temp3(2,1)+temp3(2,2))/(temp3(2,1)*w[1]+temp3(2,2)*w[2])
    alpha010 = (4*K[k])/((exp(alpha011*w[1])*Et1j(1))+(exp(alpha011*w[2])*Et1j(2)))
    alpha020 = (4*K[k])/((exp(alpha021*w[1])*Et2j(1))+(exp(alpha021*w[2])*Et2j(2)))
    ite=ite+1
    alpha11 = matrix(c(alpha010, alpha011, alpha020, alpha021), ncol = 4, byrow = TRUE)
    dis = ((alpha0[1,1]-alpha010)^2+(alpha0[1,2]-alpha011)^2+(alpha0[2,1]-alpha020)^2+(alpha0[2,2]-alpha021)^2)^0.5
    alpha0 = matrix(c(alpha010, alpha011, alpha020, alpha021), ncol = 2, byrow = TRUE)
    ite=ite+1
    alpha11 = matrix(c(alpha010, alpha011, alpha020, alpha021), ncol = 4, byrow = TRUE)
    alphaAll = rbind(alphaAll,alpha11)
    disAll = rbind(disAll,dis)
    if (ite>=200){
      break
    }
  }
  
  error1=(alpha0[1,1]-alpha[1,1])^2
  errorAll1 = rbind(errorAll1,error1)
  error2=(alpha0[1,2]-alpha[1,2])^2
  errorAll2 = rbind(errorAll2,error2)
  error3=(alpha0[2,1]-alpha[2,1])^2
  errorAll3 = rbind(errorAll3,error3)
  error4=(alpha0[2,2]-alpha[2,2])^2
  errorAll4 = rbind(errorAll4,error4)
  
  bias1 =(alpha0[1,1]-alpha[1,1])
  biasAll1 = rbind(biasAll1,bias1)
  bias2 =(alpha0[1,2]-alpha[1,2])
  biasAll2 = rbind(biasAll2,bias2)
  bias3 =(alpha0[2,1]-alpha[2,1])
  biasAll3 = rbind(biasAll3,bias3)
  bias4 =(alpha0[2,2]-alpha[2,2])
  biasAll4 = rbind(biasAll4,bias4)
  
}

errorAll1=errorAll1[errorAll1<1,]
errorAll2=errorAll2[errorAll2<1,]
errorAll3=errorAll3[errorAll3<1,]
errorAll4=errorAll4[errorAll4<1,]
biasAll1=biasAll1[abs(biasAll1)<10,]
biasAll2=biasAll2[abs(biasAll2)<10,]
biasAll3=biasAll3[abs(biasAll3)<10,]
biasAll4=biasAll4[abs(biasAll4)<10,]


mse=matrix(c(sum(errorAll1)/1000,sum(errorAll2)/1000,sum(errorAll3)/1000,sum(errorAll4)/1000), ncol = 4, byrow = TRUE)
bias1=matrix(c(sum(biasAll1)/1000,sum(biasAll2)/1000,sum(biasAll3)/1000,sum(biasAll4)/1000), ncol = 4, byrow = TRUE)
mse
bias1

