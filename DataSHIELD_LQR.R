# Developer: Marco, Chak Yan, YU
# contact: marcocyyu@gmail.com / marcocyyu@outlook.com
# License: All resources are made available under the GPL3 licence.
# Please cite our paper and 
# the github page (https://github.com/marco-c-yu/quantile-regression-in-datashield/) 
# when using this R program.

require('DSI')
require('dsBaseClient')
require('quantreg')

# DataSHIELD connection
####################
# DataSHIELD connection connects all server(s) 
# and get all project.table into DM (data for regression model)
# all tables must have the same order of columns (variables)
{
  builder <- DSI::newDSLoginBuilder()
  builder$append(server="server",url='url',user='usr',password='pwd',table='project.table',driver="OpalDriver")
  logindata <- builder$build()
  connections <- DSI::datashield.login(logins = logindata, assign = TRUE, symbol = "DM") 
}
# create Intercept column in DataSHIELD
{
  ds.make(toAssign=1,newobj='Intercept',datasources=connections)
  ds.cbind(c("DM",'Intercept'),newobj='DM',datasources=connections)
}
##########

# parameters initalization
####################
# initalize quantile regression parameters
{
  # tau: percentile to be fitted in the quantile regression
  tau = 0.2
  
  # index of column in DM to be used as the response variable y
  yc = 1
  
  # indices of columns in DM to be used as predictor variables x (including the Intercept column)
  # all variables are assumed to be continuous or binary without perfect collinearity
  # categorical variables should be coded into dummary variables with values 0/1 in DM
  # here I assumed the last column as the Intercept as we created it after loading the data DM, so that Intercept will be listed as the last column
  xcs = 2:5
  
  # the inital the column vector of regression coefficients (beta) to be used
  # a reasonable inital can be the unconditional tau-th quantile value
  # assuming the last column is the Intercept, so all previous beta elements are 0, keeping the last element to be the unconditional tau-th quantile value (assumed to be 20 in this example)
  b0 = matrix(c(0,0,0,20),ncol=1)
}
# initalize IRLS optimization parameters
{
  # IRLS iterative calculation will be ended IF
  # { max(|b1-b0|)<maxatol AND max(|b1-b0|/(b0+1e-8))<maxrtol }
  # OR { total number of iteration > maxiter }
  maxatol = 1e-4 # max absolute tolerence
  maxrtol = 1e-3 # max relative tolerence
  maxiter = 200 # max number of iteration
  
  Delta2=1e-8 # a small Delta^2 value to avoid division by 0 in calculating the W matrix
  verbose = FALSE
}
##########

# perform quantile regression from DataSHIELD
####################
# define the complete data (X,y) as DXy
# define MX,My as matrices X and y, and b as the beta vector for matrix calculation
{
  # select required variable first
  ds.subset(x="DM",subset="DXy",completeCases=FALSE,cols=c(yc,xcs),datasources=connections)
  # filter complete cases with X and y
  ds.subset(x="DXy",subset="DXy",completeCases=TRUE,datasources=connections)
  # convert data to matrix (MX,My,b) for calc
  {
    # convert X into matrix
    ds.subset(x="DXy",subset="DX",cols=(1:length(xcs))+1,datasources=connections)
    ds.asMatrix(x.name = "DX",newobj = "MX",datasources = connections)
    # conver y into matrix and vector
    ds.asMatrix(x.name = paste("DXy",colnames(dm)[yc],sep='$'),newobj = "My",datasources = connections)
    ds.matrixDimnames("My",dimnames=list(c(),colnames(dm)[yc]),newobj = "My",datasources = connections)
    ds.asNumeric("My",newobj="y",datasources=connections)
    # create b matrix
    ds.matrix(mdata = b0,from = "clientside.scalar",
              nrows.scalar = dim(b0)[1],ncols.scalar = 1,byrow = FALSE,
              newobj = "b",datasources = connections)
  }
  ds.make(tau,newobj="tau",datasources=connections)
  ds.make(Delta2,newobj="Delta2",datasources=connections)
  n.ds = unlist(ds.dim("DXy",type='combine',datasources=connections))[1]
}

# subfunction to be used inside the while loop in (1)
lqrWLS.solve <- function(XWX,XWy) {
  H = solve(XWX)
  b.wls = H %*% XWy
  b1 = b.wls
  # return
  {
    sol <- NULL
    sol$b <- b1
    sol$XWX <- XWX
    sol$XWy <- XWy
    return(sol)
  }
}
# subfunction to be used for calculating sum product statistics in (1) and (2)
ds.returnXX <- function(D,datasources=connections) {
  n = as.numeric(unlist(ds.dim(D,datasources=datasources,type="combine"))[1])
  # cov
  cov = ds.cov(x=D,datasources=datasources,type="combine")$`Variance-Covariance Matrix`
  # xbar
  {
    xcns = ds.colnames(D,datasources=datasources)[[1]]
    xbar = NULL
    for (xcn in xcns) {
      xbar = c(xbar,ds.mean(paste(D,xcn,sep="$"),datasources=datasources,type="combine")$Global.Mean[1,1])
    }
    xbar = matrix(xbar,ncol=1)
  }
  # (n-1)*cov = t(X)Y - n Xbar t(Ybar)
  XX = (n-1)*cov + n* xbar %*% t(xbar)
  # return
  {
    out = NULL
    out$XX = XX
    out$meanX = t(xbar)
    out$cov = cov
    out$n = n
  }
  return(out)
}

# (1) estimating beta in linear quantile regression using Iterative Reweighted Least Squares (IRLS)
# ref: Schnabel, S. K., & Eilers, P. H. C. (2013). Simultaneous estimation of quantile curves using quantile sheets. AStA Advances in Statistical Analysis, 97(1), 77-87. https://doi.org/10.1007/s10182-012-0198-1
# ref: Waltrup, L. S., Sobotka, F., Kneib, T., & Kauermann, G. (2015). Expectile and quantile regression-David and Goliath? Statistical Modelling, 15(5), 433-456. https://doi.org/10.1177/1471082X14561155
{
  b1 = b0
  iter <- NULL
  atol = Inf
  rtol = Inf
  r=0
  while (((atol>maxatol)|(rtol>maxrtol))&(r<maxiter)) {
    b0 = b1
    # create b matrix
    {
      ds.matrix(mdata = b0,from = "clientside.scalar",
                nrows.scalar = dim(b0)[1],ncols.scalar = 1,byrow = FALSE,
                newobj = "b",datasources = connections)
    }
    # calc sqrt(W) as vector weighting
    {
      # ifelse(y>=X%*%b,Tau,1-Tau)/(sqrt((Ys-Xs%*%b)^2+Delta2))
      ds.matrixMult(M1="MX",M2="b",newobj="MXb",datasources=connections)
      ds.asNumeric("MXb",newobj="Xb",datasources=connections)
      ds.Boole("y",'Xb',">=",newobj="w",datasources=connections)
      # ds.make("W*0.05+(1-W)*0.95",newobj="w",datasources=connections)
      ds.make("(w*tau)+(1-w)*(1-tau)",newobj="w",datasources=connections)
      # ds.make("w/( sqrt((y-Xb)^2)+1e-8 )",newobj="w",datasources=connections)
      ds.make("(w)/( sqrt(((y-Xb)^2)+Delta2) )",newobj="w",datasources=connections)
      ds.make("w",newobj="W",datasources=connections) # big-W: matrix squared weight
      ds.make("sqrt(w)",newobj="w",datasources=connections)# small-w: rowwise weight
    }
    # calc XWX and XWy
    {
      ds.make("DXy*w",newobj="DwXy",datasources=connections) # multiply with sqrt(w): rowwise weight
      # mat.ds = ds.returnXX("DwXy",datasources=connections)
      # mat.ds$XWy = matrix(mat.ds$XX[-1,1],ncol=1)
      # mat.ds$XWX = mat.ds$XX[-1,-1]
      ds.make("DXy*W",newobj="DWXy",datasources=connections) # multiply with W: matrix squared weight
      ds.subset(x="DwXy",subset="DwXy0",cols=c(1,(2:length(xcs))+1),datasources=connections) # taking out Intercept
      mat.ds = ds.returnXX("DwXy0",datasources=connections)
      mat.ds$XWX = matrix(0,nrow=length(xcs),ncol=length(xcs))
      mat.ds$XWX[2:length(xcs),2:length(xcs)] = matrix(mat.ds$XX[-1,-1],ncol=1)
      for (c in 1:length(xcs)) {
        mat.ds$XWX[1,c] = as.numeric(unlist(ds.mean(paste("DWXy",colnames(dm)[xcs[c]],sep='$'),datasources=connections,type='combine'))[1])*n.ds
        mat.ds$XWX[c,1] = mat.ds$XWX[1,c]
      }
      mat.ds$XWy = matrix(c(
        as.numeric(unlist(ds.mean(paste("DWXy",colnames(dm)[yc],sep='$'),datasources=connections,type='combine'))[1])*n.ds,
        mat.ds$XX[-1,1]),ncol=1)
    }
    # solve updated beta
    b1 = lqrWLS.solve(mat.ds$XWX,mat.ds$XWy)$b
    atol = max(abs(b0-b1))
    rtol = max(abs((b0-b1)/(b0+1e-8)))
    r = r+1
    iter <- rbind(iter,c(b1,atol,rtol))
  }
  iter <- as.data.frame(iter)
  colnames(iter) <- c(paste('b',1:dim(X)[2],sep=''),'atol','rtol')
}
# store solution as sol
{
  sol <- NULL
  sol$b <- b1
  sol$XWX <- mat.ds$XWX
  sol$XWy <- mat.ds$XWy
  sol$iter <- iter
}
# sol$b

# (2) estimating variance for hat_beta (sol$b)
# this is a recreation of quantreg::summary.rq(se='ker') in DataSHIELD
# View(quantreg::summary.rq)
# ref: Koenker, R. (2005). Quantile Regression (pp. 80-81, 116-150). Cambridge University Press. https://doi.org/10.1017/CBO9780511754098
# note that hat_beta/se~N(0,1)
# ref: Kato, K. (2012). Asymptotic normality of Powell's kernel estimator. Annals of the Institute of Statistical Mathematics, 64(2), 255-273. https://doi.org/10.1007/s10463-010-0310-9
{
  # n.ds = unlist(ds.dim("DXy",type='combine',datasources=connections))[1]
  n = n.ds
  # uhat <- c(y-X%*%b)
  {
    ds.matrix(mdata = sol$b,from = "clientside.scalar",
              nrows.scalar = dim(b0)[1],ncols.scalar = 1,byrow = FALSE,
              newobj = "b",datasources = connections)
    ds.matrixMult(M1="MX",M2="b",newobj="MXb",datasources=connections)
    ds.asNumeric("My",newobj="y",datasources=connections)
    ds.asNumeric("MXb",newobj="Xb",datasources=connections)
    ds.make("y-Xb",newobj="uhat",datasources=connections)
  }
  h <- as.numeric(quantreg::bandwidth.rq(tau, n, TRUE))
  while ((tau - h < 0) || (tau + h > 1)) h <- h/2
  # h <- (qnorm(tau + h) - qnorm(tau - h)) * min(
  #   sqrt(var(uhat)), 
  #   (quantile(uhat, 0.75) - quantile(uhat, 0.25))/1.34)
  h = (qnorm(tau+h) - qnorm(tau-h)) *
    min(
      sqrt(as.numeric(unlist(ds.var("uhat",type='combine',datasources=connections))[1])),
      diff(unlist(ds.quantileMean("uhat",type='combine',datasources=connections))[c(3,5)])/1.34)
  # D = t(X)%*%(X)/n
  D = ds.returnXX("DX",datasources=connections)$XX/n
  # kX = diag(sqrt( dnorm(uhat/h)/h )) %*% X
  {
    ds.make(h,newobj="h",datasources=connections)
    # ds.make(pi,newobj="pi",datasources=connections)
    ds.make("uhat/h",newobj="k",datasources=connections)
    ds.make("-1*(k^2)/2",newobj="k",datasources=connections)
    ds.exp("k",newobj="k",datasources=connections)
    ds.make("k/sqrt(2*pi)",newobj="k",datasources=connections)
    ds.make("k/h",newobj="k",datasources=connections)
    ds.sqrt("k",newobj="k",datasources=connections)
    ds.make("DX*k",newobj="DkX",datasources=connections)
  }
  # Hinv = solve(t(kX)%*%(kX)/n)
  H = ds.returnXX("DkX")$XX/n
  Hinv = solve(H)
  # Vb = (tau)*(1-tau)*Hinv%*%D%*%Hinv/n
  Vb = (tau)*(1-tau)*Hinv%*%D%*%Hinv/n
}
# store the variance of coefficient into sol
sol$Vb = Vb
# create summary table in sol$table
{
  sol$table = NULL
  sol$table$b = sol$b # hat_beta, the IRLS estimator of regression coefficients
  sol$table$se = sqrt(diag(sol$Vb)) # Powell's kernel estimator of the standard error of hat_beta
  sol$table = as.data.frame(sol$table)
  rownames(sol$table) = paste('column',xcs)
  sol$table$b.lcl = sol$table$b + qnorm(0.025)*sol$table$se # 95% lower confidence limit of beta
  sol$table$b.ucl = sol$table$b + qnorm(1-0.025)*sol$table$se # 95% upper confidence limit of beta
  sol$table$p = pnorm(-abs(sol$table$b/sol$table$se))*2 # two-sided Wald test p-value for beta = 0
}
# sol$table
##########
