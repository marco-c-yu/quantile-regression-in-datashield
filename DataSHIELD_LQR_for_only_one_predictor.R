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
  
  # indices of Intercept column in DM
  xcs = 5
  
  # the inital estimate of the quantile
  # a reasonable inital can be the unconditional tau-th quantile value
  b0 = 20
}
# initalize IRLS optimization parameters
{
  # IRLS iterative calculation will be ended IF
  # { max(|b1-b0|)<maxatol AND max(|b1-b0|/(b0+sigdigit))<maxrtol }
  # OR { total number of iteration > maxiter }
  maxatol = 1e-4 # max absolute tolerence
  maxrtol = 1e-3 # max relative tolerence
  sigdigit = 1e-8 # used together with maxrtol
  maxiter = 200 # max number of iteration
  
  Delta2=1e-8 # a small Delta^2 value to avoid division by 0 in calculating the W matrix
  verbose = FALSE
}
##########

# perform linear quantile regression from DataSHIELD for only one predictor
# intended predictor: Intercept (in this case, the solution, sol$b, 
# will be the unconditional quantile when pooling all datasets together)
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
    # conver y into matrix and vector
    ds.asMatrix(x.name = paste("DXy",colnames(dm)[yc],sep='$'),newobj = "My",datasources = connections)
    ds.matrixDimnames("My",dimnames=list(c(),colnames(dm)[yc]),newobj = "My",datasources = connections)
    ds.asNumeric("My",newobj="y",datasources=connections)
    # # create b matrix
    # ds.matrix(mdata = b0,from = "clientside.scalar",
    #           nrows.scalar = dim(b0)[1],ncols.scalar = 1,byrow = FALSE,
    #           newobj = "b",datasources = connections)
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
    # create scale b
    ds.make(b0,newobj="b",datasources=connections)
    # calc sqrt(W) as vector weighting
    {
      # ifelse(y>=X%*%b,Tau,1-Tau)/(sqrt((Ys-Xs%*%b)^2+Delta2))
      ds.Boole("y",'b',">=",newobj="w",datasources=connections)
      # ds.make("W*0.05+(1-W)*0.95",newobj="w",datasources=connections)
      ds.make("(w*tau)+(1-w)*(1-tau)",newobj="w",datasources=connections)
      # ds.make("w/( sqrt((y-Xb)^2)+1e-8 )",newobj="w",datasources=connections)
      ds.make("(w)/( sqrt(((y-b)^2)+Delta2) )",newobj="w",datasources=connections)
      ds.make("w",newobj="W",datasources=connections) # big-W: matrix squared weight
      ds.make("sqrt(w)",newobj="w",datasources=connections)# small-w: rowwise weight
    }
    # calc XWX and XWy
    {
      ds.make("DXy*w",newobj="DwXy",datasources=connections) # multiply with sqrt(w): rowwise weight
      mat.ds = ds.returnXX("DwXy",datasources=connections)
      mat.ds$XWX = mat.ds$XX[2,2]
      mat.ds$XWy = mat.ds$XX[2,1]
    }
    # solve updated beta
    b1 = lqrWLS.solve(mat.ds$XWX,mat.ds$XWy)$b
    atol = max(abs(b0-b1))
    rtol = max(abs((b0-b1)/(b0+sigdigit)))
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
##########
