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

# create Intercept column in DataSHIELD
{
  ds.make(toAssign=1,newobj='Intercept',datasources=connections)
  ds.cbind(c("DM",'Intercept'),newobj='DM',datasources=connections)
}
# Xy matrices from DataSHIELD
{
  # select required variable first
  ds.subset(x="DM",subset="DXy",completeCases=FALSE,cols=c(yc,xcs),datasources=connections)
  # complete cases with X and y
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
  ds.make(deltaW,newobj="deltaW",datasources=connections)
  n.ds = unlist(ds.dim("DXy",type='combine',datasources=connections))[1]
}
# lqrWLS (stepwise) from DataSHIELD
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
      # ifelse(y>=X%*%b,Tau,1-Tau)/(sqrt((Ys-Xs%*%b)^2+deltaW))
      ds.matrixMult(M1="MX",M2="b",newobj="MXb",datasources=connections)
      ds.asNumeric("MXb",newobj="Xb",datasources=connections)
      ds.Boole("y",'Xb',">=",newobj="w",datasources=connections)
      # ds.make("W*0.05+(1-W)*0.95",newobj="w",datasources=connections)
      ds.make("(w*tau)+(1-w)*(1-tau)",newobj="w",datasources=connections)
      # ds.make("w/( sqrt((y-Xb)^2)+1e-8 )",newobj="w",datasources=connections)
      ds.make("(w)/( sqrt(((y-Xb)^2)+deltaW) )",newobj="w",datasources=connections)
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
    rtol = max(abs((b0-b1)/(1+b0)))
    r = r+1
    iter <- rbind(iter,c(b1,atol,rtol))
  }
  iter <- as.data.frame(iter)
  colnames(iter) <- c(paste('b',1:dim(X)[2],sep=''),'atol','rtol')
}
# solution
{
  sol <- NULL
  sol$b <- b1
  sol$XWX <- mat.ds$XWX
  sol$XWy <- mat.ds$XWy
  sol$iter <- iter
}
# sol$b
