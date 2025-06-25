# pre-existing objects
# sol$b # solved
# yc
# xcs

# create Intercept column DataSHIELD
if (FALSE) {
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
    # conver y into matrix
    ds.asMatrix(x.name = paste("DXy",colnames(dm)[yc],sep='$'),newobj = "My",datasources = connections)
    ds.matrixDimnames("My",dimnames=list(c(),colnames(dm)[yc]),newobj = "My",datasources = connections)
  }
  ds.make(tau,newobj="tau",datasources=connections)
  n.ds = unlist(ds.dim("DXy",type='combine',datasources=connections))[1]
}
# lqr.Vb.Powell: variance estimate for beta (beta/se~t(n-nbeta))
{
  # View(quantreg::summary.rq)
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
sol$Vb = Vb
sol$n = n
# lqr.Vb.Powell(rbind(X,X.ds),rbind(y,y.ds),sol$b,tau,hs=TRUE)
