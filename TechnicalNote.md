# Technical Note: Linear Quantile Regression (LQR)


### Outline

[Conditional Quantile Function](#cqf) 

[Linear Quantile Regression (LQR)](#linear-quantile-regression-lqr) 

[Iteratively Reweighted Least Squarees method for LQR](#lqr) 

[Federated LQR Algorithm](#flqr) 

[References](#references) 

##

### <a id="cfq">Conditional Quantile Function</a>

The $\tau$-th conditional quantile function of the response variable Y given the predictors X is defined as

$$ Q_{Y|X}(\tau)=inf(y:F_{Y|X}(y)\ge\tau) $$

where $F_{Y|X}$ is the cumulative distribution function (CDF) of Y given X.

##

### Linear Quantile Regression (LQR) {#lqr}

In **Linear Quantile Regression (LQR)** model, the conditional quantile function is assumed to be a linear combination of the predictors X:

$$ Q_{Y|X}(\tau)=X\beta_\tau $$

where X is a $(n \times p)$-dimensional matrix with element $x_{ij}$ represents the observed j-th predictor value of the i-th subject,

and $\beta_\tau$ is a p-dimensional column vector of regression coefficients estimated for that specific quantile $\tau$.

$Q_{Y|X}(\tau)=inf(y:F_{Y|X}(y)\ge\tau)=X\beta_\tau$ can be solved by

$$ \beta_\tau = {arg\max} \\{ \sum_{y_i \ge X_i\beta_\tau} [\tau|y_i-X_i\beta_\tau|] + \sum_{y_i < X_i\beta_\tau} [(1-\tau)|y_i-X_i\beta_\tau|] \\} $$

### Iteratively Reweighted Least Squarees method for LQR [^1]<sup>,</sup>[^2] {#irls}

By considering 

$$ |y_i-X_i\beta_\tau|=\frac{1}{\sqrt{(y_i-X_i\beta_\tau)^2}}(y_i-X_i\beta_\tau)^2 $$

we have

$$ \beta_\tau = {arg\max} \\{ \sum_{y_i} [w_i(y_i-X_i\beta_\tau)^2] \\} $$

where

$$ w_i=\frac{ \tau I(y_i \ge X_i\beta_\tau ) + (1- \tau ) I(y_i < X_i\beta_\tau ) }{\sqrt{(y_i-X_i\beta_\tau)^2}} $$

It can be considered as an **Iteratively Reweighted Least Squares (IRLS)**, so that $\beta_\tau$ can be esimated iteratively by 
$\hat\beta_\tau = \lim_{n\to\infty} \hat\beta_{\tau;t}$, 
such that

$$ \hat\beta_{\tau;t+1} = {arg\max} \\{ \sum_{y_i} [w_{i;t}( y_i - X_i \hat\beta_{\tau;t+1} )^2] \\} = (X^TW_tX)^{-1}X^TW_ty $$

where $W_t$ is the diagonal matrix of weights

$$ w_{i;t}=\frac{ \tau I(y_i \ge X_i \hat\beta_{\tau;t} ) + (1- \tau ) I(y_i < X_i \hat\beta_{\tau;t} ) }{\sqrt{(y_i - X_i \hat\beta_{\tau;t})^2}} $$

##

### Federated LQR Algorithm {#flqr}

##

### References

[^1]: Schnabel, S. K., & Eilers, P. H. C. (2013). Simultaneous estimation of quantile curves using quantile sheets. AStA Advances in Statistical Analysis, 97(1), 77–87. https://doi.org/10.1007/s10182-012-0198-1

[^2]: Waltrup, L. S., Sobotka, F., Kneib, T., & Kauermann, G. (2015). Expectile and quantile regression—David and Goliath? Statistical Modelling, 15(5), 433–456. https://doi.org/10.1177/1471082X14561155 

[^3]: Powell, J. (1991) Estimation of Monotonic Regression Models under Quantile Restrictions, in Nonparametric and Semiparametric Methods in Econometrics, W. Barnett, J. Powell, and G Tauchen (eds.), Cambridge U. Press 

[^4]: Kato, K. (2012). Asymptotic normality of Powell’s kernel estimator. Annals of the Institute of Statistical Mathematics, 64(2), 255–273. https://doi.org/10.1007/s10463-010-0310-9 
