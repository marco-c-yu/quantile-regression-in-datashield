# Technical Note: Linear Quantile Regression (LQR)


### Outline

[Introduction to Quantile Regression (QR)](#introduction-to-quantile-regression-qr) 

[Linear Quantile Regression (LQR)](#linear-quantile-regression-lqr) 

[Federated LQR Algorithm](#federated-lqr-algorithm) 

[References](#references) 

##

### Introduction to Quantile Regression (QR)

**Quantile Regression (QR)** can be formulated as

$$ q=Q_{Y|X}(\tau)=inf(y:F_{Y|X}(y)\ge\tau) $$

where $q=Q_{Y|X}(\tau)$ represents the $\tau$-th conditional quantile of the response variable Y given the predictors X.

##

### Linear Quantile Regression (LQR)

For **Linear Quantile Regression (LQR)**,

$$ q=Q_{Y|X}(\tau)=inf(y:F_{Y|X}(y)\ge\tau)=X\beta_\tau $$

where $\beta_\tau$ are the regression coefficients estimated for that specific quantile $\tau$.

It can be solved by

$$ \beta_\tau = {arg\max} \\{ \sum_{y_i \ge X_i\beta_\tau} [\tau|y_i-X_i\beta_\tau|] + \sum_{y_i < X_i\beta_\tau} [(1-\tau)|y_i-X_i\beta_\tau|] \\} $$


By considering 
$|y_i-X_i\beta_\tau|=\frac{1}{\sqrt{(y_i-X_i\beta_\tau)^2}}(y_i-X_i\beta_\tau)^2$,

We have

$$ \beta_\tau = {arg\max} \\{ \sum_{y_i} [w_i(y_i-X_i\beta_\tau)^2] \\} $$

where

$$ w_i=\frac{ \tau I(y_i \ge X_i\beta_\tau ) + (1- \tau ) I(y_i < X_i\beta_\tau ) }{\sqrt{(y_i-X_i\beta_\tau)^2}} $$

It can be treated as an **Iteratively Reweighted Least Squares (IRLS)**[^1]<sup>,</sup>[^2], so that $\beta_\tau$ can be solved iteratively by

$$ \beta_\tau(t+1) = {arg\max} \\{ \sum_{y_i} [w_i(t)(y_i-X_i\beta_\tau)^2] \\} = (X^TW(t)X)^{-1}X^TW(t)y $$

where W(t) is the diagonal matrix of weights

$$ w_i(t)=\frac{ \tau I(y_i \ge X_i\beta_\tau(t) ) + (1- \tau ) I(y_i < X_i\beta_\tau(t) ) }{\sqrt{(y_i-q)^2}} $$

##

### Federated LQR Algorithm

##

### References

[^1]: Schnabel, S. K., & Eilers, P. H. C. (2013). Simultaneous estimation of quantile curves using quantile sheets. AStA Advances in Statistical Analysis, 97(1), 77–87. https://doi.org/10.1007/s10182-012-0198-1

[^2]: Waltrup, L. S., Sobotka, F., Kneib, T., & Kauermann, G. (2015). Expectile and quantile regression—David and Goliath? Statistical Modelling, 15(5), 433–456. https://doi.org/10.1177/1471082X14561155 

[^3]: Powell, J. (1991) Estimation of Monotonic Regression Models under Quantile Restrictions, in Nonparametric and Semiparametric Methods in Econometrics, W. Barnett, J. Powell, and G Tauchen (eds.), Cambridge U. Press 

[^4]: Kato, K. (2012). Asymptotic normality of Powell’s kernel estimator. Annals of the Institute of Statistical Mathematics, 64(2), 255–273. https://doi.org/10.1007/s10463-010-0310-9 
