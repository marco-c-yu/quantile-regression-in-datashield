# Technical Note: Linear Quantile Regression (LQR)

## 
**Quantile Regression (QR)** can be formulated as

$$ q=Q_{Y|X}(\tau)=inf(y:F_{Y|X}(y)\ge\tau) $$

where $q=Q_{Y|X}(\tau)$ represents the $\tau$-th conditional quantile of the response variable Y given the predictors X.

## 
For **Linear Quantile Regression (LQR)**,

$$ q=Q_{Y|X}(\tau)=inf(y:F_{Y|X}(y)\ge\tau)=X\beta_\tau $$

where $\beta_\tau$ are the regression coefficients estimated for that specific quantile $\tau$.

It can be solved by

$$ \beta_\tau = {arg\max} \\{ \sum_{y_i \ge X_i\beta_\tau} [\tau|y_i-X_i\beta_\tau|] + \sum_{y_i < X_i\beta_\tau} [(1-\tau)|y_i-X_i\beta_\tau|] \\} $$

## 
By considering 

$|y_i-X_i\beta_\tau|=\frac{1}{\sqrt{(y_i-X_i\beta_\tau)^2}}(y_i-X_i\beta_\tau)^2$,

$$ \beta_\tau = {arg\max} \\{ \sum_{y_i} [w_i(y_i-X_i\beta_\tau)^2] \\} $$

where

$$ w_i=\frac{ \tau I(y_i \ge X_i\beta_\tau ) + (1- \tau ) I(y_i < X_i\beta_\tau ) }{\sqrt{(y_i-X_i\beta_\tau)^2}} $$

It can also be treated as an **Iterative Weighted Least Squares (IWLS)** 

so that $\beta_\tau$ can be solved iteratively by

$$ \beta_\tau(t+1) = {arg\max} \\{ \sum_{y_i} [w_i(t)(y_i-X_i\beta_\tau)^2] \\} = (X^TW(t)X)^{-1}X^TW(t)y $$

where W(t) is the diagonal matrix of weights

$$ w_i(t)=\frac{ \tau I(y_i \ge X_i\beta_\tau(t) ) + (1- \tau ) I(y_i < X_i\beta_\tau(t) ) }{\sqrt{(y_i-q)^2}} $$


