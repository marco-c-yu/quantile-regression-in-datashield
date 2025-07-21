# Technical Note: Linear Quantile Regression (LQR)

Quantile regression can be formulated as

$$ q=Q_{Y|X}(\tau)=inf(y:F_{Y|X}(y)\ge\tau) $$

$q=Q_{Y|X}(\tau)$ represents the $\tau$-th conditional quantile of the response variable Y given the predictors X.


For linear quantile regression

$$ q=Q_{Y|X}(\tau)=inf(y:F_{Y|X}(y)\ge\tau)=X\beta_\tau $$

where $\beta_\tau$ are the regression coefficients estimated for that specific quantile $\tau$,

which can be solved by

$$ \beta_\tau = {arg\max} \\{ \sum_{y_i \ge q} [\tau|y_i-q|] + \sum_{y_i < q} [(1-\tau)|y_i-q|] \\} $$


It can also be treated as an **Iterative Weighted Least Squares (IWLS)** by considering 

$$ |y_i-q|=\frac{1}{\sqrt{(y_i-q)^2}}(y_i-q)^2 $$

so that $\beta_\tau$ can be solved iteratively by

$$ \beta_\tau(t+1) = {arg\max} \\{ \sum_{y_i} [w_i(t)(y_i-q)^2] \\} = (X^TW(t)X)^{-1}X^TW(t)y $$

where W(t) is the diagonal matrix of weights

$$ w_i(t)=\frac{ \tau I(y_i \ge X_i\beta_\tau(t) ) + (1- \tau ) I(y_i < X_i\beta_\tau(t) ) }{\sqrt{(y_i-q)^2}} $$
