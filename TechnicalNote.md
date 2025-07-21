# Technical Note: Linear Quantile Regression (LQR)

Quantile regression can be formulated as
<br>
$$ q=Q_{Y|X}(\tau)=inf(y:F_{Y|X}(y)\ge\tau) $$
<br>
$q=Q_{Y|X}(\tau)$ represents the $\tau$-th conditional quantile of the response variable Y given the predictors X.
<br>
<br>
For linear quantile regression
<br>
$$ q=Q_{Y|X}(\tau)=inf(y:F_{Y|X}(y)\ge\tau)=X\beta_\tau $$
<br>
where $\beta_\tau$ are the regression coefficients estimated for that specific quantile $\tau$,
<br>
which can be solved by
<br>
$$ \beta_\tau = {arg\max} \\{ \sum_{y_i \ge q} [\tau|y_i-q|] + \sum_{y_i < q} [(1-\tau)|y_i-q|] \\} $$
<br>
<br>
It can also be treated as an **Iterative Weighted Least Squares (IWLS)** by considering 
<br>
$$ |y_i-q|=\frac{1}{\sqrt{(y_i-q)^2}}(y_i-q)^2 $$
<br>
so that $\beta_\tau$ can be solved iteratively by
<br>
$$ \beta_\tau(t+1) = {arg\max} \\{ \sum_{y_i} [w_i(t)(y_i-q)^2] \\} = (X^TW(t)X)^{-1}X^TW(t)y $$
<br>
where W(t) is the diagonal matrix of weights
<br>
$$ w_i(t)=\frac{ \tau I(y_i \ge X_i\beta_\tau(t) ) + (1- \tau ) I(y_i < X_i\beta_\tau(t) ) }{\sqrt{(y_i-q)^2}} $$
<br>
