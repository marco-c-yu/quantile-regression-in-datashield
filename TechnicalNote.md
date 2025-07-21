# Technical Note: Linear Quantile Regression (LQR)

Quantile regression can be formulated as

$$ q=Q_{Y|X}(\tau)=inf(y:F_{Y|X}(y)\ge\tau) $$

For linear quantile regression
$$ q=Q_{Y|X}(\tau)=inf(y:F_{Y|X}(y)\ge\tau)=X\beta_\tau $$

which can be solved by
$$ \beta_\tau = {arg\max} \\{ \sum_{y_i \ge q} [\tau|y_i-q|] + \sum_{y_i < q} [(1-\tau)|y_i-q|] \\} $$

It can also be treated as **Iterative Weighted Least Squares (IWLS)** by considering 
$$ |y_i-q|=\frac{1}{\sqrt{(y_i-q)^2}}(y_i-q)^2 $$
which gives weights
