# Technical Note: Linear Quantile Regression (LQR)

Quantile regression can be formulated as

$$q=Q_{Y|X}(\tau)=inf(y:F_{Y|X}(y)\ge\tau)$$

For linear quantile regression

$$Q_{Y|X}(\tau)=inf(y:F_{Y|X}(y)\ge\tau)=X\beta_\tau$$

which can be solved by

$$\beta_\tau = {arg\max} \\{ \sum_{y_i \ge q} [\tau|y_i-q|] + \sum_{y_i < q} [(1-\tau)|y_i-q|] \\} $$
