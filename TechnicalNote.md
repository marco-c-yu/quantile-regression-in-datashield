# Technical Note: Horizontal Federated Linear Quantile Regression (LQR)


### Contents

[1.1 Conditional Quantile Function](#11-conditional-quantile-function) 

[1.2 Linear Quantile Regression (LQR)](#12-linear-quantile-regression-lqr) 

[1.3 Iteratively Reweighted Least Squares (IRLS) method for LQR](#13-iteratively-reweighted-least-squares-irls-method-for-lqr) 

[1.4 Asymptotic normality of the regression coefficients estimator of LQR](#14-asymptotic-normality-of-the-regression-coefficients-estimator-of-lqr)

[2.1 Horizontal Federated Learning](#21-horizontal-federated-learning) 

[2.2 Horizontal Federated LQR Algorithm](#22-horizontal-federated-lqr-algorithm) 

[References](#references) 

##

### 1.1 Conditional Quantile Function

The $\tau$-th conditional quantile function of the response variable Y given the predictors X is defined as

$$ Q_{Y|X}(\tau)=inf(y:F_{Y|X}(y)\ge\tau) $$

where $F_{Y|X}$ is the cumulative distribution function (CDF) of Y given X.

##

### 1.2 Linear Quantile Regression (LQR)

In **Linear Quantile Regression (LQR)**[^1]<sup>,</sup>[^2] model, the conditional quantile function is assumed to be a linear combination of the predictors X:

$$ Q_{Y|X}(\tau)=X\beta_\tau $$

where X is a $(n \times p)$-dimensional matrix with element $x_{ij}$ representing the observed j-th predictor value of the i-th subject,

y is a $(n \times 1)$-dimentional column vector with element $y_i$ representing the observed response value of the i-th subject,

and $\beta_\tau$ is a $(p \times 1)$-dimensional column vector of regression coefficients estimated for that specific quantile $\tau$.

$Q_{Y|X}(\tau)=inf(y:F_{Y|X}(y)\ge\tau)=X\beta_\tau$ can be solved by

$$ \beta_\tau = {arg\max} \\{ \sum_{y_i \ge X_i\beta_\tau} [\tau|y_i-X_i\beta_\tau|] + \sum_{y_i < X_i\beta_\tau} [(1-\tau)|y_i-X_i\beta_\tau|] \\} $$

where $X_i=[x_{i1},...,x_{ip}]$ is a $(1 \times p)$-dimensional row vector of the predictor values for the i-th subject.

### 1.3 Iteratively Reweighted Least Squares (IRLS) method for LQR

By considering 

$$ |y_i-X_i\beta_\tau|=\frac{1}{\sqrt{(y_i-X_i\beta_\tau)^2}}(y_i-X_i\beta_\tau)^2 $$

we have

$$ \beta_\tau = {arg\max} \\{ \sum_{y_i} [w_i(y_i-X_i\beta_\tau)^2] \\} $$

where

$$ w_i=\frac{ \tau I(y_i \ge X_i\beta_\tau ) + (1- \tau ) I(y_i < X_i\beta_\tau ) }{\sqrt{(y_i-X_i\beta_\tau)^2}} $$

It can be considered as an **Iteratively Reweighted Least Squares (IRLS)**[^3], 

so that $\beta_\tau$ can be esimated iteratively by 
$\hat\beta_\tau = \lim_{n\to\infty} \hat\beta_{\tau;t}$, 
such that

$$ \hat\beta_{\tau;t+1} = {arg\max} \\{ \sum_{y_i} [w_{i;t}( y_i - X_i \hat\beta_{\tau;t+1} )^2] \\} = (X^TW_tX)^{-1}X^TW_ty $$

where $W_t$ is the diagonal matrix of weights

$$ w_{i;t}=\frac{ \tau I(y_i \ge X_i \hat\beta_{\tau;t} ) + (1- \tau ) I(y_i < X_i \hat\beta_{\tau;t} ) }{\sqrt{(y_i - X_i \hat\beta_{\tau;t})^2}} $$

In practice, it is suggested to approximate the weights by adding a small value $\Delta^2$ in the denominator, such that 

$$ w_{i;t}=\frac{ \tau I(y_i \ge X_i \hat\beta_{\tau;t} ) + (1- \tau ) I(y_i < X_i \hat\beta_{\tau;t} ) }{\sqrt{(y_i - X_i \hat\beta_{\tau;t})^2 + \Delta^2}} $$

for numerical stability.

### 1.4 Asymptotic normality of the regression coefficients estimator of LQR

Under suitable regularity conditions[^4], $\hat\beta_\tau$ is a $\sqrt{n}$-consistent estimator of $\beta_\tau$ and 

$$\sqrt{n}(\hat\beta_\tau - \beta_\tau) \rightarrow N(0,\Sigma_{\beta,\tau})$$

where $\Sigma_{\beta,\tau} = J_\tau^{-1} \Sigma_\tau J_\tau^{-1}$, 
$J_\tau = E(f(X\beta_\tau|X) X^TX)$ and 
$\Sigma_\tau = \tau(1-\tau) E(X^TX)$

$J_\tau$ can be estimated by the Powell's kernel estimator

$\hat{J}_\tau = \frac{1}{nh} \sum_{i}^{n} K(\frac{y_i - X_i \hat\beta_\tau}{h}) X^TX$

[^5]<sup>,</sup>[^6]

##

### 2.1 Horizontal Federated Learning

A **federated learning system** is a learning process trained on data sets distributed across multiple parties while preventing data leakage.

**Horizontal Federated Learning (a.k.a. sample-based federated learning)** refers to the scenarios that all data sets share the same feature space but with different sample space.[^7]

There are other categories of federated learning, including **Vertical Federated Learning** and **Federated Transfer Learning**, which will not be discussed here.

##

### 2.2 Horizontal Federated LQR Algorithm

In [1.3](#13-iteratively-reweighted-least-squares-irls-method-for-lqr), we showed that the regression coefficients of LQR can be solved by IRLS,

it is similar to the algorithm for solving federated generalized linear model (GLM).[^8]

The proposed algorithms for regression coefficients and variance of coefficients estimation of horizontal federated LQR are listed as follow:

> #### Algorithm 1: regression coefficients estimation of horizontal federated LQR <br>
> 1. initialize a global inital estimation, $\beta_{\tau;0}$, in server
> 2. in each party node, compute
> 3. $w_{i;t}=\frac{ \tau I(y_i \ge X_i \hat\beta_{\tau;t} ) + (1- \tau ) I(y_i < X_i \hat\beta_{\tau;t} ) }{\sqrt{(y_i - X_i \hat\beta_{\tau;t})^2 + \Delta^2}}$
> 4. 



##

### References

[^1]: Furno, M., & Vistocco, D. (2018). Quantile regression: Estimation and simulation. Wiley. https://doi.org/10.1002/9781118863718

[^2]: Waltrup, L. S., Sobotka, F., Kneib, T., & Kauermann, G. (2015). Expectile and quantile regression—David and Goliath? Statistical Modelling, 15(5), 433–456. https://doi.org/10.1177/1471082X14561155 

[^3]: Schnabel, S. K., & Eilers, P. H. C. (2013). Simultaneous estimation of quantile curves using quantile sheets. AStA Advances in Statistical Analysis, 97(1), 77–87. https://doi.org/10.1007/s10182-012-0198-1

[^4]: Koenker, R. (2005). Quantile Regression (pp. 116-150). Cambridge University Press. https://doi.org/10.1017/CBO9780511754098

[^5]: Powell, J. (1991) Estimation of Monotonic Regression Models under Quantile Restrictions, in Nonparametric and Semiparametric Methods in Econometrics, W. Barnett, J. Powell, and G Tauchen (eds.), Cambridge U. Press 

[^6]: Kato, K. (2012). Asymptotic normality of Powell’s kernel estimator. Annals of the Institute of Statistical Mathematics, 64(2), 255–273. https://doi.org/10.1007/s10463-010-0310-9 

[^7]: Yang, Q., Liu, Y., Chen, T., & Tong, Y. (2019). Federated Machine Learning: Concept and Applications. ACM Transactions on Intelligent Systems and Technology, 10(2), 19. https://doi.org/10.1145/3298981

[^8]: Cellamare, M., van Gestel, A. J., Alradhi, H., Martin, F., & Moncada-Torres, A. (2022). A Federated Generalized Linear Model for Privacy-Preserving Analysis. Algorithms 2022, Vol. 15, Page 243, 15(7), 243. https://doi.org/10.3390/A15070243
