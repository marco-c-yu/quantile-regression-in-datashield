# Technical Note: <br> Horizontal Federated Linear Quantile Regression (LQR) <br> based on Iteratively Reweighted Least Squares (IRLS) method

**Content of this repository is part of a working paper from the [Asian Eye Epidemiology Consortium (AEEC)](https://www.snec.com.sg/research-innovation/research-groups-platforms/research-groups/ocular-epidemiology), <br>affiliate with [Singapore Eye Research Institute](https://www.snec.com.sg/research-innovation/about-seri).**

#### Developer (for mathematical theory and R program DataSHIELD implementation)

> [Marco, Chak Yan, YU](https://www.linkedin.com/in/marcocyyu/) <br>
> contact: [marcocyyu@gmail.com](mailto:marcocyyu@gmail.com) / [marcocyyu@outlook.com](mailto:marcocyyu@outlook.com) <br>

#### License

> All resources are made available under the GPL3 licence. <br>
> Please cite our paper and this github page when using any resource here. <br>

##

### Contents

[1.1 Conditional Quantile Function](#11-conditional-quantile-function) 

[1.2 Linear Quantile Regression (LQR)](#12-linear-quantile-regression-lqr) 

[1.3 Iteratively Reweighted Least Squares (IRLS) method for LQR](#13-iteratively-reweighted-least-squares-irls-method-for-lqr) 

[1.4 Asymptotic normality of the regression coefficients estimator of LQR](#14-asymptotic-normality-of-the-regression-coefficients-estimator-of-lqr)

[2.1 Horizontal Federated Learning](#21-horizontal-federated-learning) 

[2.2 Horizontal Federated LQR Algorithm](#22-horizontal-federated-lqr-algorithms) 

[3.1 Extension of LQR: Allowing for nonlinearity](#31-extension-of-lqr-allowing-for-nonlinearity)

[3.2 Extension of LQR: Multiple non-crossing quantile estimation](#32-extension-of-lqr-multiple-non-crossing-quantile-estimation)

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

It can be considered as an **Iteratively Reweighted Least Squares (IRLS)**[^3]<sup>,</sup>[^4], 

so that $\beta_\tau$ can be esimated iteratively by 
$\hat\beta_\tau = \lim_{t\to\infty} \hat\beta_{\tau;t}$, 
such that

$$ \hat\beta_{\tau;t+1} = {arg\max} \\{ \sum_{y_i} [w_{i;t}( y_i - X_i \hat\beta_{\tau;t+1} )^2] \\} = (X^TW_tX)^{-1}X^TW_ty $$

where $W_t$ is the diagonal matrix of weights

$$ w_{i;t}=\frac{ \tau I(y_i \ge X_i \hat\beta_{\tau;t} ) + (1- \tau ) I(y_i < X_i \hat\beta_{\tau;t} ) }{\sqrt{(y_i - X_i \hat\beta_{\tau;t})^2}} $$

In practice, it is suggested to approximate the weights by adding a small value $\Delta^2$ in the denominator, such that 

$$ w_{i;t}=\frac{ \tau I(y_i \ge X_i \hat\beta_{\tau;t} ) + (1- \tau ) I(y_i < X_i \hat\beta_{\tau;t} ) }{\sqrt{(y_i - X_i \hat\beta_{\tau;t})^2 + \Delta^2}} $$

for numerical stability.[^4]

### 1.4 Asymptotic normality of the regression coefficients estimator of LQR

Under suitable **regularity conditions**[^5], $\hat\beta_\tau$ is a $\sqrt{n}$-consistent estimator of $\beta_\tau$ and 

$$\sqrt{n}(\hat\beta_\tau - \beta_\tau) \rightarrow N(0,\Sigma_{\beta,\tau})$$

where $\Sigma_{\beta,\tau} = J_\tau^{-1} \Sigma_\tau J_\tau^{-1}$, 
$J_\tau = E(f(X\beta_\tau|X) X^TX)$ and 
$\Sigma_\tau = \tau(1-\tau) E(X^TX)$

$J_\tau$ can be estimated by the **Powell's kernel (PK) estimator**

$$ \hat J_\tau = \frac{1}{nh} \sum_{i=1}^{n} K(\frac{y_i - X_i \hat\beta_\tau}{h}) X_i^T X_i $$

where $K(\cdot)$ is a kernel density function.[^6]<sup>,</sup>[^7]

The $(1-\alpha) \times 100\\%$ confidence interval of $\beta_{\tau,j}$ for the j-th predictor variable can be estimated by

$$\hat\beta_{\tau,j} \pm \Phi^{-1}(1-\alpha/2) (\sigma_{\beta,\tau,j}^2 /n)^{1/2} $$

where $\hat\beta_{\tau,j}$ is the j-th element of $\hat\beta_\tau$ and $\sigma_{\beta,\tau,j}^2$ is the $(j,j)$ element of $\Sigma_{\beta,\tau}$.

##

### 2.1 Horizontal Federated Learning

A **federated learning system** is a learning process trained on data sets distributed across multiple parties while preventing data leakage.

**Horizontal Federated Learning (a.k.a. sample-based federated learning)** refers to the scenarios that all data sets share the same feature space but with different sample space.[^8]

There are other categories of federated learning, including **Vertical Federated Learning** and **Federated Transfer Learning**, which will not be discussed here.

### 2.2 Horizontal Federated LQR Algorithms

In [1.3](#13-iteratively-reweighted-least-squares-irls-method-for-lqr), we showed that the regression coefficients of LQR can be solved by IRLS,

it is similar to the algorithm for solving federated generalized linear model (GLM).[^9]

Estimation for variance of coefficient follows the implementation of Powell's kernel estimator in the **quantreg**[^10] package in **R**.

The proposed algorithms for regression coefficients and variance of coefficients estimation of horizontal federated LQR are summarized as follow:

> #### Algorithm 1: IRLS estimator of the regression coefficients of horizontal federated LQR[^4]<sup>,</sup>[^9] <br>
> In server, <br>
> 1: initialize an inital estimatior, $\beta_{\tau;0}$ <br>
> <br>
> In each party node, m, for $(t \ge 0)$, <br>
> 2: compute $w_{m,i;t}=\frac{ \tau I(y_m,i \ge X_m,i \hat\beta_{\tau;t} ) + (1- \tau ) I(y_m,i < X_m,i \hat\beta_{\tau;t} ) }{\sqrt{(y_m,i - X_m,i \hat\beta_{\tau;t})^2 + \Delta^2}}$ <br>
> 3: compute the weighted predictors $(W_{m,t}^{1/2}X_m)$ and weighted response $(W_{m,t}^{1/2}y_m)$ <br>
> 4: compute $(X_m^TW_{m,t}X_m)$ and $(X_m^TW_{m,t}y_m)$ <br>
> <br>
> In server, <br>
> 5: comppute $(X^TW_tX) = \sum_m (X_m^TW_{m,t}X_m)$ and $(X^TW_ty) = \sum_m (X_m^TW_{m,t}y_m)$ <br>
> 6: solve the matrix inverse $(X^TW_tX)^{-1}$ <br>
> 7: compute $\beta_{\tau;t+1} = (X^TW_tX)^{-1}X^TW_ty$ <br>
> 8: replace $t$ by $t+1$ <br>
> <br>
> 9: Repeat steps 2-8 until $|\beta_{\tau;t}-\beta_{\tau;t-1}|<\delta$, for some pre-specified small tolerence value $\delta$. <br>
> 10: After the end loop of step 9, $\hat\beta_\tau=\beta_{\tau;t}$ is the IRLS estimator of the regression coefficients of LQR. <br>
> <br>
> ###### Remarks:
> * $y_m$ and $X_m$ represents the observed response and predictors in the m-th party. <br>
> * $(X_m^TW_{m,t}X_m)$ and $(X_m^TW_{m,t}y_m)$ can be derived from the mean and covariance matrix of $(W_{m,t}^{1/2}X_m)$ and $(W_{m,t}^{1/2}y_m)$ using the property $Cov(X,Y)=E(XY)-E(X)E(Y)$ <br>
> * This federated LQR IRLS estimator will give the same estimate as the pooled LQR with access to all Individual Participant Data (IPD). <br>

> #### Algorithm 2: Powell's kernel (PK) estimator of the variance of coefficients of horizontal federated LQR[^6]<sup>,</sup>[^7] <br>
> 1: Obtain the IRLS estimator of regression coefficients, $\hat\beta_\tau$, in Algorithm 1. <br>
> <br>
> In server, <br>
> 2: compute $h = n^{-1/3} \times \Phi^{-1}(1 - \alpha/2)^{2/3} \times \\{ \frac{ 1.5 \times [\phi(\Phi^{-1}(\tau))]^2 }{ 2 \times [\Phi^{-1}(\tau)]^2 + 1 } \\} ^{1/3}$, where $\alpha=0.05$. <br>
> 3: while $(\tau-h <0)$ OR $(\tau+h >1)$, replace $h$ by $h/2$ until both while-conditions in this step are false. <br>
> <br>
> In each party node, m, <br>
> 4: compute $u_m=(y_m - X_m \hat\beta_\tau)$ <br>
> 5: compute the sample mean $\bar u_m$ and variance $v(u_m)$ of $u_m$ <br>
> 6: compute the sum of squares of $u_m$: $\sum_{i} u_{m,i}^2 = (n_m-1) \times v(u_m) + n_m \times \bar u_m^2$ <br>
> <br>
> In server, <br>
> 7: compute $v(u)=\frac{1}{n-1} [\sum_{} u^2 - (\sum_{} u)^2/n]=\frac{1}{n-1} [\sum_{m} \sum_{i} u_{m,i}^2 - (\sum_{m} n_m \times \bar u_m)^2/n]$ <br>
> 8: replace $h$ by $h=[\Phi^{-1}(\tau+h)-\Phi^{-1}(\tau-h)] \times min(\sqrt{v(u)},(Q_u(0.75)-Q_u(0.25))/1.34)$ <br>
> <br>
> In each party node, m, <br>
> 9: compute $X_m^T X_m$ <br>
> 10: compute the diagonal matrix $k_m=diag(\sqrt{\phi(u_m/h)/h})$ <br>
> 11: compute the weighted predictor matrix $k_m X_m$ <br>
> 12: calculate $H_m=(k_m X_m)^T (k_m X_m)$ <br>
> <br>
> In server, <br>
> 13: compute $H=\sum_m H_m$ and solve the matrix inverse $H^{-1}$ <br>
> 14: compute $X^T X = \sum_m X_m^T X_m$ <br>
> 15: estimate the variance of coefficients by $\hat\Sigma_{\beta,\tau}/n = \tau(1-\tau) H^{-1} (X^T X) H^{-1}/n$ <br>
> <br>
> ###### Remarks:
> This algorithm implemented the default selection of kernel density function $K(\cdot)$ and bandwidth $h$ implemented in **summary.rq** and **bandwidth.rq** functions in the **quantreg**[^10] package in **R**. <br>
> * $\Phi(\cdot)$ is the cumulative distribution function (CDF) of a standard normal distribution and <br>
> * $\phi(\cdot)$ is the probability density function (PDF) of a standard normal distribution.<br>
> * $Q_u(\cdot)$ is the quantile function of u. <br>
> * This federated LQR PK estimator will give the same variance of coefficients estimate as the pooled LQR PK estimator with access to all Individual Participant Data (IPD). <br>

##

### 3.1 Extension of LQR: Allowing for nonlinearity

Nonlinearity can be modelled by extending the LQR using **Restricted Cubic Splines (RCS)** functions.[^11]

For example a $K$ knots $(t_1<t_2<...<t_K)$ RCS quantile regression can be formulated as

$$ Q_{Y|x}(\tau)=\beta_{\tau,0}+\sum_{j=1}^{K-1} \beta_{\tau,j} \times s_j $$

where $s_1 = x$ and

$s_j = (x-t_{j-1})^3I(x-t_{j-1}>0) - (x-t_{K-1})^3I(x-t_{K-1}>0)\frac{(t_{K}-t_{j-1}>0)}{(t_{K}-t_{K-1}>0)} + (x-t_{K})^3I(x-t_{K}>0)\frac{(t_{K-1}-t_{j-1}>0)}{(t_{K}-t_{K-1}>0)}$ for $j=2,...,K-1$

By storing $s_1$ to $s_{K-1}$ as variables in the datasets, it reduced to a LQR problem of fitting the quantile of Y on $s_1$, $s_2$,..., $s_{K-1}$.

### 3.2 Extension of LQR: Multiple non-crossing quantile estimation

In some situations, we might want to estimate multiple quantiles simultaneously and impose the **non-crossing constraints**: 
$Q_{Y|x}(\tau_i) < Q_{Y|x}(\tau_j) \iff \tau_i < \tau_j$

In LQR, the non-crossing constraints are equivalent to $X\beta_{\tau_i} < X\beta_{\tau_j} \iff \tau_i < \tau_j$,

which can be solved by **Inequality Constrained Least-Squares (ICLS)**[^12]<sup>,</sup>[^13] in combining with **Iteratively Reweighted Least Squares (IRLS)**[^3]<sup>,</sup>[^4]:

> In each iteration, <br>
> 1: Solve the unconstrainted regression coefficients by **Weighted Least Squares (WLS)**: $\hat\beta=(X^TW_tX)^{-1}X^TW_ty$ <br>
> 2: Estimate $\hat\beta_c$ with constraints $(A_1 \hat\beta_c \gg c_1)$ and $(A_2 \hat\beta_c = c_2)$ by **ICLS**: 
> $\hat\beta_c = \hat\beta + (Z^T Z)^{-1} A_2^T(A_2(Z^T Z)^{-1}A_2^T)^{-1} (c_2-A_2\hat\beta)$, 
> where $Z=W^{1/2}X$. <br>
> 3: Update the diagonal elements of the diagonal matrix $W$ by
> $w_i=[ \tau I(y_i \ge X_i \hat\beta_c ) + (1- \tau ) I(y_i < X_i \hat\beta_c ) ] / [(y_i - X_i \hat\beta_c)^2 + \Delta^2]^{1/2}$. <br>
> Repeat steps 1-3 until $\hat\beta_c$ converged, and $\hat\beta_c$ will be the constrained regression coefficient estimator. <br>

By combining the horizontal federated LQR IRLS algorithm together with the ICLS method, simultaneous non-crossing LQR can also be solved iteratively by matrix calculation in horizontal federated learning.

##

### References

[^1]: Koenker, R., & Bassett, G. (1978). Regression Quantiles. Econometrica, 46(1), 33. https://doi.org/10.2307/1913643

[^2]: Furno, M., & Vistocco, D. (2018). Quantile regression: Estimation and simulation. Wiley. https://doi.org/10.1002/9781118863718

[^3]: Schlossmacher, E. J. (1973). An iterative technique for absolute deviations curve fitting. Journal of the American Statistical Association, 68(344), 857–859. https://doi.org/10.1080/01621459.1973.10481436;PAGE:STRING:ARTICLE/CHAPTER

[^4]: Schnabel, S. K., & Eilers, P. H. C. (2013). Simultaneous estimation of quantile curves using quantile sheets. AStA Advances in Statistical Analysis, 97(1), 77–87. https://doi.org/10.1007/s10182-012-0198-1

[^5]: Koenker, R. (2005). Quantile Regression (pp. 116-150). Cambridge University Press. https://doi.org/10.1017/CBO9780511754098

[^6]: Koenker, R. (2005). Quantile Regression (pp. 80-81). Cambridge University Press. https://doi.org/10.1017/CBO9780511754098

[^7]: Kato, K. (2012). Asymptotic normality of Powell’s kernel estimator. Annals of the Institute of Statistical Mathematics, 64(2), 255–273. https://doi.org/10.1007/s10463-010-0310-9 

[^8]: Yang, Q., Liu, Y., Chen, T., & Tong, Y. (2019). Federated Machine Learning: Concept and Applications. ACM Transactions on Intelligent Systems and Technology, 10(2), 19. https://doi.org/10.1145/3298981

[^9]: Cellamare, M., van Gestel, A. J., Alradhi, H., Martin, F., & Moncada-Torres, A. (2022). A Federated Generalized Linear Model for Privacy-Preserving Analysis. Algorithms 2022, Vol. 15, Page 243, 15(7), 243. https://doi.org/10.3390/A15070243

[^10]: Roger Koenker (2022). quantreg: Quantile Regression. R package version 5.94. https://CRAN.R-project.org/package=quantreg

[^11]: Marrie, R. A., Dawson, N. V., & Garland, A. (2009). Quantile regression and restricted cubic splines are useful for exploring relationships between continuous variables. Journal of Clinical Epidemiology, 62(5), 511-517.e1. https://doi.org/10.1016/J.JCLINEPI.2008.05.015

[^12]: Liew, C. K. (1976). Inequality Constrained Least-Squares Estimation. Journal of the American Statistical Association, 71(355), 746. https://doi.org/10.2307/2285614

[^13]: Goldfarb, D., & Idnani, A. (1983). A numerically stable dual method for solving strictly convex quadratic programs. Mathematical Programming, 27(1), 1–33. https://doi.org/10.1007/BF02591962
