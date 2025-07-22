# Linear Quantile Regression (LQR) in DataSHIELD <br>- A Horizontal Federated Quantile Regression Method

**Content of this repository is part of a working paper affiliate with Singapore Eye Research Institute.**

##
The R scripts were developed under R version 4.1.2 for performing linear quantile regression (LQR) using DataSHIELD.

**[LQR_in_DataSHIELD.R](LQR_in_DataSHIELD.R)** recorded the script for estimating the regression coefficients and the variance of coefficients for linear quantile regression model.

Regression coefficients were estimated by Iterative Weighted Least Squares (IWLS). [^1]<sup>,</sup>[^2] 

Variance of regression coefficients were estimated by Powell’s kernel estimator. [^3]<sup>,</sup>[^4]

##
The mathematical details of this horizontal federated quantile regression can refer to **[TechnicalNote.md](TechnicalNote.md)**

##
**Reference:**

[^1]: Schnabel, S. K., & Eilers, P. H. C. (2013). Simultaneous estimation of quantile curves using quantile sheets. AStA Advances in Statistical Analysis, 97(1), 77–87. https://doi.org/10.1007/s10182-012-0198-1

[^2]: Waltrup, L. S., Sobotka, F., Kneib, T., & Kauermann, G. (2015). Expectile and quantile regression—David and Goliath? Statistical Modelling, 15(5), 433–456. https://doi.org/10.1177/1471082X14561155 

[^3]: Koenker, R. (2005). Quantile Regression. Cambridge University Press. https://doi.org/10.1017/CBO9780511754098 

[^4]: Kato, K. (2012). Asymptotic normality of Powell’s kernel estimator. Annals of the Institute of Statistical Mathematics, 64(2), 255–273. https://doi.org/10.1007/s10463-010-0310-9 
