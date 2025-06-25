# quantile-regression-in-datashield
The R scripts were developed under R version 4.1.2 for performing linear quantile regression using DataSHIELD.


**DataSHIELD_LQR_subfunc_solveB.R** recorded the script for estimating the regression coefficients of linear quantile regression,
where regression coefficients were estimated by Iterative Weighted Least Squares (IWLS) (Schnabel and Eilers; Waltrup, et al) 

Reference:

Schnabel, S. K., & Eilers, P. H. C. (2013). Simultaneous estimation of quantile curves using quantile sheets. AStA Advances in Statistical Analysis, 97(1), 77–87. https://doi.org/10.1007/s10182-012-0198-1

Waltrup, L. S., Sobotka, F., Kneib, T., & Kauermann, G. (2015). Expectile and quantile regression—David and Goliath? Statistical Modelling, 15(5), 433–456. https://doi.org/10.1177/1471082X14561155 


**DataSHIELD_LQR_subfunc_solveVb.R** recorded the script for estimating variance of regression coefficient estimators following Powell’s kernel estimator.

Reference:

Powell, J. (1991) Estimation of Monotonic Regression Models under Quantile Restrictions, in Nonparametric and Semiparametric Methods in Econometrics, W. Barnett, J. Powell, and G Tauchen (eds.), Cambridge U. Press 

Kato, K. (2012). Asymptotic normality of Powell’s kernel estimator. Annals of the Institute of Statistical Mathematics, 64(2), 255–273. https://doi.org/10.1007/s10463-010-0310-9 
