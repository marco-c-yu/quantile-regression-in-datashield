# Linear Quantile Regression (LQR) in DataSHIELD <br>- A Horizontal Federated Quantile Regression Method

**Content of this repository is part of a working paper affiliate with Singapore Eye Research Institute.**

#### Developer (for mathematical theory and R program DataSHIELD implementation)

> [Marco, Chak Yan, YU](https://www.linkedin.com/in/marcocyyu/) <br>
> contact: [marcocyyu@gmail.com](mailto:marcocyyu@gmail.com) / [marcocyyu@outlook.com](mailto:marcocyyu@outlook.com) <br>

#### License

> All resources are made available under the GPL3 licence. <br>
> Please cite our paper and this github page when using any resource here. <br>

##

## Mathematical details

Mathematical details of this horizontal federated quantile regression can refer to **[TechnicalNote.md](TechnicalNote.md)**

##

## R program

The R program was developed for performing **linear quantile regression (LQR)** using **[DataSHIELD](https://datashield.org/)**.

##

### Analyst side

Tested in R version 4.1.2

**Packages required:**

> **DSI** (tested version: 1.5.0) <br>
> **dsBaseClient** (tested version: 6.2.0) <br>
> **quantreg** (tested version: 5.94) <br>

**[DataSHIELD_LQR.R](DataSHIELD_LQR.R)** recorded the script for estimating the regression coefficients and the variance of coefficients for linear quantile regression model. 
Regression coefficients were estimated by Iterative Weighted Least Squares (IWLS). [^1]<sup>,</sup>[^2] 
Variance of regression coefficients were estimated by Powell’s kernel estimator. [^3]<sup>,</sup>[^4] 

##

### Server side

DataSHIELD servers were setup using **[OBiBa/Opal docker](https://github.com/obiba/docker-opal)**.

Tested server setting was recorded in the docker-compose file: **[docker-compose-opal-4.5.yml](docker-compose-opal-4.5.yml)**,

which is based on **obiba/opal:4.5** with **datashield/rock-base:6.2-R4.1**.

##

## Reference:

[^1]: Schnabel, S. K., & Eilers, P. H. C. (2013). Simultaneous estimation of quantile curves using quantile sheets. AStA Advances in Statistical Analysis, 97(1), 77–87. https://doi.org/10.1007/s10182-012-0198-1

[^2]: Waltrup, L. S., Sobotka, F., Kneib, T., & Kauermann, G. (2015). Expectile and quantile regression—David and Goliath? Statistical Modelling, 15(5), 433–456. https://doi.org/10.1177/1471082X14561155 

[^3]: Koenker, R. (2005). Quantile Regression. Cambridge University Press. https://doi.org/10.1017/CBO9780511754098 

[^4]: Kato, K. (2012). Asymptotic normality of Powell’s kernel estimator. Annals of the Institute of Statistical Mathematics, 64(2), 255–273. https://doi.org/10.1007/s10463-010-0310-9 
