# SpatialWavePredict Toolbox
## 📌 Overview

## Overview

**SpatialWavePredict Toolbox** is a **user-friendly MATLAB-based toolkit** designed for fitting and forecasting epidemic trajectories using the **ensemble spatial wave sub-epidemic framework**. This framework captures diverse wave dynamics by aggregating multiple asynchronous growth processes and has demonstrated superior performance in short-term forecasts of various infectious disease outbreaks, including SARS, Ebola, and early waves of the COVID-19 pandemic in the US.

Key functionalities include:

- **Fitting models to time series data using a five-parameter epidemic wave model that aggregates overlapping sub-epidemics.**
- **Estimating parameters with quantified uncertainty through parametric bootstrapping.**
- **Plotting fits and AICc values of top-ranked models.**
- **Generating forecasts and ensemble forecasts based on top-ranked models.**
- **Quantifying forecasting performance using metrics that evaluate point and distributional forecasts, including the weighted interval score.**

Additional features:

- **Support for different parameter estimation approaches (least-squares, maximum likelihood estimation).**
- **Ability to assume different error structures (normal, Poisson, negative binomial).**
- **Selection of underlying functions for the sub-epidemic building block (generalized-logistic model, Richards model, and the generalized Richards model).**
- **Option to choose between two decline functions for sub-epidemic sizes: Exponential and Power-law.**


## 📺 Tutorial Resources

- **Tutorial Paper:** [SpatialWavePredict: A tutorial-based primer and toolbox for forecasting growth trajectories using the ensemble spatial wave sub-epidemic modeling framework](https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-024-02241-2)
- **Video Tutorial:** [YouTube Series on SpatialWavePredict Toolbox](https://www.youtube.com/watch?v=qxuF_tTzcR8)

---

    
# Installation requirements

The n-subepidemic framework toolbox requires a MATLAB installation.

# Fitting the model to your data

To use the toolbox to fit the spatial wave sub-epidemic framework to your data, you just need to:

<ul>
    <li>download the code </li>
    <li>create 'input' folder in your working directory where your data is located </li>
    <li>create 'output' folder in your working directory where the output files will be stored</li>      
    <li>open a MATLAB session </li>
    <li>define the model parameter values and time series parameters by editing <code>options.m</code> </li>
    <li>run the function <code>Run_SW_subepidemicFramework.m</code> </li>
</ul>
  
# Plotting the fits of the top-ranked models and parameter estimates

After fitting the model to your data, you can use the toolbox to plot the model fits and parameter estimates as follows:

<ul>
    <li>define the model parameter values and time series parameters by editing <code>options.m</code></li>
    <li>run the function <code>plotFit_SW_subepidemicFramework.m</code> </li>
</ul>
    
# Plotting the top-ranked subepidemic model profiles and the corresponding AIC values

After fitting the model to your data, you can use the toolbox to plot the subepidemic profiles and AICc values as follows:

<ul>
    <li>define the model parameter values and time series parameters by editing <code>options.m</code></li>
    <li>run the function <code>plotRankings_SW_subepidemicFramework.m</code></li>
</ul>
    
# Generating and plotting forecasts of the top-ranked and ensemble subepidemic models

After fitting the model to your data, you can use the toolbox to plot forecasts derived from the top-ranked and ensemble subepidemic models as follows:

<ul>
    <li>define the model parameter values and time series parameters by editing <code>options.m</code> and <code>options_forecast.m</code></li>
    <li>run the function <code>plotForecast_subepidemicFramework.m</code></li>
</ul>
    
# Publications

<ul>
    
 <li>Chowell, G., Tariq, A., & Hyman, J. M. (2019). A novel sub-epidemic modeling framework for short-term forecasting epidemic waves. BMC medicine, 17(1), 1-18.  </li>

 <li>Chowell, G., Rothenberg, R., Roosa, K., Tariq, A., Hyman, J. M., & Luo, R. (2022). Sub-epidemic Model Forecasts During the First Wave of the COVID-19 Pandemic in the USA and European Hotspots. In Mathematics of Public Health (pp. 85-137). Springer, Cham. </li>

</ul>

# Disclaimer

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3 of the License.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.  
