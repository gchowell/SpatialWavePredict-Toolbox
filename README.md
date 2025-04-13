# SpatialWavePredict Toolbox
 A Matlab toolbox for fitting and forecasting epidemic trajectories using the spatial wave sub-epidemic framework 

<p> <a href="https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-024-02241-2" target="_blank">SpatialWavePredict Tutorial</a><\p> 
<a href="https://www.youtube.com/watch?v=qxuF_tTzcR8" target="_blank">Video Tutorial</a>


<p> It carries out the following tasks: </p> 
<ul>
    <li>fitting models to time series data,</li>
    <li>estimation of parameters with quantified uncertainty,</li>
    <li>plotting the fits of the top-ranked models,</li>
    <li>plotting the AICc values of the top-ranked models,</li>
    <li>generates forecasts of the top-ranked models</li>
    <li>gnerates ensemble forecasts based on the top-ranked models.</li>
</ul>

<p> Additional features include:</p>

<ul>
    <li>fitting models using different parameter estimation approaches (least-squares, maximum likelihood estimation),</li>
    <li>fitting models using assuming different error structures (normal, Poisson, negagive binomial),</li>
    <li>user can select the underlying function for the sub-epidemic building block (generalized-logistic model, Richards model, and the generalized Richards model (GRM),</li>
    <li>user can select two different decline functions for sub-epidemic sizes: Exponential and Power-law.</li>
    
</ul>
    
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
