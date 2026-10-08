# SpatialWavePredict Toolbox

**Fit epidemic curves, compare sub-epidemic models, and generate short-term ensemble forecasts in MATLAB.**

SpatialWavePredict represents an epidemic wave as the aggregate of overlapping sub-epidemics. This structure can describe trajectories that are difficult to capture with a single growth curve, including multiple peaks, prolonged plateaus, and damped oscillations. The toolbox combines parameter estimation, model ranking, parametric bootstrap uncertainty, forecast visualization, and performance evaluation.

[Quick start](#quick-start) · [Input data](#input-data) · [Configuration](#configuration) · [Outputs](#outputs) · [Citation](#citation)

> **Research software:** validate the numerical results and uncertainty estimates for your dataset and configuration. The [current limitations](#interpretation-and-current-limitations) apply to the implementation, not just to the interpretation of its plots.

## Overview

The workflow fits candidate models, selects and refits top candidates, ranks the accepted refits using the corrected Akaike information criterion (AICc), and generates individual-model and ensemble forecasts. Available components include generalized-logistic and Richards-type growth functions, exponential or power-law decline in successive sub-epidemic sizes, alternative estimation objectives and observation-error distributions, and bootstrap-based summaries.

**What “spatial wave” means:** sub-epidemics are latent components of an aggregate trajectory. They are not automatically identified geographic locations, transmission chains, or measured subpopulations. A data matrix may contain several regions, but the standard workflow fits one selected column at a time; it does not jointly estimate a spatial transmission network.

For the methodology and illustrated applications, see the [2024 tutorial][tutorial] and the [original framework paper][framework]. Published comparisons apply to the datasets, methods, and forecast horizons evaluated in those studies; they do not establish universal forecasting superiority.

### Default model and parameters

For the generalized-logistic building block (`flag1=1`) with exponential decline in successive sub-epidemic sizes (`typedecline2=[1]`), the cumulative size of an **active** sub-epidemic follows

```math
\frac{dC_j(t)}{dt} = rC_j(t)^p\left(1-\frac{C_j(t)}{K_j}\right), \qquad K_j = K_0 e^{-q(j-1)}.
```

| Symbol | Interpretation |
| --- | --- |
| $C_j(t)$ | Cumulative size of sub-epidemic $j$ at time $t$. |
| $r$ | Growth-rate coefficient shared across sub-epidemics. |
| $p$ | Early-growth scaling exponent; $p=1$ gives approximately exponential early growth when $C_j(t)\ll K_j$. |
| $K_0$ | Carrying capacity of the first sub-epidemic. |
| $q$ | Decline parameter across successive sub-epidemics: $q=0$ gives equal carrying capacities, while $q>0$ gives progressively smaller capacities. |

With threshold-triggered onset (`onset_fixed=0`), each subsequent component activates when the preceding component reaches the onset threshold, $C_{\mathrm{thr}}$ (`onset_thr` in the model function). Inactive components do not grow. See [the implemented model equations][ode] for the activation rule and alternative growth functions.

## Requirements and installation

The fitting and bootstrap workflow uses MATLAB and the following products:

| Product | Functions used by the toolbox |
| --- | --- |
| Optimization Toolbox | `fmincon`, `optimoptions` |
| Global Optimization Toolbox | `MultiStart`, `createOptimProblem`, `CustomStartPointSet` |
| Statistics and Machine Learning Toolbox | Random sampling functions such as `normrnd`, `poissrnd`, `nbinrnd`, and `datasample` |

See the MathWorks documentation for [constrained optimization][fmincon-doc], [MultiStart][multistart-doc], and [normal random sampling][normrnd-doc]. The current serial workflow does not require Parallel Computing Toolbox. A tested minimum MATLAB release and GNU Octave compatibility are not documented; record your environment using `ver`.

Download the repository using **Code → Download ZIP**, or clone it:

```bash
git clone https://github.com/gchowell/SpatialWavePredict-Toolbox.git
```

Open the repository root in MATLAB, then run:

```matlab
cd('spatialWave_subepidemicFramework code');
addpath(pwd);

if ~exist('input', 'dir'), mkdir('input'); end
if ~exist('output', 'dir'), mkdir('output'); end

ver
which options -all
which plotForecast_SW_subepidemicFramework -all
```

**Keep MATLAB's current folder set to `spatialWave_subepidemicFramework code`.** The functions use relative `./input/` and `./output/` paths. The repository also contains a root-level copy of `plotForecast_SW_subepidemicFramework.m`; use the implementation inside the code directory. Avoid adding multiple versions of this project recursively to the MATLAB path.

## Quick start

This historical example uses the bundled U.S. COVID-19 data and the configuration shown below. It is not a live data feed or a guarantee of exact reproduction of published numerical results.

### 1. Configure the example

Open the two configuration functions:

```matlab
edit options
edit options_forecast
```

In **`options.m`**, replace the corresponding assignments with these values. Retain the existing function declaration, global declarations, and other configuration logic; do not replace the whole file with this excerpt.

```matlab
% Data and dates
cumulative1 = 1;
outbreakx = 52;                 % U.S. aggregate in this example matrix
caddate1 = '05-11-2020';         % Calibration snapshot / forecast origin
cadregion = 'USA';
caddisease = 'coronavirus';
datatype = 'cases';
DT = 1;                        % Daily observations
datevecfirst1 = [2020, 2, 27];
datevecend1 = [2021, 5, 31];     % Later snapshot used for verification

% Fitting and uncertainty
smoothfactor1 = 7;
calibrationperiod1 = 90;
method1 = 0;                   % Least squares
dist1 = 0;                     % Normal observation-error draws
numstartpoints = 20;
B = 300;

% Candidate models
npatches_fixed = 3;
topmodelsx = 4;
flag1 = 1;                     % Generalized-logistic building block
onset_fixed = 0;
typedecline2 = [1];             % Exponential decline
```

In **`options_forecast.m`**, use:

```matlab
getperformance = 1;            % Evaluate against later observations
deletetempfiles = 0;
forecastingperiod = 30;         % 30 observation intervals: here, 30 days
weight_type1 = 1;              % Akaike relative-likelihood weights
```

These settings use two different files already included in `input/`:

```text
cumulative-daily-coronavirus-cases-USA-05-11-2020.txt
cumulative-daily-coronavirus-cases-USA-05-31-2021.txt
```

The first supplies calibration data; the second supplies later observations for forecast evaluation. Keep the fitting settings unchanged between fitting, plotting, and forecasting. Assigning a value only in MATLAB's Command Window is not a substitute for editing the configuration functions, which are called again by the workflow.

For an installation check, `B` can be temporarily reduced, for example to `20`. Such a run checks setup only: it is not sufficient evidence of stable tail quantiles, optimizer convergence, or interval coverage. Restore and assess the computational settings before substantive analysis. The smoothing and bootstrap caveat below also applies to this default-style example.

### 2. Fit, inspect, and forecast

Run these commands from the code directory:

```matlab
rng(12345, 'twister');
Run_SW_subepidemicFramework(52, '05-11-2020');

plotRankings_SW_subepidemicFramework(52, '05-11-2020');
plotFit_SW_subepidemicFramework(52, '05-11-2020');

rng(67890, 'twister');
plotForecast_SW_subepidemicFramework(52, '05-11-2020', 30, 1);
```

The forecast arguments are the selected data column, calibration date, horizon in observation intervals, and ensemble-weighting option. Other settings are read from the configuration functions. The four entry points can also be called without arguments to use their configured defaults.

**Run the fitting step with your current code before plotting.** Current plotting and forecasting functions require a matching `FinalFit-original-*.mat` manifest and its associated rank files. Older bundled results may not contain this information; an old `ABC-original-*.mat` search result is not a replacement. See the [ranking loader][ranking-loader]. Archive earlier output before rerunning an analysis, because files can be overwritten.

## Input data

### Matrix layout

Place a plain numeric, whitespace-delimited `.txt` file in `input/`. Rows represent successive observation intervals; columns represent distinct regions or groups. Do not include a header, dates, region names, or a separate time-index column. Select one column using the one-based integer `outbreakx`.

For example, this illustrates the layout of two incidence series, not a complete fitting example:

```text
12  4
18  7
27  9
35  14
43  20
```

### Filenames and dates

The [runner][runner] constructs the input filename from the options:

```text
Incidence:  <frequency>-<disease>-<datatype>-<region>-<mm-dd-yyyy>.txt
Cumulative: cumulative-<frequency>-<disease>-<datatype>-<region>-<mm-dd-yyyy>.txt
```

`<frequency>` is `daily` for `DT=1`, `weekly` for `DT=7`, or `yearly` for `DT=365`. The remaining fields correspond to `caddisease`, `datatype`, `cadregion`, and `caddate1`. Use two-digit months/days and a four-digit year consistently. Set `datevecfirst1` to the date represented by the first row.

**A filename date does not truncate the matrix.** The calibration file must actually end at the intended forecast origin. For retrospective evaluation, provide a separate later snapshot identified by `datevecend1`, with the same start date, column ordering, and observation frequency. That file must cover the full forecast horizon; it is read by [the verification-data loader][data-loader].

### Data preparation

Use regularly spaced, finite, nonnegative observations. Resolve missing intervals, reporting corrections, and changes in definitions before fitting, and document the handling rather than silently replacing missing values with zeros. Count-likelihood analyses require integer incidence observations.

Set `cumulative1=0` for incidence or `cumulative1=1` for cumulative counts. Cumulative input is converted as `[y(1); diff(y)]`, so the first cumulative value is treated as the first interval's incidence. Account for any pre-existing cumulative baseline when preparing a truncated series; cumulative decreases create negative incidence and require reconciliation.

The runner retains the most recent calibration window and removes leading zeros within that window. Check the resulting start date and usable sample size. All-zero inputs are not supported by the current preprocessing path. For AICc, the actual number of fitted observations must exceed the number of estimated parameters by more than one; eight observations alone do not guarantee an adequate dataset.

## Configuration

### Fitting and model selection: `options.m`

| Setting | Meaning |
| --- | --- |
| `calibrationperiod1` | Maximum number of recent observations used for fitting, before leading-zero removal. |
| `smoothfactor1` | Moving-mean span in observations; `1` disables smoothing. |
| `numstartpoints` | Candidate-search initialization setting; later refit/bootstrap stages use their own start-point procedures. |
| `B` | Requested number of parametric bootstrap realizations per selected model. |
| `npatches_fixed` | Maximum candidate sub-epidemic count during fitting, not necessarily the number simulated during forecasting. |
| `topmodelsx` | Requested number of selected models to refit; limited by available candidates. A one-component setting forces this to one. |
| `flag1` | Growth-function selector; `1` is generalized logistic, `2` is a generalized Richards-type variant, and `4` is Richards. |
| `onset_fixed` | `0`: threshold-triggered onset; `1`: simultaneous onset at time zero. See the implementation cautions. |
| `typedecline2` | `[1]`: exponential size decline; `[2]`: power-law decline; `[1 2]`: search both families. |

The implementation also includes generalized-growth (`flag1=0`) and classical-logistic (`flag1=3`) branches. The “Linear Model” comment for `flag1=3` in `options.m` does not match the implemented logistic equation. Consult [the model equations][ode] when using additional branches rather than relying on abbreviated labels.

### Estimation objective versus observation error

`method1` selects the fitting objective; `dist1` selects the observation-error distribution used in uncertainty calculations. They are distinct settings and must be consistent with the analysis.

| `method1` | Objective |
| --- | --- |
| `0` | Least squares |
| `1` | Poisson maximum likelihood |
| `2` | Pearson chi-squared objective |
| `3` | Negative-binomial likelihood: variance = mean + alpha × mean |
| `4` | Negative-binomial likelihood: variance = mean + alpha × mean² |
| `5` | Negative-binomial likelihood: variance = mean + alpha × mean^d |
| `6` | Sum of absolute deviations / Laplace branch; see limitations before use |

For `dist1`, `0` denotes normal errors, `1` Poisson errors, `2` a variance-to-mean negative-binomial formulation, `3`–`5` the corresponding negative-binomial formulations above, and `6` Laplace errors. Listing an option does not imply that every combination has been validated. Do not apply a count likelihood to moving-average values as though they were integer observations.

### Forecasting and ensembles: `options_forecast.m`

`forecastingperiod` counts **observations**, not always days: a value of `4` means four days for daily data or four weeks for weekly data. `getperformance=1` enables evaluation using the later snapshot. `deletetempfiles=0` retains intermediate forecast files.

| `weight_type1` | Ensemble weights |
| --- | --- |
| `-1` | Equal weights |
| `1` | Akaike relative-likelihood weights; the default |
| `0` | Legacy inverse-AICc weights; requires positive AICc and is not standard Akaike weighting |
| `2` | Inverse weighted interval score (WIS) weights calculated on calibration data; based on in-sample fit, not held-out forecasts |

The [AICc-weight helper][weights] uses the final selected-refit scores. With more than one selected model, the driver constructs ensembles of the top two, top three, and so on, up to the available requested count. A prior-forecast-WIS branch exists in the code, but the driver does not provide a complete rolling-origin weighting workflow; do not interpret that branch as validated historical-performance weighting.

## Outputs

Results are written to `output/`. Exact filenames encode many, but not all, analysis settings.

| File or pattern | Contents |
| --- | --- |
| `ABC-original-*.mat` | Candidate-search results; not the authoritative final refit ranking. |
| `FinalFit-original-*.mat` | Final ranking manifest for selected, accepted refits and their identities. |
| `modifiedLogisticPatch-original-*-rank-*.mat` | Rank-specific fitted parameters, trajectories, and bootstrap results. |
| `Forecast-*.mat` | Intermediate model trajectories and observation-scale forecast draws. |
| `AICc-topRanked-*.csv`, `parameters-topRanked-*.csv` | Ranking and parameter summaries. |
| `ranked(...)-*.csv`, `Ensemble(...)-*.csv` | Dated trajectories with observations, median, and lower/upper interval bounds. |
| `quantile-*.csv` | Predictive quantile tables. |
| `performance-*.csv` | Calibration or forecast performance summaries. |
| `doublingTimes-*.csv` | Derived doubling-time summaries. |

The dated trajectory tables use `year`, `month`, `day`, `data`, `median`, `LB`, and `UB`. They can contain both calibration and forecast rows; distinguish these using the forecast origin. Quantile tables must be aligned to the calibration length and target dates rather than assumed to contain forecast-only rows.

Performance routines calculate point-error and interval-based measures, including mean absolute error (MAE), mean squared error (MSE), prediction-interval coverage, and weighted interval score (WIS). Some horizon-indexed results summarize **the first h forecast observations together**, rather than the single observation at horizon h; check [the scoring routine][metrics] before combining results across origins. Assess interval coverage together with interval width and score, not coverage alone.

## Reproducibility and troubleshooting

Record the code revision, MATLAB/toolbox versions, both configuration files, raw data snapshots, preprocessing decisions, and random seeds. Keep complete output from each run together; do not mix rank files and manifests from different runs. A fixed seed controls random draws within a given workflow but does not guarantee identical results across releases or changed execution paths.

| Symptom | Check |
| --- | --- |
| Input file cannot be found | Current folder, filename fields, cumulative prefix, and four-digit year. |
| `MultiStart`, `fmincon`, or a sampling function is unavailable | Installed/licensed MATLAB products and path resolution. |
| Too many outputs from `options_forecast` | Use the forecast implementation in the code directory, not the duplicate at the repository root. |
| `SpatialWave:FinalFitResultsRequired` or fit-identity mismatch | Rerun fitting for the exact settings and regenerate all downstream results together. |
| Forecast horizon cannot be evaluated | Later snapshot coverage, `datevecfirst1`, `datevecend1`, selected column, and sampling frequency. |
| Invalid fits, nonfinite values, or too few observations | Missing/negative data, effective calibration length, model complexity, likelihood choice, and optimizer diagnostics. |

For bug reports, open a [GitHub issue][issues] with the code revision, `ver` output, relevant settings, exact command, and complete error message. Include a small shareable or synthetic dataset that reproduces the issue; do not upload confidential records.

## Tutorials and citation

The [tutorial paper][tutorial] provides the methodological background and illustrated examples. A [video walkthrough][video] is also available. Older materials may use different function names or output formats; follow the commands in this README for this repository layout.

### Citation

When using the toolbox, cite the tutorial and the methodological papers relevant to your analysis, and report the code revision used.

**Toolbox tutorial**  
Chowell G, Tariq A, Dahal S, Bleichrodt A, Luo R, Hyman JM. SpatialWavePredict: a tutorial-based primer and toolbox for forecasting growth trajectories using the ensemble spatial wave sub-epidemic modeling framework. *BMC Medical Research Methodology*. 2024;24:131. [doi:10.1186/s12874-024-02241-2][tutorial]

**Original framework**  
Chowell G, Tariq A, Hyman JM. A novel sub-epidemic modeling framework for short-term forecasting epidemic waves. *BMC Medicine*. 2019;17:164. [doi:10.1186/s12916-019-1406-6][framework]

**COVID-19 forecasting application**  
Chowell G, Rothenberg R, Roosa K, Tariq A, Hyman JM, Luo R. Sub-epidemic Model Forecasts During the First Wave of the COVID-19 Pandemic in the USA and European Hotspots. In: *Mathematics of Public Health*. 2022:85–137. [doi:10.1007/978-3-030-85053-1_5][application]

## License

Distributed under the **GNU General Public License, version 3**. See [LICENSE](LICENSE) for the terms. The software is provided without warranty.

[tutorial]: https://doi.org/10.1186/s12874-024-02241-2
[framework]: https://doi.org/10.1186/s12916-019-1406-6
[application]: https://doi.org/10.1007/978-3-030-85053-1_5
[video]: https://www.youtube.com/watch?v=qxuF_tTzcR8
[issues]: https://github.com/gchowell/SpatialWavePredict-Toolbox/issues
[fmincon-doc]: https://www.mathworks.com/help/optim/ug/fmincon.html
[multistart-doc]: https://www.mathworks.com/help/gads/multistart.html
[normrnd-doc]: https://www.mathworks.com/help/stats/normrnd.html
[runner]: spatialWave_subepidemicFramework%20code/Run_SW_subepidemicFramework.m
[data-loader]: spatialWave_subepidemicFramework%20code/getData.m
[ranking-loader]: spatialWave_subepidemicFramework%20code/loadSpatialWaveFinalRanking.m
[weights]: spatialWave_subepidemicFramework%20code/getSpatialWaveAICcWeights.m
[bootstrap]: spatialWave_subepidemicFramework%20code/fittingModifiedLogisticFunctionPatchMultiple.m
[ode]: spatialWave_subepidemicFramework%20code/modifiedLogisticGrowthPatch.m
[forecast]: spatialWave_subepidemicFramework%20code/plotForecast_SW_subepidemicFramework.m
[metrics]: spatialWave_subepidemicFramework%20code/computeforecastperformance.m

