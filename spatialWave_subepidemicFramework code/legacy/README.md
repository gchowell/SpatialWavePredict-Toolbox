# Legacy code

These scripts come from earlier versions of the toolbox. None of the current drivers or tests calls them, and they are not on the MATLAB path when you run the toolbox from the code folder. They read output filenames and option layouts that the current version no longer uses, so they are kept for reference only.

| File | What it was |
|---|---|
| `bagging_forecasts.m` | Earlier routine that combined ranked forecasts by resampling their trajectories and scored the result. The current ensembles use weighted quantiles instead. |
| `combineforecasts.m` | Script that reads the active options and calls `bagging_forecasts.m`. It clears the workspace and closes all figures. |

The earlier candidate-search function `fittingModifiedLogisticFunctionPatchABC_old_09042025.m` was removed; it remains in the repository's Git history.
