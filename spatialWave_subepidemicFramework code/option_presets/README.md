# Option presets

Earlier configuration files, kept for reference. They are **not on the MATLAB
path** and the toolbox does not read them: it reads `options.m` and
`options_forecast.m` in the code directory.

To use one, copy its assignments into `options.m` (or copy the file over
`options.m`). Each of these files defines a function called `options` (or
`options_forecast`), whatever the file is called, which is why they must not
sit next to the real `options.m`.

| File | Note |
| --- | --- |
| `options_SARS.m`, `options_SARS2003Singapore.m`, `options_covid19.m`, `options_fmd.m`, `options_mpx.m`, `options_sw.m`, `options_tutorial.m`, `options_yearly.m` | Current layout of `options.m` (first output `cumulative1`). Check the settings against the current `options.m` before use; settings added later may be missing. |
| `options_mpox.m`, `options_mpx2.m`, `optionsx.m` | **Older layout** (no `cumulative1` output). They do not work with the current code as they are. |
| `options_forecast_sw.m` | A copy of `options_forecast.m`. |
