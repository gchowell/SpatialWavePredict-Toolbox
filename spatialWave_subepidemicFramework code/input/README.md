# Input data

Each file is a header-free numeric matrix: one row per regularly spaced observation, one column per region or group. The drivers read files from this folder by name, so the name must match the settings in `options.m`.

## Naming convention

```
[cumulative-]<frequency>-<caddisease>-<datatype>-<cadregion>-<mm-dd-yyyy>.txt
```

- `cumulative-` is present when `cumulative1 = 1` (running totals) and absent for interval counts.
- `<frequency>` is `daily`, `weekly` or `yearly` (`DT = 1`, `7` or `365`).
- The date always has a four-digit year. The fitting file uses `caddate1`; the verification file uses `datevecend1`.
- Names are case-sensitive on Linux.

## Bundled datasets

| Files | Used by |
|---|---|
| `cumulative-daily-coronavirus-cases-USA-<date>.txt`: 43 snapshots, April 21 – July 6, 2020 | Quick start (`05-11-2020`), rolling-origin example (`04-27-2020`, `05-11-2020`, `05-25-2020`), presets `options_covid19.m`, `options_sw.m`, `options_tutorial.m` |
| `cumulative-daily-coronavirus-cases-USA-05-31-2021.txt` | Verification file for all U.S. COVID-19 examples |
| `cumulative-daily-coronavirus-cases-World13-05-31-2021.txt` | Additional example (13 world regions) |
| `cumulative-daily-FMD-cases-England-04-18-2001.txt`, `...-06-20-2001.txt` | Preset `options_fmd.m` (foot-and-mouth disease, England 2001) |
| `daily-SARS-cases-Singapore-02-25-2003.txt` | Preset `options_SARS.m` |
| `cumulative-weekly-monkeypox-cases-USA-<date>.txt`: five snapshots, July 2022 – July 2023 | Preset `options_mpx.m` (`11-10-2022`, verified against `01-26-2023`) |
| `cumulative-yearly-monkeypox-cases-USA-11-16-2049.txt`, `...-07-20-2079.txt` | Preset `options_yearly.m`. These are the weekly mpox files relabeled to demonstrate `DT = 365`; the dates are placeholders. |
| `yearly-HIV-cases-World_Regions-01-01-2019.txt` | Additional example |

## Data vintages

All 43 U.S. daily snapshots dated 2020 are the first rows of `...-05-31-2021.txt`, cut at each date. They contain counts as revised by May 31, 2021, not the counts that were available on each date. Forecasts scored with them, including the rolling-origin example, therefore use revised data. For a real-time evaluation, replace them with archived snapshots and rerun with `'refit', true`.
