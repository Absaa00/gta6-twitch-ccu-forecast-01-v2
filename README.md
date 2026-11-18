[![make.fun](https://make.fun/badge/gta6-twitch-ccu-forecast-01.svg)](https://make.fun/bounty/gta6-twitch-ccu-forecast-01)
[![rank](https://make.fun/badge/gta6-twitch-ccu-forecast-01/Absaa00/gta6-twitch-ccu-forecast-01-v2/rank.svg)](https://make.fun/bounty/gta6-twitch-ccu-forecast-01)
[![reward](https://make.fun/badge/gta6-twitch-ccu-forecast-01/Absaa00/gta6-twitch-ccu-forecast-01-v2/reward.svg)](https://make.fun/bounty/gta6-twitch-ccu-forecast-01)

# GTA VI launch-day Twitch CCU quantile forecaster (v2)

Entry for the [make.fun bounty](https://make.fun/bounty/gta6-twitch-ccu-forecast-01)
"Forecast GTA VI's launch-day Twitch viewers hour by hour" ($250, top 20).

`main.py` (Python 3.13, standard library only) reads every `*.in` case in
the input directory (schema `makefun/twitch-ccu-case/v1`) and writes one
`<id>.out` forecast per case: five quantiles (5/25/50/75/95%) of the Twitch
category's hourly viewer count for release-day hours 2–24.

## Method

Each past launch in the case's `comparables` contributes its release-day
hourly curve as log-ratios to its own hour-1 print. Per hour, the forecast
is the weighted empirical quantiles of those log-ratios, scaled back by the
case's hour-1 print. Weights combine:

- **scale** — a Gaussian kernel (bandwidth 1.0 in log space) on hour-1
  print similarity, so big launches shape the forecast of a big launch;
- **unlock pattern** — launches that unlocked at rolling midnight local
  time (like GTA VI) count double, because they ride the same timezone
  unlock wave (New Zealand → Europe → UK → US East → US West);
- **floor** — a small weight on every comparable keeps the pool from
  collapsing onto a couple of launches (0.05 normally, 0.01 once the
  big-launch gate below is active).

v2 improvements over v1 (validated by leave-one-out over the kit's 37
usable comparables, pinball loss exactly as scored in the rules):

- **big-launch gate** — when the hour-1 print is ≥ 50,000 viewers, only
  comparables whose own hour-1 print was ≥ 50,000 train the forecast;
- **interpolated weighted quantiles** — linear interpolation on Hazen
  plotting positions instead of stepping to the next order statistic;
- **normal blend** — quantiles blended 25% toward the normal quantiles
  implied by the same median and interquartile range;
- **hour smoothing** — the five log-quantile curves smoothed with a
  3-hour moving average.

Hours with too few comparables fall back to the house baseline shape
(flat median, spread widening over the day).

## Local validation

Leave-one-out over the kit's 37 usable comparables (same protocol as v1):

| | ALL | big launches (hour-1 ≥ 50k) |
|---|---|---|
| v2 | **0.446** | **0.201** |
| v1 | 0.466 | 0.202 |
| house baseline | 0.930 | 0.251 |

90% coverage 0.89 on big launches. Practice outputs are all valid and
byte-identical across repeated runs (~0.3 s for the four practice cases).

## Run

```sh
MAKEFUN_INPUT_DIR=input MAKEFUN_OUTPUT_DIR=out python3 main.py
```
