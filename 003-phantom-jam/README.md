# 003: why do traffic jams appear out of nowhere?

Put cars on an empty loop, evenly spaced, with nobody braking and nothing in the way. Does a jam still form? How many cars does a kilometre take before it does, and why does the jam move backwards?

`sim.py` answers with the Nagel–Schreckenberg traffic model (1992). One lane is cut into 7.5 m boxes, time moves in one-second steps, and the top speed is 5 boxes a second (135 km/h). Every second, all drivers at once speed up by one, slow to the gap ahead, ease off by one with chance 0.2 for no reason, and move. The 0.2 is my choice, not a measurement.

Two exact answers check the code, and the run stops if either fails:

- With no easing off, the flow must equal min(5c, 1 − c) cars per box per second. It matched within 1e-9 for every count from 1 to 133 cars.
- With top speed 1, the flow has an exact formula even with easing off (Schadschneider and Schreckenberg, arXiv cond-mat/9511037, eq. 2). All nine densities landed within 2 standard errors of it (`check-vmax1.csv`; the gate allows 4).

## Results

Run 1: a loop of 133 boxes (997.5 m). 50,000 loops for each count of 10 to 29 cars, 1,000,000 loops in all, five minutes each. Does any car stop dead?

| Cars on the loop | Loops with a stopped car | Share | Median first stop (s) |
|---|---|---|---|
| 10 | 1 | 1 in 50,000 | 203 |
| 15 | 1,291 | 2.6% | 193 |
| 16 | 4,774 | 9.5% | 175 |
| 17 | 15,214 | 30.4% | 160 |
| 18 | 35,091 | 70.2% | 129 |
| 19 | 48,671 | 97.3% | 77 |
| 20 | 49,993 | 99.99% | 44 |
| 24 | 50,000 | 100% | 14 |

Every count is in `loops.csv`. With no easing off, no car stopped at any count.

Run 2: cars an hour past one point on the loop (`flow.csv`).

| | Most cars an hour | At cars on the loop | With 30 cars |
|---|---|---|---|
| No easing off | 2,977 | 22 | 2,788 (92.7 km/h average) |
| Easing off 1 time in 5 | 2,188 | 18 | 1,859 (61.8 km/h average) |

Run 3: on a 7.5 km loop with 300 cars, stopped boxes line up best with stopped boxes 20 s later 14 boxes further back: the jam moves backwards at 18.9 km/h (the next-best shifts give ±1.4 km/h). A 2008 experiment with 22 real cars on a 230 m circuit measured about 20 km/h (Sugiyama et al., New Journal of Physics 10, 033001).

Run 4: ten minutes of one 24-car loop, every car every second, with and without easing off, for the picture in the video. It uses the first seed tried, not a picked one. Its first car stops dead at 17 s; across the 50,000 loops with 24 cars the median is 14 s.

Before the first run I guessed half of all loops would jam at about 18 cars, the jam would move backwards at about 10 km/h, and the most cars an hour would be about 2,000, at 20 cars. Half falls between 17 and 18 cars, the wave came back at 18.9 km/h, and the best flow was 2,188 at 18.

These results assume one lane, identical drivers and no lane changes. Real roads also jam at slip roads, behind lorries and at lane changes.

## Run it

```
python3 sim.py
```

It uses only the standard library. Running it rewrites the data files below. On the machine that made the video (Python 3.15.1, macOS arm64) it took about 300 s.

The seed is 20261103 plus a fixed offset per run, so a rerun on Python 3.15.1 draws the same numbers. Python keeps a seeded `random()` sequence the same across versions but not `sample()`, which places the cars for runs 2 and 3 and the top-speed-1 check, so another version may draw different numbers there. The checks should still pass.

## Files

- `sim.py`: the simulation and the exact-answer checks.
- `loops.csv`: run 1, one row per car count.
- `convergence.csv`: the 18-car share of stopped loops at log-spaced loop counts, used for the chart in the video.
- `flow.csv`: run 2, cars an hour and average speed for every car count, with and without easing off, and the exact value where one exists.
- `check-vmax1.csv`: the top-speed-1 check, nine densities against the exact formula.
- `wave.csv`: run 3, how many stopped boxes match 20 s later at each backward or forward shift.
- `showcase.csv`: run 4, every car's box and speed each second, with (`cell`, `speed`) and without (`calm_cell`, `calm_speed`) easing off.
- `results.csv`: the headline numbers.
- `run.txt`: the Python version and platform of the run shown in the video.
- `answer.md`: the short answer in plain text.
