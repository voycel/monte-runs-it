# 002: is the other queue really faster?

There are three tills, each with a line of six people, and you're last in yours. Does your line finish first? And does the way a shop runs its queue change how long people wait?

`sim.py` answers both. The first run races equal lines 1,000,000 times. The second puts the same 1,000,000 customers, with the same baskets, through four queue rules at three tills. Each result is checked against an exact answer, and a run stops if it lands more than 4 standard errors away.

## Results

Run 1: three equal lines, six people each, checkouts random with a two-minute average. Both numbers are my assumption, not a measurement.

| Lines | Your line first | One named line beat you | Some other line beat you | Exact (some other line) |
|---|---|---|---|---|
| 2 | 499,263 | 500,737 | 500,737 | 1/2 |
| 3 | 332,766 | 500,830 | 667,234 | 2/3 |
| 5 | 199,983 | 500,459 | 800,017 | 4/5 |
| 10 | 99,668 | 500,213 | 900,332 | 9/10 |

Run 2: three tills, 90% busy, 1.35 customers a minute, 1,000,000 counted customers after a 10,000-customer warm-up.

| Rule | Mean wait (min) | Exact mean | 1-in-100 wait (min) | Overtaken |
|---|---|---|---|---|
| Random line, no looking | 17.677 | 18.0 | 86.72 | 629,611 |
| Shortest line | 6.026 | 5.965 | 33.884 | 402,161 |
| One shared line | 5.512 | 5.447 | 30.675 | 0 |
| Shortest line, everyone switches | 5.511 | 5.447 | 32.816 | 370,782 |

"Overtaken" counts customers who saw someone who arrived after them reach a till first. Shortest line has no closed form, so its exact mean comes from a Markov chain solved numerically.

Before the first run I guessed a shared line would beat the shortest line by 15% on average and 25% on the 1-in-100 wait. It came back about 9% on both.

These results assume random checkout times and that nobody gives up and leaves. The waits belong to a shop at 90% busy. Change those rules and the answers can change. The 2/3 only needs equal lines.

## Run it

```
python3 sim.py
```

It uses only the standard library. Running it rewrites the data files below. On the machine that made the video (Python 3.15.1, macOS arm64) it took about 80 s.

The seed is 20261030 plus a fixed offset per run, so a rerun on Python 3.15.1 draws the same numbers. Python only guarantees the same seeded `expovariate` within one version, so another version may draw different numbers. The checks should still pass.

## Files

- `sim.py`: the simulation and the exact-answer checks.
- `race.csv`: run 1, for 2, 3, 5 and 10 lines, with counts, exact answers and z-scores.
- `results.csv`: run 2, one row per rule.
- `convergence.csv`: the three-line share at log-spaced trip counts, used for the chart in the video.
- `first-trips.csv`: the first 20 trips with three lines, in minutes per line.
- `lines-window.csv`: two hours of line lengths, minute by minute, per rule. The window was chosen to run through a 14-person line, so it shows a busy stretch and not a typical one.
- `run.txt`: the Python version and platform of the run shown in the video.
- `answer.md`: the short answer in plain text.
