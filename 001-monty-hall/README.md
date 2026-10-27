# 001: the Monty Hall game, 1,000,000 times

There are three doors, with a car behind one and goats behind the other two. You pick a door. The host knows where the car is, so he opens one of the other doors to show a goat. Then he offers you the last closed door. Should you stay or switch?

`sim.py` plays the game for real, a million times per variant. It also works out the exact answer separately, by counting every equally likely case with fractions, and checks each run against it. A run fails if the simulated share is more than 4 standard errors from the exact one.

## Results

| Variant | Seed | Stay wins | Switch wins | Car revealed | Switch share | Exact |
|---|---|---|---|---|---|---|
| Classic: 3 doors, host knows | 20261027 | 333,038 | 666,962 | 0 | 66.70% | 2/3 |
| Host forgets: opens a door at random | 20261028 | 332,606 | 333,662 | 333,732 | 50.08% of 666,268 | 1/2 |
| 100 doors: host opens 98 goat doors | 20261029 | 10,083 | 989,917 | 0 | 98.99% | 99/100 |

In the forgetful version, games where the host reveals the car aren't counted, because once you've seen the car there is nothing left to choose.

## Run it

```
python3 sim.py
```

It uses only the standard library. Running it rewrites the four data files below. On the machine that made the video (Python 3.15.1, macOS arm64) all three variants took about 12 s together.

The seed is 20261027 plus the variant's index, so a rerun on Python 3.15.1 plays the same games. Python only guarantees the same seeded `randrange` and `sample` within one version, so another version may draw different games. The shares should still pass the 4-standard-error check.

## Files

- `sim.py`: the simulation and the exact-answer check.
- `results.csv`: one row per variant, with the seed, counts, switch share, exact answer, z-score and whether the run passed the 4-standard-error check.
- `convergence.csv`: the switch share at log-spaced game counts from 1 to 1,000,000, used for the chart in the video.
- `first-games.csv`: the first 30 games of the classic run, door by door (doors are numbered 1 to 3).
- `run.txt`: the Python version and platform of the run shown in the video.

## Sources

- Marilyn vos Savant, Ask Marilyn, Parade, September 1990.
- John Tierney, "Behind Monty Hall's Doors: Puzzle, Debate and Answer?", The New York Times, 21 July 1991.
