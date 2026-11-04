# 004: how many people until two share a birthday?

A room fills one person at a time. How many people does it take before two share a birthday, and does it change with real birthdays, which aren't spread evenly over the year?

`sim.py` runs four times, with 1,000,000 rooms each time. Each room stops at its first match, so one room answers every room size at once: a room that first matched at person 19 also has a match at 23, 40 and 80.

- Run 1, pairs: 365 equally likely birthdays; stop when two people share one.
- Run 2, my birthday: stop when someone has my birthday.
- Run 3, triples: stop when three people share one birthday.
- Run 4, real birthdays: each birthday is drawn with that calendar date's share of US births in 2000–2014 (366 dates, 29 February included); stop when two share.

Every run has an exact formula (they are in the docstring of `sim.py`), and the run stops if a simulated share lands more than 4 standard errors from it at any room size from 1 to 100 (1 to 2,000 for my birthday). The real-birthday formula is exact too: n! times the n-th elementary symmetric polynomial of the 366 date shares, worked out in fractions.

Run it with `python3 -I sim.py`. Python 3.15.1, standard library only, 20.4 s on a Mac (`run.txt`). Seed 20261104 for run 1, plus 1, 2 and 3 for runs 2–4.

## Results

| Question | Answer | File |
|---|---|---|
| Rooms of 23 with a shared birthday | 507,116 of 1,000,000 (50.7%; exact 50.73%) | `results.csv`, `curve.csv` |
| Rooms of 22 | exact 47.57% | `curve.csv` |
| People for 99% and 99.9% | 57 and 70 | `curve.csv` |
| Most common person to make the first match | exact: the 20th (32,320 rooms per million); the run put the 19th on top, but 19 to 21 sit inside the noise | `first_match.csv` |
| People until someone shares my birthday, for even odds | 253 (50.1% at 253; 6.1% at 23) | `mine.csv` |
| People until three share a birthday, for even odds | 88 (1.26% at 23) | `curve.csv` |
| Real US birthdays at 23 | run 50.75% vs 50.71% even (a gap smaller than the noise); exact 50.79% vs 50.73%, a gap of 0.06 of a point | `results.csv` |
| Still 23 for even odds with real birthdays? | yes | `results.csv` |

23 people make 253 pairs, and each pair is a chance to match. You are in only 22 of them, which is why someone sharing your own birthday takes 253 people.

The births data: 62,187,024 births. The busiest date is 12 September (191,024), the rarest without 29 February is 25 December (96,568), and 29 February has 41,869. 352 of the 365 dates sit within 10% of the average (`dates.csv`).

Before the first run I guessed that real birthdays would add about 1 point at 23, that the busiest date would be in mid-September and the rarest 25 December, and that the most common first match would be the 20th person. The dates and the 20th person were right. The 1-point guess was about 16 times too big.

These rooms assume independent birthdays and no twins.

## Data

`data/US_births_2000-2014_SSA.csv`: daily US births 2000–2014 from the Social Security Administration, compiled by FiveThirtyEight (https://github.com/fivethirtyeight/data/tree/master/births), CC BY 4.0. Their README is in `data/README-fivethirtyeight.md`. This file keeps its CC BY 4.0 licence; the MIT licence of this repository does not cover it.

## Files

`curve.csv` (share with a match at each room size, simulated and exact, all three questions), `mine.csv`, `first_match.csv`, `dates.csv`, `convergence.csv`, `rooms.csv` (the first 20 rooms of run 1, person by person: the first seed, not picked ones), `hist_steps.csv`, `flicker.csv`, `results.csv`, `run.txt`.
