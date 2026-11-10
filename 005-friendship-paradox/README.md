# 005: why do your friends have more friends than you?

Pick a random person, then one of their friends. Does the friend usually have more friends? And is it the same rule that makes you wait for buses and see elevators going the wrong way?

`sim.py` runs every case 1,000,000 times from one seed, 20261106.

- **Towns.** 1,000 people, each with a sociability score averaging 10. Two people are friends with chance min(1, w_i w_j / W), where W is the town's total score (a Chung–Lu random graph).
  - Even town: every score is 10.
  - Sociable town: log-normal scores (sigma 1), rescaled to average 10.
  - Clockwork town: everyone has exactly 10 friends.
  - One run picks a random person with at least one friend, then one of their friends at random, and compares friend counts.
- **Buses.** One every 10 minutes on average, and a rider turns up at a random moment. Clockwork: exactly every 10. Bunched: gaps of 2 and 18 minutes, alternating. Random: exponential gaps, mean 10.
- **Elevators.** Floors 1–10, and I wait on floor 2. Each elevator rides 1 → 10 → 1 nonstop, one floor per time unit (an 18-unit loop), starting at a random point of the loop. 1, 3 and 10 elevators.

Every simulated number is checked against its exact value, and the run stops if one sits more than 4 standard errors away. The checks are in the docstring of `sim.py`.

Run it with `python3 -I sim.py`. Python 3.15.1, standard library only, 4.1 s on a Mac (`run.txt`).

## Results

| Question | Answer | File |
|---|---|---|
| Even town: friends of the person picked, and of their friend | 10.15 and 11.09 | `results.csv` |
| Even town: the friend has more friends | 54.5% (exact 54.46%); fewer 37.3%, tie 8.2% | `results.csv` |
| Even town: friend counts | 2 to 22, median 10; 624 of 1,000 people have friends who average more than they do | `results.csv`, `degrees.csv` |
| Sociable town: friends of the person picked, and of their friend | 10.55 and 26.41 (exact 26.42) | `results.csv` |
| Sociable town: the friend has more friends | 74.5%; 869 of the 959 people with a friend have friends who average more | `results.csv` |
| Sociable town: friend counts | 0 to 131, median 6; 41 people have none | `results.csv`, `degrees.csv` |
| Clockwork town | 1,000,000 ties | `results.csv` |
| Clockwork buses: average wait | 5.00 minutes | `results.csv` |
| Bunched buses: average wait | 8.19 minutes (exact 8.2); 90.0% of riders land in an 18-minute gap | `results.csv` |
| Random buses: average wait, and the gap a rider lands in | 10.01 and 20.01 minutes, against an average gap of 10.00 | `results.csv` |
| First elevator at floor 2 going down: 1, 3, 10 elevators | 88.9%, 73.5%, 54.0% (exact 8/9, 73.53%, 54.05%) | `results.csv`, `elevator_runs.csv` |

Before the first run I guessed that in the even town the friend would have about 11 friends and more friends than the picker 45% of the time; in the sociable town, about 25 and 75%; and that the first of 3 elevators would be going down 70% of the time, and of 10, 55%. The friend counts and the elevators were close. The even-town 45% had the direction wrong: equal scores still give unequal friend counts, and chance alone is enough.

For comparison, Ugander, Karrer, Backstrom and Marlow found that 92.7% of about 721 million active Facebook users in May 2011 had fewer friends than their friends' average (https://arxiv.org/abs/1111.4503). In the sociable town it is 869 of 959, about 91%.

## Why the elevator formula is 1/2 + 1/2 (7/9)^k

An elevator passes floor 2 going down, then 2 time units later going up, then 16 units later going down again. So I arrive in the long stretch before a down-pass 16 times in 18, which gives 8/9 for one elevator.

With k elevators, the first to arrive is going up only if one elevator's up-pass comes at some time x in the first 2 units, and no other elevator passes before x. An elevator avoids both of its passes in [0, x) with probability (18 − 2x)/18. So

P(up) = k ∫₀² (1/18) ((18 − 2x)/18)^(k−1) dx = ½ (1 − (7/9)^k),

and P(down) = ½ + ½ (7/9)^k: 8/9 for one elevator, 73.53% for 3, 54.05% for 10, and ½ in the limit. Gamow and Stern's 1958 reasoning didn't hold for more than one elevator; Knuth worked it out in 1969.

These towns are made up, friendships are formed independently, and my elevators never stop to load anyone.

## Files

`towns.csv` (each person's friend count and place in the video's picture), `edges.csv` (every friendship), `degrees.csv` (people and friend-list spots at each friend count), `picks.csv` (the first picks in each town, in order), `town_conv.csv` and `elevators.csv` (running totals as the runs add up), `buses.csv` and `riders.csv` (the first buses and riders of each timetable), `elevator_runs.csv` (the first waits), `results.csv`, `run.txt`.
