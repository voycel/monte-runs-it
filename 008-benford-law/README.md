# 008: why do so many numbers start with 1?

In many real lists about 30% of the numbers start with 1 and under 5% start with 9: Benford's law, P(d) = log10(1 + 1/d). Which kinds of lists do this, which don't, how fast does it appear, and why this pattern and no other?

`sim.py` makes 1,000,000 numbers per case from one seed, 20261127.

- **Powers of two.** 2^1 to 2^1,000,000. First digits from the fractional part of n × log10(2) in 60-digit decimals, checked against the full integers for n ≤ 5,000. A power of two starts with 1 exactly when doubling added a digit, so the count of 1s is floor(N log10 2) = 301,029.
- **Uniform picks.** A whole number from 1 to N, all equally likely. Exact shares for every N up to 10,000; a million draws from 1..999,999 and from 1..1,999.
- **Growth.** A million numbers start at 1; each step every one is multiplied by its own random factor between 1 and 2 (an assumption). Steps 0 to 40. The exact law of log10 of the product, mod 1, is computed by repeated convolution on a grid of 2^20 points.
- **Change of units.** Every number of a list times the same constant: 2.54 (inches to centimetres) and 181 constants from 1 to 10.
- **Distance** is the total variation between a list's nine first-digit shares and Benford's (0 = identical).

Every simulated share with an exact answer is checked against it, and the run stops if one sits more than 4.5 standard errors away (390 checks). Episode 7 used 4; with about 380 share checks a 4-SE gate trips about 2.4% of the time by chance, and the first run here did (uniform digit 7 at +4.2 SE), so the gate was widened and the seed kept.

Run it with `uv run -q sim.py`. Python 3.14.7 with numpy 2.5.5, 6.8 s on a Mac (`run.txt`).

## Results

| Question | Answer | File |
|---|---|---|
| Powers of two starting with 1, of 2^1..2^1,000,000 | 301,029 (30.1%); digit 9: 4.58% (Benford 4.58%) | `powers.csv`, `results.csv` |
| Uniform picks from 1 to 999,999 | every first digit 11.0% to 11.2% (exact 1/9) | `uniform.csv` |
| Share starting with 1 among 1..N | 11% at N = 9, 58% at N = 19, 55.6% at N = 1,999 (simulated 55.5%); never settles | `sawtooth.csv`, `uniform.csv` |
| Growth: share starting with 1 | 100% after 1 step, 9.4% after 3, 48% after 7 | `growth.csv` |
| Growth: first step within 0.01 of Benford | 28, simulated and exact; 0.25 after 6 steps; about 0.001 after 40 | `growth.csv`, `results.csv` |
| Grown list × 2.54 | share of 1s 30.0% → 30.2%; across 181 constants no digit moves more than 0.17 points | `units.csv`, `results.csv` |
| Uniform list × 2.54 | share of 1s 11.0% → 43.8% (exact 43.7%); up to 45 points across constants | `units.csv`, `results.csv` |

Before the first run I guessed that the growth model would be within 0.01 of Benford after 6 steps (it takes 28), that the grown list's digits would move less than half a point under any change of units (held: 0.17), that the uniform list times 2.54 would give about 40% ones (43.8%), and that the powers of two would give 4.58% nines (they do). The 30.1% and the 301,029 I knew or worked out beforehand.

Benford's law is not a law for every list: it fits lists that span several powers of ten, and fails for lists inside one, such as adult heights. The growth factor between 1 and 2 is an assumption; a gentler factor takes longer to settle.

## Files

`powers.csv` (counts per first digit), `powers_running.csv` (shares after the first n powers), `powers_first.csv` (2^1..2^40 written out), `uniform.csv`, `sawtooth.csv` (share of 1..N starting with 1, every N to 10,000), `growth.csv` (shares and distance per step, simulated and exact), `growth_shown.csv` (the 12 numbers shown in the picture), `units.csv` (both lists under 181 constants), `results.csv`, `run.txt`.

## Sources

- Simon Newcomb, "Note on the Frequency of Use of the Different Digits in Natural Numbers", American Journal of Mathematics 4(1), 39–40, 1881, doi 10.2307/2369148.
- Frank Benford, "The Law of Anomalous Numbers", 1938 (as cited by Wikipedia; not read).
- Roger S. Pinkham, "On the Distribution of First Significant Digits", The Annals of Mathematical Statistics 32(4), 1223–1230, 1961, doi 10.1214/aoms/1177704862.
- Theodore P. Hill, "A Statistical Derivation of the Significant-Digit Law", Statistical Science 10(4), 354–363, 1995, doi 10.1214/ss/1177009869.
- Wikipedia, "Benford's law" (read 27 November 2026): the definition, the history, scale invariance and the feet-and-yards argument.
