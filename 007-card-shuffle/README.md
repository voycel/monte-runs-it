# 007: how many shuffles make a deck of cards random?

The famous answer is seven riffle shuffles (Bayer and Diaconis, 1992). What does seven really buy, how does a single card compare with the whole deck, and how many passes of the overhand shuffle most people use does it take to match?

`sim.py` shuffles 1,000,000 decks per case from one seed, 20261124.

- **Deck.** 52 cards, every deck starting in the same order.
- **Riffle.** The Gilbert–Shannon–Reeds model, done card by card: cut after a binomial(52, ½) number of cards, then drop cards one at a time, from a packet with chance (cards left in it) / (cards left in both). Twelve riffles.
- **Overhand.** Each of the 51 gaps between cards is a break with chance p, and the packets are restacked in reverse order. Three hands: p = 0.15 (packets of about 6 cards), 0.05 (about 15) and 0.5 (about 2). The three p values are assumptions.
- **Three tests,** each read as a total-variation distance from a truly random deck (0 = the same, 1 = nothing in common): where the starting top card ends up; how many rising runs the deck has; how far apart two starting neighbours end up.

For riffles, the chance of a deck depends only on its number of rising runs (Bayer and Diaconis), so the run test gives the distance for the whole deck, and it has an exact value. Every riffle figure with an exact answer is checked against it, and the run stops if one sits more than 4 standard errors away (45 checks; the formulas are in `sim.py`). The overhand has exact checks for the first pass only; after that its figures are simulated only.

Run it with `uv run -q sim.py`. Python 3.14.7 with numpy 2.5.5, 754.2 s on a Mac (`run.txt`).

## Results

| Question | Answer | File |
|---|---|---|
| Distance from random after 1–4 riffles | 1.00 every time: at most 2, 4, 8, 16 runs, and a random deck has 16 or fewer about once in 2 million | `riffle.csv`, `results.csv` |
| After 5, 6, 7 riffles | 0.92, 0.61, 0.33 (at seven: simulated 0.335, exact 0.334) | `riffle.csv`, `riffle_rs.csv` |
| After 8, 9, 10 riffles | 0.17, 0.09, 0.04; under 0.01 takes 13 | `riffle.csv`, `results.csv` |
| Starting top card still on top | 49.9% after one riffle (exact 50%), 6.5% after four, 2.3% after seven; fair is 1 in 52, 1.9% | `riffle.csv`, `riffle_top.csv` |
| Riffles until the top card's place is within 1% of flat | 10 | `results.csv` |
| Starting neighbours still together | 25.5 of 51 after one riffle, 1.19 after seven; a random deck has 0.98 | `riffle.csv` |
| Overhand, six-card packets: passes to match seven riffles | 15 on the run test, 20 on the top card, 39 on neighbours | `overhand.csv`, `results.csv` |
| Overhand, fifteen-card packets | about 115 passes (tested every 5) | `overhand.csv`, `results.csv` |
| Overhand, two-card packets | about 350 passes (tested every 10) | `overhand.csv`, `results.csv` |

Before the first run I guessed that a single card would be lost after six riffles (it takes ten) and that the overhand shuffle would need about 150 passes to match seven riffles (39 with six-card packets). The seven-riffle result itself I knew from the paper.

Three tests are a floor, not a proof of randomness: a deck has about 8 × 10^67 orders and a million decks cannot see most of them. Robin Pemantle (1989) proved that for this overhand model the passes needed grow with the square of the deck size, and his abstract says a 52-card deck needs more than 1,000 by the full measure. Real hands are not the model: a riffle that drops cards in clumps needs more than seven. Nothing here is advice for any card game.

## Files

`riffle.csv` (every figure per riffle, simulated and exact), `riffle_top.csv` (where the top card ended, 52 places per riffle), `riffle_rs.csv` (decks by number of rising runs per riffle, with the exact share), `random.csv` (the same for 1,000,000 truly random decks), `overhand.csv` (the three tests per pass and hand), `overhand_top.csv`, `shown_riffle.csv` and `shown_overhand.csv` (the decks shown in the picture), `results.csv`, `run.txt`.

## Sources

- Dave Bayer and Persi Diaconis, "Trailing the Dovetail Shuffle to its Lair", The Annals of Applied Probability 2(2), 294–313, 1992, doi 10.1214/aoap/1177005705.
- Robin Pemantle, "Randomization time for the overhand shuffle", Journal of Theoretical Probability 2, 37–49, 1989.
- Johan Jonasson, "The overhand shuffle mixes in Θ(n² log n) steps", The Annals of Applied Probability 16(1), 2006, arXiv:math/0501401.
