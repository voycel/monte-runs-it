# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy==2.5.5"]
# ///
"""How many shuffles make a deck of cards random? 1,000,000 decks, riffle and overhand.

Deck: 52 cards, numbered 0 (top) to 51 (bottom), every deck starting in that order.

Riffle (the Gilbert-Shannon-Reeds model), done card by card: cut the deck after a binomial(52, 1/2) number of
cards, then drop the cards one at a time, taking the next card from the left packet with chance
(cards left in it) / (cards left in both). 1,000,000 decks, 12 riffles each.

Overhand, a simple model with one assumed number p: each of the 51 gaps between cards is a break point with chance
p, independently, and the packets are stacked in reverse order. Three hands: "neat" p = 0.15 (packets of about 6
cards), "lazy" p = 0.05 (about 15 cards) and "tiny" p = 0.5 (about 2 cards). 1,000,000 decks each. The p values
are assumptions.

Random baseline: 1,000,000 decks drawn uniformly at random, to show what each test reads on a truly random deck.

Three tests, each read as a total-variation distance from what a random deck gives (0 = same, 1 = nothing alike):
  top card     where the original top card ends up; random = each of the 52 places 1 time in 52
  rising runs  the number of rising sequences (runs of cards still in their starting order); random = Eulerian
               numbers A(52, r) / 52!. For riffles the chance of a deck depends only on this count (Bayer and
               Diaconis 1992), so this distance is the full distance between the shuffled deck and a random one.
               For the overhand it is one test, so it is only a floor.
  neighbours   how far apart two cards that started next to each other end up; random = 2(52 - d) / (52 x 51)

Checks (the run stops if one fails; a simulated number must sit within 4 standard errors of its exact value):
  riffle k = 1..12   share of decks with the top card still on top; mean rising sequences; mean number of
                     starting neighbours still together; and |simulated distance - exact distance| <= 0.005
  random decks       the same three means against 1/52, 26.5 and 51/52
  overhand pass 1    mean packets = 1 + 51 p; share of decks whose top card drops to the very bottom = p
There is no exact answer for the overhand after pass 1: those figures are simulated only.

Run: uv run -q sim.py    Writes: riffle.csv, riffle_top.csv, riffle_rs.csv, random.csv, overhand.csv,
overhand_top.csv, shown_riffle.csv, shown_overhand.csv, results.csv, run.txt
"""
import csv
import math
import platform
import time
from fractions import Fraction
from pathlib import Path

import numpy as np

SEED = 20261124
HERE = Path(__file__).parent
RUNS = 1_000_000
N = 52
RIFFLES = 12
HANDS = {"neat": (0.15, 80, 1), "lazy": (0.05, 600, 5), "tiny": (0.5, 1000, 10)}  # p, passes, test every this many passes (and passes 1-10)
TOP_SHOW = {"neat": (1, 7, 20, 50), "lazy": (7, 50, 150, 400), "tiny": (7, 50, 200, 600)}
IDX = np.arange(N, dtype=np.int8)
ROW = np.arange(RUNS)
FACT = math.factorial(N)


def gate(name, sim, exact, se):
    assert abs(sim - exact) <= 4 * se + 1e-12, (name, sim, exact, se)
    print(f"ok {name}: sim {sim:.6f} exact {exact:.6f} ({(sim - exact) / se if se else 0:+.2f} SE)")


def write(name, rows):
    with open(HERE / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


# ---- exact answers ----
def eulerian(n):
    """A[r-1] = number of orderings of n cards with exactly r rising sequences."""
    a = [1]
    for m in range(2, n + 1):
        a = [(k + 1) * (a[k] if k < len(a) else 0) + (m - k) * (a[k - 1] if k else 0) for k in range(m)]
    return a


EUL = eulerian(N)
assert sum(EUL) == FACT
U_RS = np.array([float(Fraction(a, FACT)) for a in EUL])  # random deck: chance of r = 1..52 rising sequences
U_NB = np.array([2 * (N - d) / (N * (N - 1)) for d in range(1, N)])  # random deck: neighbours end d apart


def rs_exact(k):
    """Chance of r = 1..52 rising sequences after k riffles (Bayer-Diaconis): A(n, r) C(2^k + n - r, n) / 2^(kn)."""
    a = 2 ** k
    return [Fraction(EUL[r - 1] * math.comb(a + N - r, N), a ** N) for r in range(1, N + 1)]


def tv_exact(k):
    return float(sum(abs(p - Fraction(e, FACT)) for p, e in zip(rs_exact(k), EUL)) / 2)


def top_exact(k):
    """Chance the top card ends in place j = 0..51 after k riffles: it is the first card to draw the lowest of
    52 random k-bit numbers in the inverse shuffle."""
    a = 2 ** k
    return np.array([float(sum(Fraction((a - d - 1) ** j * (a - d) ** (N - 1 - j), a ** N) for d in range(a)))
                     for j in range(N)])


def adj_exact(k):
    """Mean number of starting neighbours (card c directly above card c + 1) still together after k riffles."""
    a = 2 ** k
    return (N - 1) * float(Fraction(1, a) + sum(Fraction((a - g) ** (N - 1), a ** N) for g in range(1, a)))


# ---- the shuffles ----
def riffle(deck, rng, trace=None):
    cut = rng.binomial(N, 0.5, len(deck))
    left, right = np.zeros(len(deck), np.int64), cut.astype(np.int64)
    out = np.empty_like(deck)
    side = np.empty(N, np.int8)
    for t in range(N):
        n_left, n_right = cut - left, N - right
        take_left = rng.random(len(deck)) * (n_left + n_right) < n_left
        out[:, t] = deck[ROW, np.where(take_left, left, right)]
        left += take_left
        right += ~take_left
        side[t] = 0 if take_left[0] else 1
    if trace is not None:
        trace.append((int(cut[0]), side.copy()))
    return out


def overhand(deck, p, rng, trace=None):
    brk = rng.random((len(deck), N - 1), dtype=np.float32) < p  # a break between place i and place i + 1
    start = np.zeros(deck.shape, np.int8)
    start[:, 1:] = np.where(brk, IDX[1:], 0)
    start = np.maximum.accumulate(start, 1)  # first place of the packet each card sits in
    end = np.full(deck.shape, N - 1, np.int8)
    end[:, :-1] = np.where(brk, IDX[:-1], N - 1)
    end = np.minimum.accumulate(end[:, ::-1], 1)[:, ::-1]  # last place of that packet
    new = (N - 1 - end + IDX - start).astype(np.int8)  # packets restacked in reverse order
    out = np.empty_like(deck)
    np.put_along_axis(out, new, deck, 1)
    if trace is not None:
        trace.append(np.concatenate([[0], np.cumsum(brk[0])]).astype(int))
    return out, brk.sum(1) + 1


def tests(deck):
    pos = np.empty_like(deck)
    np.put_along_axis(pos, deck, np.broadcast_to(IDX, deck.shape), 1)  # pos[d, c] = place of card c
    top = np.bincount(pos[:, 0], minlength=N)
    rs = 1 + (pos[:, 1:] < pos[:, :-1]).sum(1)
    adj = (deck[:, 1:] == deck[:, :-1] + 1).sum(1)
    gap = np.abs(pos[:, 1:].astype(np.int16) - pos[:, :-1])
    nb = np.bincount(gap.ravel(), minlength=N)[1:]
    return {"top": top, "rs_hist": np.bincount(rs, minlength=N + 1)[1:], "rs": rs, "adj": adj, "nb": nb,
            "top_tv": 0.5 * np.abs(top / RUNS - 1 / N).sum(),
            "rs_tv": 0.5 * np.abs(np.bincount(rs, minlength=N + 1)[1:] / RUNS - U_RS).sum(),
            "nb_tv": 0.5 * np.abs(nb / nb.sum() - U_NB).sum(), "nb_mean": gap.mean()}


def fresh():
    return np.broadcast_to(IDX, (RUNS, N)).copy()


t0 = time.time()
rng = np.random.default_rng(SEED)
results = {"cards": N, "runs per case": RUNS, "seed": SEED, "riffles": RIFFLES,
           "random deck mean rising sequences exact": (N + 1) / 2,
           "random deck mean neighbours together exact": (N - 1) / N,
           "random deck mean neighbour distance exact": (N + 1) / 3,
           "random deck top card on top exact": 1 / N}

# ---- random baseline: what the tests read on truly random decks ----
t = tests(rng.permuted(fresh(), axis=1))
gate("random top on top", t["top"][0] / RUNS, 1 / N, math.sqrt((1 / N) * (1 - 1 / N) / RUNS))
gate("random mean rising", t["rs"].mean(), (N + 1) / 2, t["rs"].std() / math.sqrt(RUNS))
gate("random mean together", t["adj"].mean(), (N - 1) / N, t["adj"].std() / math.sqrt(RUNS))
FLOOR = {"top": t["top_tv"], "rs": t["rs_tv"], "nb": t["nb_tv"]}
results |= {"random decks top-card distance (noise floor)": t["top_tv"],
            "random decks rising distance (noise floor)": t["rs_tv"],
            "random decks neighbour distance (noise floor)": t["nb_tv"],
            "random decks fewest rising sequences": int(t["rs"].min()),
            "random decks most rising sequences": int(t["rs"].max())}
write("random.csv", [{"rising": r, "decks": int(t["rs_hist"][r - 1]), "exact_share": U_RS[r - 1]} for r in range(1, N + 1)])
for cap in (2, 4, 8, 16, 32):  # how often a random deck has this few rising sequences
    results[f"random deck chance of {cap} or fewer rising sequences exact"] = float(Fraction(sum(EUL[:cap]), FACT))

# ---- riffles ----
deck, trace, rows, top_rows, rs_rows, shown = fresh(), [], [], [], [], []
shown += [{"riffle": 0, "place": j + 1, "card": j + 1, "from": ""} for j in range(N)]
for k in range(1, RIFFLES + 1):
    deck = riffle(deck, rng, trace)
    t = tests(deck)
    ex_rs, ex_top = rs_exact(k), top_exact(k)
    mean_rs = float(sum(r * p for r, p in zip(range(1, N + 1), ex_rs)))
    tv, top_tv = tv_exact(k), 0.5 * np.abs(ex_top - 1 / N).sum()
    gate(f"riffle {k} top on top", t["top"][0] / RUNS, ex_top[0], math.sqrt(ex_top[0] * (1 - ex_top[0]) / RUNS))
    gate(f"riffle {k} mean rising", t["rs"].mean(), mean_rs, max(t["rs"].std(), 1e-9) / math.sqrt(RUNS))
    gate(f"riffle {k} mean together", t["adj"].mean(), adj_exact(k), t["adj"].std() / math.sqrt(RUNS))
    assert abs(t["rs_tv"] - tv) <= 0.005, (k, t["rs_tv"], tv)
    rows.append({"riffles": k, "distance_exact": tv, "distance_sim": t["rs_tv"], "top_distance_exact": top_tv,
                 "top_distance_sim": t["top_tv"], "top_on_top_exact": ex_top[0], "top_on_top_sim": t["top"][0] / RUNS,
                 "mean_rising_exact": mean_rs, "mean_rising_sim": t["rs"].mean(), "most_rising_possible": min(2 ** k, N),
                 "together_exact": adj_exact(k), "together_sim": t["adj"].mean(),
                 "neighbour_distance_sim": t["nb_tv"], "neighbour_gap_sim": t["nb_mean"]})
    top_rows += [{"riffles": k, "place": j + 1, "decks": int(t["top"][j]), "exact_share": ex_top[j]} for j in range(N)]
    rs_rows += [{"riffles": k, "rising": r, "decks": int(t["rs_hist"][r - 1]), "exact_share": float(ex_rs[r - 1])}
                for r in range(1, N + 1)]
    cut, side = trace[-1]
    shown += [{"riffle": k, "place": j + 1, "card": int(deck[0, j]) + 1, "from": "LR"[side[j]]} for j in range(N)]
    results[f"shown deck cut at riffle {k}"] = cut
    print(f"riffle {k}: distance sim {t['rs_tv']:.4f} exact {tv:.4f}; top card {t['top_tv']:.4f} / {top_tv:.4f}; "
          f"neighbours {t['nb_tv']:.4f}")
write("riffle.csv", rows)
write("riffle_top.csv", top_rows)
write("riffle_rs.csv", rs_rows)
write("shown_riffle.csv", shown)
R7 = rows[6]
results |= {"riffle 7 distance exact": R7["distance_exact"], "riffle 7 distance sim": R7["distance_sim"],
            "riffle 7 top-card distance sim": R7["top_distance_sim"],
            "riffle 7 neighbour distance sim": R7["neighbour_distance_sim"],
            "first riffle with top-card distance under 0.01 exact":
                next(r["riffles"] for r in rows if r["top_distance_exact"] < 0.01),
            "first riffle with distance under 0.5 exact": next(r["riffles"] for r in rows if r["distance_exact"] < 0.5),
            "first riffle with distance under 0.01 exact":
                next(k for k in range(1, 30) if tv_exact(k) < 0.01)}

# ---- overhand ----
orows, otop, oshown = [], [], []
for hand, (p, passes, every) in HANDS.items():
    deck, trace = fresh(), []
    crossed = {}
    for m in range(1, passes + 1):
        deck, packets = overhand(deck, p, rng, trace if m <= 8 else None)
        if m <= 8:
            oshown += [{"hand": hand, "pass": m, "place": j + 1, "card": int(deck[0, j]) + 1,
                        "packet_before": int(trace[-1][j]) + 1} for j in range(N)]
        if m == 1:
            gate(f"{hand} pass 1 mean packets", packets.mean(), 1 + (N - 1) * p, packets.std() / math.sqrt(RUNS))
            bottom = (deck[:, N - 1] == 0).mean()
            gate(f"{hand} pass 1 top card to the bottom", bottom, p, math.sqrt(p * (1 - p) / RUNS))
            results |= {f"{hand} break chance per gap (assumed)": p, f"{hand} mean packets per pass exact": 1 + (N - 1) * p,
                        f"{hand} mean packets per pass sim": packets.mean(),
                        f"{hand} mean cards per packet exact": N / (1 + (N - 1) * p)}
        if m <= 10 or m % every == 0:
            t = tests(deck)
            orows.append({"hand": hand, "passes": m, "top_distance_sim": t["top_tv"], "rising_distance_sim": t["rs_tv"],
                          "neighbour_distance_sim": t["nb_tv"], "mean_rising_sim": t["rs"].mean(),
                          "together_sim": t["adj"].mean(), "neighbour_gap_sim": t["nb_mean"]})
            for name, key, target in (("top-card", "top_tv", R7["top_distance_sim"]), ("rising", "rs_tv", R7["distance_sim"]),
                                      ("neighbour", "nb_tv", R7["neighbour_distance_sim"])):
                if name not in crossed and t[key] <= target:
                    crossed[name] = m
                    results[f"{hand} first tested pass matching 7 riffles on the {name} test"] = m
            for name, key in (("top-card", "top_tv"), ("rising", "rs_tv"), ("neighbour", "nb_tv")):
                if "flat " + name not in crossed and t[key] <= 2 * FLOOR[{"top_tv": "top", "rs_tv": "rs", "nb_tv": "nb"}[key]]:
                    crossed["flat " + name] = m
                    results[f"{hand} first tested pass within twice the noise floor on the {name} test"] = m
            if m in TOP_SHOW[hand]:
                otop += [{"hand": hand, "passes": m, "place": j + 1, "decks": int(t["top"][j])} for j in range(N)]
            if m in (7, 20, 50, 100, 150, 200, 400, 600, 1000):
                results |= {f"{hand} {m} passes top-card distance sim": t["top_tv"],
                            f"{hand} {m} passes rising distance sim": t["rs_tv"],
                            f"{hand} {m} passes neighbour distance sim": t["nb_tv"]}
            if len([k for k in crossed if not k.startswith("flat")]) == 3 and "all" not in crossed:
                crossed["all"] = m
                results[f"{hand} first tested pass matching 7 riffles on all three tests"] = m
            print(f"{hand} {m}: top {t['top_tv']:.4f} rising {t['rs_tv']:.4f} neighbours {t['nb_tv']:.4f}")
write("overhand.csv", orows)
write("overhand_top.csv", otop)
write("shown_overhand.csv", oshown)

write("results.csv", [{"measure": k, "value": round(v, 6) if isinstance(v, float) and v > 1e-4 else v}
                      for k, v in results.items()])
secs = time.time() - t0
(HERE / "run.txt").write_text(f"python {platform.python_version()} numpy {np.__version__} on {platform.platform()}\n"
                              f"seed {SEED}\nseconds {secs:.1f}\n")
print(f"done in {secs:.1f} s")
