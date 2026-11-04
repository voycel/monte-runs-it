"""How many people until two share a birthday? 1,000,000 rooms, four ways, plus real US birthdays.

Every room fills one person at a time and stops at the first match, so one room answers every room size at
once: a room that first matched at person 19 also has a match at 23, 40 and 80.

Run 1 (pairs): 365 equally likely birthdays; stop when two people share one. 1,000,000 rooms.
Run 2 (mine): stop when someone has MY birthday (day 0). 1,000,000 rooms.
Run 3 (triples): stop when three people share one birthday. 1,000,000 rooms.
Run 4 (real): birthdays drawn with the real share of US births on each calendar date, 2000-2014 (SSA via
FiveThirtyEight, 366 dates, 29 February included); stop when two share. 1,000,000 rooms.
The first 20 rooms of run 1, person by person, are kept for the picture (and the first 300 as day lists for the
running ring): the first seed, not picked ones.

Checks (the run stops if one fails; the simulated share must sit within 4 standard errors of the exact one
at every room size from 1 to 100, and for run 2 from 1 to 2,000):
  pairs     P(no match among n) = 365 x 364 x ... x (365 - n + 1) / 365^n
  mine      P(nobody has my birthday among n) = (364/365)^n
  triples   P(no day has three among n) = sum over k days with two people of
            365! / (k! (n - 2k)! (365 - n + k)!) x n! / 2^k / 365^n
  real      P(no match among n) = n! x e_n(p_1, ..., p_366), e_n the n-th elementary symmetric polynomial
            of the date shares

Run: python3 -I sim.py    Writes: curve.csv, mine.csv, first_match.csv, dates.csv, convergence.csv,
rooms.csv, hist_steps.csv, flicker.csv, results.csv, run.txt
"""
import bisect
import csv
import math
import platform
import random
import time
from fractions import Fraction
from itertools import accumulate
from pathlib import Path

SEED = 20261104
HERE = Path(__file__).parent
ROOMS = 1_000_000
DAYS = 365
UPTO = 100  # room sizes on the curve
MINE_UPTO = 2000
SHOW_ROOMS = 20
FLICKER_ROOMS = 300  # the year ring in M04 runs through these
MONTHS = "Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec".split()


def real_dates():
    """Total births on each calendar date, 2000-2014: [(month, day, births)], 366 dates."""
    total = {}
    with open(HERE / "data" / "US_births_2000-2014_SSA.csv") as f:
        for r in csv.DictReader(f):
            key = (int(r["month"]), int(r["date_of_month"]))
            total[key] = total.get(key, 0) + int(r["births"])
    assert len(total) == 366, len(total)
    return [(m, d, b) for (m, d), b in sorted(total.items())]


def first_match(draw, need=2):
    """People in the room when some birthday first reaches `need` people."""
    seen = {}
    n = 0
    while True:
        n += 1
        day = draw()
        c = seen.get(day, 0) + 1
        if c == need:
            return n
        seen[day] = c


def run(rooms, room_size):
    """Histogram of the room size at which each room stopped."""
    hist = {}
    for _ in range(rooms):
        n = room_size()
        hist[n] = hist.get(n, 0) + 1
    return hist


def cumulative(hist, upto):
    """Share of rooms that had stopped by each room size 1..upto."""
    out, c = [], 0
    for n in range(1, upto + 1):
        c += hist.get(n, 0)
        out.append(c / ROOMS)
    return out


def exact_pairs(upto):
    q, out = 1.0, []
    for n in range(1, upto + 1):
        q *= (DAYS - n + 1) / DAYS
        out.append(1 - q)
    return out


def exact_triples(upto):
    out = []
    for n in range(1, upto + 1):
        ways = sum(math.factorial(DAYS) // (math.factorial(k) * math.factorial(n - 2 * k) * math.factorial(DAYS - n + k))
                   * (math.factorial(n) // 2 ** k)
                   for k in range(0, n // 2 + 1) if n - k <= DAYS)
        out.append(float(1 - Fraction(ways, DAYS ** n)))
    return out


def exact_real(shares, upto):
    e = [1.0] + [0.0] * upto
    for p in shares:
        for k in range(upto, 0, -1):
            e[k] += e[k - 1] * p
    return [1 - math.factorial(n) * e[n] for n in range(1, upto + 1)]


def gate(name, sim, exact):
    for n, (s, x) in enumerate(zip(sim, exact), 1):
        se = math.sqrt(x * (1 - x) / ROOMS)
        assert abs(s - x) <= 4 * se + 1e-12, (name, n, s, x, se)


def first_at(curve, share):
    return next(n for n, p in enumerate(curve, 1) if p >= share)


def write(name, rows):
    with open(HERE / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main():
    t0 = time.time()
    dates = real_dates()
    births = sum(b for _, _, b in dates)
    shares = [b / births for _, _, b in dates]
    cum = list(accumulate(shares))

    # Run 1, keeping the first rooms person by person and the running share for the live counter.
    rand = random.Random(SEED).random
    draw = lambda: int(rand() * DAYS)
    show, conv, hist, by23, steps, flicker = [], [], {}, 0, [], []
    for room in range(1, ROOMS + 1):
        if room <= FLICKER_ROOMS:  # same draws as first_match, kept person by person
            seen, n, days = {}, 0, []
            while True:
                n += 1
                day = draw()
                match = day in seen
                days.append(day)
                if room <= SHOW_ROOMS:
                    show.append({"room": room, "person": n, "day": day, "matches_person": seen.get(day, "")})
                if match:
                    break
                seen[day] = n
            flicker.append({"room": room, "days": " ".join(map(str, days))})
        else:
            n = first_match(draw)
        hist[n] = hist.get(n, 0) + 1
        by23 += n <= 23
        if room < 1000 or room % 1000 == 0:
            conv.append({"rooms": room, "matched_by_23": by23})
        if room in (1000, 10_000, 100_000, ROOMS):  # the histogram as it stood, for the filling chart
            steps += [{"rooms": room, "person": k, "count": v} for k, v in sorted(hist.items())]
    pairs = cumulative(hist, UPTO)

    rand = random.Random(SEED + 1).random
    mine_hist = run(ROOMS, lambda: first_match_mine(rand))
    mine = cumulative(mine_hist, MINE_UPTO)

    rand = random.Random(SEED + 2).random
    triples = cumulative(run(ROOMS, lambda: first_match(lambda: int(rand() * DAYS), 3)), UPTO)

    rand = random.Random(SEED + 3).random
    real_hist = run(ROOMS, lambda: first_match(lambda: min(bisect.bisect(cum, rand()), 365)))
    real = cumulative(real_hist, UPTO)

    ex_pairs, ex_triples, ex_real = exact_pairs(UPTO), exact_triples(UPTO), exact_real(shares, UPTO)
    ex_mine = [1 - (364 / 365) ** n for n in range(1, MINE_UPTO + 1)]
    gate("pairs", pairs, ex_pairs)
    gate("mine", mine, ex_mine)
    gate("triples", triples, ex_triples)
    gate("real", real, ex_real)

    r6 = lambda x: round(x, 6)
    write("curve.csv", [{"people": n, "pairs": pairs[n - 1], "pairs_exact": r6(ex_pairs[n - 1]),
                         "triples": triples[n - 1], "triples_exact": r6(ex_triples[n - 1]),
                         "real": real[n - 1], "real_exact": r6(ex_real[n - 1])} for n in range(1, UPTO + 1)])
    write("mine.csv", [{"people": n, "mine": mine[n - 1], "mine_exact": r6(ex_mine[n - 1])}
                       for n in range(1, MINE_UPTO + 1)])
    write("first_match.csv", [{"person": n, "rooms": hist.get(n, 0), "real_rooms": real_hist.get(n, 0)}
                              for n in range(2, max(max(hist), max(real_hist)) + 1)])
    ranked = sorted(range(len(dates)), key=lambda i: -dates[i][2])
    rank = {i: r for r, i in enumerate(ranked, 1)}
    write("dates.csv", [{"month": m, "day": d, "births": b, "share": f"{s:.8f}", "rank": rank[i]}
                        for i, ((m, d, b), s) in enumerate(zip(dates, shares))])
    write("convergence.csv", conv)
    write("rooms.csv", show)
    write("hist_steps.csv", steps)
    write("flicker.csv", flicker)

    label = lambda i: f"{dates[i][1]} {MONTHS[dates[i][0] - 1]}"
    no_leap = [i for i in ranked if dates[i][:2] != (2, 29)]
    mode = max(hist, key=hist.get)
    # Exact chance the first match is person n: nobody matched among n - 1, then n hits one of their n - 1 days.
    first_exact = {n: (1 - ex_pairs[n - 2]) * (n - 1) / DAYS for n in range(2, UPTO + 1)}
    mode_exact = max(first_exact, key=first_exact.get)
    write("results.csv", [
        {"measure": "rooms per run", "value": ROOMS},
        {"measure": "pairs share at 23", "value": pairs[22]},
        {"measure": "pairs exact at 23", "value": r6(ex_pairs[22])},
        {"measure": "pairs exact at 22", "value": r6(ex_pairs[21])},
        {"measure": "pairs people for 50%", "value": first_at(ex_pairs, 0.5)},
        {"measure": "pairs people for 99%", "value": first_at(ex_pairs, 0.99)},
        {"measure": "pairs people for 99.9%", "value": first_at(ex_pairs, 0.999)},
        {"measure": "pairs most common first match person", "value": mode},
        {"measure": "pairs rooms first matching at mode", "value": hist[mode]},
        {"measure": "pairs exact most common first match person", "value": mode_exact},
        {"measure": "pairs exact rooms per million at exact mode", "value": round(first_exact[mode_exact] * ROOMS)},
        {"measure": "pairs rooms first matching at exact mode", "value": hist[mode_exact]},
        {"measure": "pairs largest room", "value": max(hist)},
        {"measure": "mine share at 23", "value": mine[22]},
        {"measure": "mine share at 253", "value": mine[252]},
        {"measure": "mine people for 50%", "value": first_at(ex_mine, 0.5)},
        {"measure": "mine largest room", "value": max(mine_hist)},
        {"measure": "triples share at 23", "value": triples[22]},
        {"measure": "triples people for 50%", "value": first_at(ex_triples, 0.5)},
        {"measure": "real share at 23", "value": real[22]},
        {"measure": "real exact at 23", "value": r6(ex_real[22])},
        {"measure": "real people for 50%", "value": first_at(ex_real, 0.5)},
        {"measure": "births counted", "value": births},
        {"measure": "most common date", "value": label(ranked[0])},
        {"measure": "most common date births", "value": dates[ranked[0]][2]},
        {"measure": "rarest date without 29 Feb", "value": label(no_leap[-1])},
        {"measure": "rarest date without 29 Feb births", "value": dates[no_leap[-1]][2]},
        {"measure": "29 Feb births", "value": dates[[i for i in ranked if dates[i][:2] == (2, 29)][0]][2]},
        {"measure": "seed", "value": SEED},
    ])
    (HERE / "run.txt").write_text(
        f"python {platform.python_version()} on {platform.platform()}\nseed {SEED}\n"
        f"seconds {time.time() - t0:.1f}\n")
    print("done", round(time.time() - t0, 1), "s")


def first_match_mine(rand):
    n = 1
    while int(rand() * DAYS) != 0:
        n += 1
    return n


if __name__ == "__main__":
    main()
