"""Is the other queue really faster? 1,000,000 checkout trips and 1,000,000 Saturday customers.

Run 1 (race): k equal lines, 6 people each, you are last in line 0. Who gets through first? The exact
answers are symmetry: your line is first 1 time in k, and any one other line beats you 1 time in 2.

Run 2 (queues): the same customers and baskets, three tills, 90% busy, four queue rules. Three rules have
an exact average: random line (M/M/1), shared line (Erlang C, M/M/3) and shortest line with everyone
switching (M/M/3 too, because no till ever stands idle while someone waits). Shortest line without switching
has no closed form, so it is checked against the Markov chain solved numerically.

Run: python3 sim.py    Writes: results.csv, race.csv, convergence.csv, first-trips.csv, lines-window.csv, run.txt
"""
import csv
import heapq
import math
import platform
import random
import time
from collections import deque
from fractions import Fraction
from pathlib import Path

SEED = 20261030
HERE = Path(__file__).parent
TRIPS = 1_000_000
AHEAD = 6  # people in each line, you included
SERVICE = 2.0  # minutes per checkout on average
TILLS = 3
BUSY = Fraction(9, 10)
CUSTOMERS = 1_000_000
WARM_UP = 10_000  # the shop starts empty; these first customers are not counted
RATE = float(BUSY * TILLS / Fraction(SERVICE))  # customers a minute
WINDOW_START = 499_905  # about 70 minutes before a 14-person line, so the two hours drawn in the video run through it


def race(k, seed):
    rng = random.Random(seed)
    exp = rng.expovariate
    you_first = beaten_by_left = beaten_by_any = 0
    first = []
    conv = []
    marks = {round(10 ** (i / 4)) for i in range(25)}
    for trip in range(1, TRIPS + 1):
        times = [sum(exp(1 / SERVICE) for _ in range(AHEAD)) for _ in range(k)]
        best_other = min(times[1:])
        you_first += times[0] < best_other
        beaten_by_left += times[1] < times[0]
        beaten_by_any += best_other < times[0]
        if len(first) < 20:
            first.append(times)
        if trip in marks:
            conv.append((trip, you_first, beaten_by_any))
    return you_first, beaten_by_left, beaten_by_any, first, conv


def arrivals(seed):
    rng = random.Random(seed)
    t = 0.0
    out = []
    for _ in range(CUSTOMERS + WARM_UP):
        t += rng.expovariate(RATE)
        out.append((t, rng.expovariate(1 / SERVICE), rng.randrange(TILLS)))
    return out


def window(cs, starts, tills):
    """Line lengths once a simulated minute for two hours, starting at customer WINDOW_START."""
    t0 = cs[WINDOW_START][0]
    near = [(cs[i][0] - t0, starts[i] + cs[i][1] - t0, tills[i]) for i in range(WINDOW_START - 2000, WINDOW_START + 4000)]
    return [[sum(a <= m < b and till == i for a, b, till in near) for i in range(TILLS)] for m in range(121)]


def shared(cs):
    free = [0.0] * TILLS
    starts = []
    for t, s, _ in cs:
        at = heapq.heappop(free)
        start = max(t, at)
        heapq.heappush(free, start + s)
        starts.append(start)
    return starts, [0] * len(cs)  # one line, whichever till serves it


def random_line(cs):
    free = [0.0] * TILLS
    starts = []
    for t, s, pick in cs:
        start = max(t, free[pick])
        free[pick] = start + s
        starts.append(start)
    return starts, [c[2] for c in cs]


def shortest(cs, seed):
    rng = random.Random(seed)
    free = [0.0] * TILLS
    leaving = [deque() for _ in range(TILLS)]  # finish times of everyone in each line
    starts, tills = [], []
    for t, s, _ in cs:
        for q in leaving:
            while q and q[0] <= t:
                q.popleft()
        low = min(len(q) for q in leaving)
        pick = rng.choice([i for i, q in enumerate(leaving) if len(q) == low])
        start = max(t, free[pick])
        free[pick] = start + s
        leaving[pick].append(free[pick])
        starts.append(start)
        tills.append(pick)
    return starts, tills


def everyone_switches(cs, seed):
    """Shortest line, and whenever two lines differ by 2 the last person in the longer one moves."""
    rng = random.Random(seed)
    lines = [deque() for _ in range(TILLS)]  # customer indices; lines[i][0] is at the till
    done = []  # (finish time, till)
    starts = [0.0] * len(cs)

    def begin(i, now):
        c = lines[i][0]
        starts[c] = now
        heapq.heappush(done, (now + cs[c][1], i))

    def finish(now):
        _, i = heapq.heappop(done)
        lines[i].popleft()
        if lines[i]:
            begin(i, now)
        long, short = max(range(TILLS), key=lambda j: len(lines[j])), min(range(TILLS), key=lambda j: len(lines[j]))
        if len(lines[long]) - len(lines[short]) >= 2:
            lines[short].append(lines[long].pop())
            if len(lines[short]) == 1:
                begin(short, now)

    for c, (t, _, _) in enumerate(cs):
        while done and done[0][0] <= t:
            finish(done[0][0])
        low = min(len(q) for q in lines)
        pick = rng.choice([i for i, q in enumerate(lines) if len(q) == low])
        lines[pick].append(c)
        if len(lines[pick]) == 1:
            begin(pick, t)
    while done:
        finish(done[0][0])
    return starts, None


def erlang_c():
    a = BUSY * TILLS
    top = a ** TILLS / math.factorial(TILLS) / (1 - BUSY)
    return top / (sum(a ** n / math.factorial(n) for n in range(TILLS)) + top)


def exact_shared():
    c = erlang_c()
    rate = float((1 - BUSY) * TILLS) / SERVICE  # how fast the wait's tail falls away, per minute
    return float(c) / rate, math.log(float(c) / 0.01) / rate


def exact_random():
    rate = float(1 - BUSY) / SERVICE
    return float(BUSY) / rate, math.log(float(BUSY) / 0.01) / rate


def exact_shortest(limit=60):
    """Average wait for join-the-shortest-line: solve the chain on sorted line lengths (a <= b <= c)."""
    lam, mu = RATE, 1 / SERVICE
    states = [(a, b, c) for a in range(limit + 1) for b in range(a, limit + 1) for c in range(b, limit + 1)]
    index = {s: i for i, s in enumerate(states)}
    into = [[] for _ in states]  # (from state, rate)
    out = [0.0] * len(states)
    for s in states:
        a, b, c = s
        moves = [] if a == limit else [(tuple(sorted((a + 1, b, c))), lam)]
        for i, n in enumerate(s):
            if n:
                moves.append((tuple(sorted((*s[:i], n - 1, *s[i + 1:]))), mu))
        for to, r in moves:
            into[index[to]].append((index[s], r))
            out[index[s]] += r
    p = [1.0 / len(states)] * len(states)
    for sweep in range(20000):
        change = 0.0
        for i in range(len(states)):
            new = sum(p[j] * r for j, r in into[i]) / out[i]
            change = max(change, abs(new - p[i]))
            p[i] = new
        if change < 1e-13:
            break
    total = sum(p)
    queue = sum(x / total * sum(max(n - 1, 0) for n in s) for x, s in zip(p, states))
    assert sum(x for x, s in zip(p, states) if s[2] == limit) / total < 1e-6, "chain truncated too early"
    return queue / lam, sweep


def score(cs, starts):
    waits = [st - c[0] for c, st in zip(cs, starts)][WARM_UP:]
    latest = math.inf
    overtaken = 0
    for c_start in reversed(starts[WARM_UP:]):
        overtaken += latest < c_start
        latest = min(latest, c_start)
    order = sorted(waits)
    batches = [sum(waits[i:i + 10_000]) / 10_000 for i in range(0, CUSTOMERS, 10_000)]
    mean = sum(waits) / CUSTOMERS
    se = math.sqrt(sum((b - mean) ** 2 for b in batches) / (len(batches) * (len(batches) - 1)))
    return {"mean_wait_min": mean, "batch_se": se, "p99_wait_min": order[int(0.99 * CUSTOMERS)],
            "max_wait_min": order[-1], "no_wait": sum(w == 0 for w in waits), "overtaken": overtaken}


def main():
    t0 = time.perf_counter()
    rows = []
    for k in (2, 3, 5, 10):
        you, left, anyone, first, conv = race(k, SEED + k)
        for name, got, exact in (("you first", you, Fraction(1, k)), ("one named line beat you", left, Fraction(1, 2)),
                                 ("some other line beat you", anyone, Fraction(k - 1, k))):
            se = math.sqrt(float(exact * (1 - exact)) / TRIPS) * TRIPS
            z = (got - float(exact) * TRIPS) / se
            assert abs(z) < 4, (k, name, z)
            rows.append({"lines": k, "measure": name, "count": got, "trips": TRIPS, "exact": str(exact),
                         "exact_count": round(float(exact) * TRIPS), "z": round(z, 2)})
        if k == 3:
            with open(HERE / "first-trips.csv", "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(("trip", "your_line_min", "left_line_min", "right_line_min"))
                w.writerows((i, *(round(x, 2) for x in ts)) for i, ts in enumerate(first, 1))
            with open(HERE / "convergence.csv", "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(("trips", "you_first", "some_other_line_beat_you"))
                w.writerows(conv)
    with open(HERE / "race.csv", "w", newline="") as f:
        w = csv.DictWriter(f, rows[0].keys())
        w.writeheader()
        w.writerows(rows)

    cs = arrivals(SEED + 100)
    chain, sweeps = exact_shortest()
    rand, rand99 = exact_random()
    pool, pool99 = exact_shared()
    runs = (("random line", *random_line(cs), rand, rand99), ("shortest line", *shortest(cs, SEED + 101), chain, None),
            ("shared line", *shared(cs), pool, pool99), ("everyone switches", *everyone_switches(cs, SEED + 102), pool, None))
    out = []
    traces = []
    for name, starts, tills, exact, exact99 in runs:
        r = score(cs, starts)
        z = (r["mean_wait_min"] - exact) / r["batch_se"]
        assert abs(z) < 4, (name, z)
        out.append({"rule": name, **{k: round(v, 3) for k, v in r.items()}, "exact_mean_min": round(exact, 3),
                    "z_batches": round(z, 2), "exact_p99_min": round(exact99, 3) if exact99 else ""})
        if tills:
            traces += [(name, m, *row) for m, row in enumerate(window(cs, starts, tills))]
    with open(HERE / "lines-window.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(("rule", "minute", "line_1", "line_2", "line_3"))
        w.writerows(traces)
    with open(HERE / "results.csv", "w", newline="") as f:
        w = csv.DictWriter(f, out[0].keys())
        w.writeheader()
        w.writerows(out)
    (HERE / "run.txt").write_text(f"python {platform.python_version()} on {platform.platform()}\nseed {SEED}\n"
                                  f"markov chain sweeps {sweeps}\n")
    print(f"{time.perf_counter() - t0:.1f} s")


if __name__ == "__main__":
    main()
