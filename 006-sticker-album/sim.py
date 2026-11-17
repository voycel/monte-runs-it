# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy==2.5.5"]
# ///
"""How many packets does it take to fill a 980-sticker album? 1,000,000 albums, plus a warm-up and two friends.

Album: 980 different stickers, all equally likely. A packet holds 7 different stickers chosen at random (no repeat
inside a packet). Buy packets until every sticker is stuck in.

Warm-up: a 10-sticker album bought one sticker at a time. 1,000,000 albums, sticker by sticker.
One collector: 1,000,000 albums. Each packet's count of new stickers is drawn from the hypergeometric distribution
(7 drawn from "still missing" and "already have"), which is exactly what opening 7 different random stickers gives;
the album is tracked packet by packet. Recorded: packets to reach 490 (half), 882 (90%), 970 (all but 10) and 980,
and the packets spent waiting with exactly 979.
Two friends who swap perfectly: they pool every packet until both albums are full, i.e. every sticker has turned up
at least twice. State per packet: stickers seen 0 times, seen once. 1,000,000 pairs.
Shown albums: 20 more albums opened sticker by sticker (real sticker numbers) for the picture.

Checks (the run stops if one fails; every simulated number must sit within 4 standard errors of its exact value):
  warm-up      mean stickers = 10 x H_10; the wait for the k-th new sticker = 10 / (11 - k)
  one album    exact Markov chain on stickers collected: E[m] = (1 + sum_{h>=1} p_m(h) E[m+h]) / (1 - p_m(0)),
               p_m(h) = C(980-m, h) C(m, 7-h) / C(980, 7); first passage to >= g the same way; the whole
               distribution of packets by stepping the chain forward; the wait at 979 is geometric, mean 140
  two friends  exact Markov chain on (seen 0 times, seen once), solved the same way

Run: uv run -q sim.py    Writes: warmup.csv, warmup_stair.csv, albums.csv, album_hist.csv, album_conv.csv,
shown_albums.csv, album1_packets.csv, stages.csv, friends.csv, results.csv, run.txt
"""
import csv
import math
import platform
import random
import time
from pathlib import Path

import numpy as np

SEED = 20261117
HERE = Path(__file__).parent
RUNS = 1_000_000
N, K = 980, 7
SMALL = 10
SHOW = 20
GOALS = (490, 882, 970)  # half, 90%, all but the last 10
BIN = 25
COMB = [[float(math.comb(n, k)) for k in range(K + 1)] for n in range(N + 1)]
CNK = float(math.comb(N, K))


def gate(name, sim, exact, se):
    assert abs(sim - exact) <= 4 * se + 1e-12, (name, sim, exact, se)
    print(f"ok {name}: sim {sim:.6f} exact {exact:.6f} ({(sim - exact) / se if se else 0:+.2f} SE)")


def write(name, rows):
    with open(HERE / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def p_new(m):
    """Chance a packet brings h = 0..7 new stickers when m are already stuck in."""
    return [COMB[N - m][h] * COMB[m][K - h] / CNK for h in range(K + 1)]


P = [p_new(m) for m in range(N + 1)]


def first_passage(goal):
    """Exact expected packets from 0 stickers until at least `goal` are stuck in."""
    F = [0.0] * (N + K + 1)
    for m in range(goal - 1, -1, -1):
        p = P[m]
        F[m] = (1 + sum(p[h] * F[m + h] for h in range(1, K + 1) if m + h < goal)) / (1 - p[0])
    return F[0]


def exact_distribution():
    """P(album full after exactly t packets), stepping the chain forward until the tail is below 1e-13."""
    T = np.array(P)  # T[m, h]
    dist = np.zeros(N + 1)
    dist[0] = 1.0
    out, t = [0.0], 0
    while dist[:N].sum() > 1e-13:
        new = np.zeros(N + 1)
        for h in range(K + 1):
            new[h:] += (dist * T[:, h])[: N + 1 - h]
        t += 1
        out.append(new[N] - dist[N])
        dist = new
    return np.array(out)


def exact_friends():
    """Exact expected packets for two friends pooling: state (a seen 0 times, b seen once)."""
    E = [[0.0] * (N + 1 - a) for a in range(N + 1)]
    for a in range(N + 1):
        for b in range(N + 1 - a):
            if a == 0 and b == 0:
                continue
            r = N - a - b
            acc = 1.0
            self_p = COMB[r][K] / CNK if r >= K else 0.0
            for x in range(min(a, K) + 1):
                ca = COMB[a][x]
                for y in range(min(b, K - x) + 1):
                    if x == 0 and y == 0:
                        continue
                    z = K - x - y
                    if z > r:
                        continue
                    p = ca * COMB[b][y] * COMB[r][z] / CNK
                    acc += p * E[a - x][b - y + x]
            E[a][b] = acc / (1 - self_p)
    return E[N][0]


def warmup(res):
    rng = random.Random(f"{SEED}-warmup")
    tot = tot2 = 0
    step_sum = [0] * SMALL
    rows = []
    for r in range(1, RUNS + 1):
        have, n, last = set(), 0, 0
        while len(have) < SMALL:
            n += 1
            s = int(rng.random() * SMALL)
            if s not in have:
                have.add(s)
                step_sum[len(have) - 1] += n - last
                last = n
        tot += n
        tot2 += n * n
        if r <= SHOW:
            rows.append({"run": r, "stickers": n})
    mean = tot / RUNS
    se = math.sqrt((tot2 / RUNS - mean * mean) / RUNS)
    exact = SMALL * sum(1 / i for i in range(1, SMALL + 1))
    gate("warm-up mean", mean, exact, se)
    stair = []
    for k in range(1, SMALL + 1):
        p = (SMALL - k + 1) / SMALL
        ex = 1 / p
        sim = step_sum[k - 1] / RUNS
        gate(f"warm-up step {k}", sim, ex, math.sqrt((1 - p) / p / p / RUNS))
        stair.append({"new_sticker": k, "sim_wait": round(sim, 6), "exact_wait": round(ex, 6)})
    res += [{"measure": "warm-up album size", "value": SMALL},
            {"measure": "warm-up sim stickers", "value": round(mean, 6)},
            {"measure": "warm-up exact stickers", "value": round(exact, 6)}]
    write("warmup.csv", rows)
    write("warmup_stair.csv", stair)


def albums(res):
    rng = np.random.default_rng(SEED)
    m = np.zeros(RUNS, dtype=np.int64)
    t = np.zeros(RUNS, dtype=np.int64)
    reach = {g: np.zeros(RUNS, dtype=np.int64) for g in GOALS}
    at979 = np.zeros(RUNS, dtype=np.int64)  # packets bought while holding exactly 979
    seen979 = np.zeros(RUNS, dtype=bool)
    act = np.arange(RUNS)
    while act.size:
        mm = m[act]
        h = rng.hypergeometric(N - mm, mm, K)
        on979 = mm == N - 1
        at979[act[on979]] += 1
        seen979[act[on979]] = True
        m[act] = mm + h
        t[act] += 1
        for g in GOALS:
            hit = (mm < g) & (mm + h >= g)
            reach[g][act[hit]] = t[act[hit]]
        act = act[m[act] < N]
    T = t.astype(float)
    mean, se = T.mean(), T.std() / math.sqrt(RUNS)
    E0 = first_passage(N)
    gate("album mean packets", mean, E0, se)
    dist = exact_distribution()
    cdf = np.cumsum(dist)
    ex_mean_dist = float((np.arange(len(dist)) * dist).sum())
    assert abs(ex_mean_dist - E0) < 1e-6, (ex_mean_dist, E0)  # the two exact routes agree
    q = {p: int(np.searchsorted(cdf, p)) for p in (0.1, 0.5, 0.9)}
    for p, v in q.items():
        ex = cdf[v]
        gate(f"album share done by {v}", float((t <= v).mean()), ex, math.sqrt(ex * (1 - ex) / RUNS))
    ex1000 = cdf[1000]
    sim1000 = float((t <= 1000).mean())
    gate("album done by 1000", sim1000, ex1000, math.sqrt(ex1000 * (1 - ex1000) / RUNS))
    ex_pass = {g: first_passage(g) for g in GOALS}
    for g in GOALS:
        r = reach[g].astype(float)
        gate(f"album packets to {g}", r.mean(), ex_pass[g], r.std() / math.sqrt(RUNS))
    w = at979[seen979].astype(float)
    p = K / N
    gate("wait at 979", w.mean(), 1 / p, math.sqrt((1 - p) / p / p / w.size))
    # exact chance of ever holding exactly 979: forward first-visit over the chain
    visit = [0.0] * (N + 1)
    visit[0] = 1.0
    for s in range(N):
        pr = P[s]
        out = visit[s] / (1 - pr[0])  # leaves s once, to s+h with h >= 1
        for hh in range(1, K + 1):
            if s + hh <= N:
                visit[s + hh] += out * pr[hh]
    gate("albums that sat at 979", float(seen979.mean()), visit[N - 1],
         math.sqrt(visit[N - 1] * (1 - visit[N - 1]) / RUNS))

    hist, rows = np.bincount(t // BIN), []
    for i, c in enumerate(hist):
        lo = i * BIN
        rows.append({"from": lo, "to": lo + BIN - 1, "albums": int(c),
                     "exact_share": round(float(dist[lo:lo + BIN].sum()), 8)})
    write("album_hist.csv", rows)
    run_sum = np.cumsum(T)
    conv = [{"albums": r, "sum_packets": int(run_sum[r - 1])} for r in range(1, RUNS + 1)
            if r < 1000 or r % 1000 == 0]
    write("album_conv.csv", conv)
    write("albums.csv", [{"album": i + 1, "packets": int(t[i])} for i in range(SHOW)])
    stages = [{"stage": "first half (0 to 490)", "sim": reach[490].mean(), "exact": ex_pass[490]},
              {"stage": "half to 90% (490 to 882)", "sim": (reach[882] - reach[490]).mean(),
               "exact": ex_pass[882] - ex_pass[490]},
              {"stage": "next 88 (882 to 970)", "sim": (reach[970] - reach[882]).mean(),
               "exact": ex_pass[970] - ex_pass[882]},
              {"stage": "last 10 (970 to 980)", "sim": (t - reach[970]).mean(), "exact": E0 - ex_pass[970]}]
    for s in stages:
        s["sim"], s["exact"] = round(float(s["sim"]), 6), round(float(s["exact"]), 6)
    write("stages.csv", stages)
    r6 = lambda x: round(float(x), 6)
    res += [{"measure": "album stickers", "value": N},
            {"measure": "stickers per packet", "value": K},
            {"measure": "perfect-luck packets", "value": N // K},
            {"measure": "album sim mean packets", "value": r6(mean)},
            {"measure": "album exact mean packets", "value": r6(E0)},
            {"measure": "album exact mean stickers opened", "value": r6(K * E0)},
            {"measure": "album exact mean duplicates", "value": r6(K * E0 - N)},
            {"measure": "album exact median packets", "value": q[0.5]},
            {"measure": "album exact 10% luckiest packets", "value": q[0.1]},
            {"measure": "album exact 90% packets", "value": q[0.9]},
            {"measure": "album sim median packets", "value": int(np.median(t))},
            {"measure": "album sim fewest packets", "value": int(t.min())},
            {"measure": "album sim most packets", "value": int(t.max())},
            {"measure": "album sim done by 1000", "value": r6(sim1000)},
            {"measure": "album exact done by 1000", "value": r6(ex1000)},
            {"measure": "album exact packets to 490", "value": r6(ex_pass[490])},
            {"measure": "album exact packets to 882", "value": r6(ex_pass[882])},
            {"measure": "album exact packets to 970", "value": r6(ex_pass[970])},
            {"measure": "album sim packets to 490", "value": r6(reach[490].mean())},
            {"measure": "album sim packets to 882", "value": r6(reach[882].mean())},
            {"measure": "album sim packets to 970", "value": r6(reach[970].mean())},
            {"measure": "album sim wait at 979", "value": r6(w.mean())},
            {"measure": "album exact wait at 979", "value": r6(1 / p)},
            {"measure": "album sim share that sat at 979", "value": r6(seen979.mean())},
            {"measure": "album exact share that sat at 979", "value": r6(visit[N - 1])}]
    return E0


def friends(res, E0):
    rng = np.random.default_rng([SEED, 2])
    a = np.full(RUNS, N, dtype=np.int64)
    b = np.zeros(RUNS, dtype=np.int64)
    t = np.zeros(RUNS, dtype=np.int64)
    act = np.arange(RUNS)
    while act.size:
        aa, bb = a[act], b[act]
        x = rng.hypergeometric(aa, N - aa, K)
        rest = K - x
        y = np.zeros_like(x)
        need = rest > 0
        y[need] = rng.hypergeometric(bb[need], N - aa[need] - bb[need], rest[need])
        a[act] = aa - x
        b[act] = bb - y + x
        t[act] += 1
        act = act[(a[act] > 0) | (b[act] > 0)]
    T = t.astype(float)
    exact = exact_friends()
    gate("two friends packets", T.mean(), exact, T.std() / math.sqrt(RUNS))
    write("friends.csv", [{"pair": i + 1, "packets_together": int(t[i])} for i in range(SHOW)])
    res += [{"measure": "two friends sim packets together", "value": round(float(T.mean()), 6)},
            {"measure": "two friends exact packets together", "value": round(exact, 6)},
            {"measure": "two friends exact packets each", "value": round(exact / 2, 6)},
            {"measure": "two friends exact saving each", "value": round(E0 - exact / 2, 6)}]


def shown():
    """20 albums opened sticker by sticker, for the picture. Album 1 is kept packet by packet."""
    rng = random.Random(f"{SEED}-shown")
    rows, log = [], []
    for album in range(1, SHOW + 1):
        have, n = set(), 0
        while len(have) < N:
            n += 1
            for s in rng.sample(range(N), K):
                new = s not in have
                have.add(s)
                if album == 1:
                    log.append({"packet": n, "sticker": s + 1, "new": int(new), "have": len(have)})
        rows.append({"album": album, "packets": n})
    write("shown_albums.csv", rows)
    write("album1_packets.csv", log)


def main():
    t0 = time.time()
    res = []
    warmup(res)
    E0 = albums(res)
    friends(res, E0)
    shown()
    res.append({"measure": "runs per case", "value": RUNS})
    res.append({"measure": "seed", "value": SEED})
    write("results.csv", res)
    (HERE / "run.txt").write_text(
        f"python {platform.python_version()} numpy {np.__version__} on {platform.platform()}\nseed {SEED}\n"
        f"seconds {time.time() - t0:.1f}\n")
    print("done", round(time.time() - t0, 1), "s")


if __name__ == "__main__":
    main()
