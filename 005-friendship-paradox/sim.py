"""Why do your friends have more friends than you? 1,000,000 picks in three towns, plus buses and elevators.

Towns: 1,000 people. Each has a sociability score averaging 10; people i and j are friends with chance
min(1, w_i w_j / W), W the town's total score (Chung-Lu). Even town: every score is 10. Sociable town: log-normal
scores (sigma 1) rescaled to average 10. Clockwork town: everyone has exactly 10 friends (a ring where each person
knows the 5 nearest on each side, then 100,000 random swaps that keep everyone at 10).
One run: pick a random person with at least one friend, then one of their friends at random; compare friend counts.
1,000,000 runs per town.

Buses: one every 10 minutes on average; 1,000,000 riders arrive at random moments on a timeline of 1,000,000 buses.
Clockwork: every 10 minutes. Bunched: gaps alternate 2 and 18 minutes. Random: exponential gaps, mean 10.

Elevators: floors 1-10, I wait on floor 2. Each elevator rides 1 -> 10 -> 1 nonstop, one floor per time unit (an
18-unit loop), starting at a random point of the loop. 1,000,000 waits each with 1, 3 and 10 elevators.

Checks (the run stops if one fails; every simulated number must sit within 4 standard errors of its exact value):
  towns     exact over the drawn town: mean over people of (own count), of (1/k_i x sum of friends' counts),
            of (share of friends with more / fewer / equal friends)
  buses     mean wait = sum g^2 / (2 sum g) over the drawn timeline; mean gap you land in = sum g^2 / sum g;
            bunched: share landing in an 18-minute gap = 18/20
  elevators P(first elevator is going down) = 1/2 + 1/2 x (7/9)^k  (derivation in the README)

Run: python3 -I sim.py    Writes: towns.csv, edges.csv, picks.csv, degrees.csv, town_conv.csv, buses.csv,
riders.csv, elevators.csv, elevator_runs.csv, results.csv, run.txt
"""
import bisect
import csv
import math
import platform
import random
import time
from pathlib import Path

SEED = 20261106
HERE = Path(__file__).parent
RUNS = 1_000_000
PEOPLE = 1000
MEAN_SCORE = 10
SHOW = 20  # first runs kept one by one for the picture: the first seed, not picked ones
FLOORS, MY_FLOOR = 10, 2
LOOP = 2 * (FLOORS - 1)


def chung_lu(scores, rng):
    total = sum(scores)
    friends = [[] for _ in scores]
    for i in range(len(scores)):
        for j in range(i + 1, len(scores)):
            if rng.random() < min(1.0, scores[i] * scores[j] / total):
                friends[i].append(j)
                friends[j].append(i)
    return friends


def clockwork(rng, swaps=100_000):
    """Everyone exactly 10 friends: a ring of 5-each-side, then degree-keeping double-edge swaps."""
    edges = {(i, (i + d) % PEOPLE) for i in range(PEOPLE) for d in range(1, 6)}
    edges = {tuple(sorted(e)) for e in edges}
    lst = list(edges)
    done = 0
    while done < swaps:
        x, y = rng.randrange(len(lst)), rng.randrange(len(lst))
        (a, b), (c, d) = lst[x], lst[y]
        if rng.random() < 0.5:
            c, d = d, c
        e1, e2 = tuple(sorted((a, d))), tuple(sorted((c, b)))
        if len({a, b, c, d}) < 4 or e1 in edges or e2 in edges:
            continue
        edges -= {lst[x], lst[y]}
        edges |= {e1, e2}
        lst[x], lst[y] = e1, e2
        done += 1
    friends = [[] for _ in range(PEOPLE)]
    for a, b in edges:
        friends[a].append(b)
        friends[b].append(a)
    assert all(len(f) == 10 for f in friends)
    return friends


def gate(name, sim, exact, se):
    assert abs(sim - exact) <= 4 * se + 1e-12, (name, sim, exact, se)


def mean_sd(total, total_sq, n):
    m = total / n
    return m, math.sqrt(max(total_sq / n - m * m, 0) / n)


def write(name, rows):
    with open(HERE / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def town(name, friends, rng, res, picks, conv):
    k = [len(f) for f in friends]
    who = [i for i in range(PEOPLE) if k[i]]
    # exact, over the drawn town
    ex_own = sum(k[i] for i in who) / len(who)
    ex_friend = sum(sum(k[j] for j in friends[i]) / k[i] for i in who) / len(who)
    ex_more = sum(sum(k[j] > k[i] for j in friends[i]) / k[i] for i in who) / len(who)
    ex_fewer = sum(sum(k[j] < k[i] for j in friends[i]) / k[i] for i in who) / len(who)
    feld = sum(sum(k[j] for j in friends[i]) / k[i] > k[i] for i in who)  # people whose friends average more
    listed = sum(x * x for x in k) / sum(k)  # every friend of every person in one list
    s_own = s_own2 = s_f = s_f2 = more = fewer = 0
    for r in range(1, RUNS + 1):
        i = who[int(rng.random() * len(who))]
        j = friends[i][int(rng.random() * k[i])]
        s_own += k[i]; s_own2 += k[i] ** 2; s_f += k[j]; s_f2 += k[j] ** 2
        more += k[j] > k[i]; fewer += k[j] < k[i]
        if r <= SHOW:
            picks.append({"town": name, "run": r, "person": i, "person_friends": k[i], "friend": j,
                          "friend_friends": k[j]})
        if r < 1000 or r % 1000 == 0:
            conv.append({"town": name, "runs": r, "sum_person": s_own, "sum_friend": s_f, "friend_more": more})
    own, se_own = mean_sd(s_own, s_own2, RUNS)
    fr, se_fr = mean_sd(s_f, s_f2, RUNS)
    p_more, p_fewer = more / RUNS, fewer / RUNS
    gate(name + " own", own, ex_own, se_own)
    gate(name + " friend", fr, ex_friend, se_fr)
    gate(name + " more", p_more, ex_more, math.sqrt(ex_more * (1 - ex_more) / RUNS))
    gate(name + " fewer", p_fewer, ex_fewer, math.sqrt(ex_fewer * (1 - ex_fewer) / RUNS))
    r6 = lambda x: round(x, 6)
    res += [
        {"measure": f"{name} people", "value": PEOPLE},
        {"measure": f"{name} people with no friends", "value": PEOPLE - len(who)},
        {"measure": f"{name} friendships", "value": sum(k) // 2},
        {"measure": f"{name} average friends per person", "value": r6(sum(k) / PEOPLE)},
        {"measure": f"{name} most friends", "value": max(k)},
        {"measure": f"{name} fewest friends", "value": min(k)},
        {"measure": f"{name} median friends", "value": sorted(k)[PEOPLE // 2]},
        {"measure": f"{name} sim picker friends", "value": r6(own)},
        {"measure": f"{name} exact picker friends", "value": r6(ex_own)},
        {"measure": f"{name} sim picked friend friends", "value": r6(fr)},
        {"measure": f"{name} exact picked friend friends", "value": r6(ex_friend)},
        {"measure": f"{name} sim friend has more", "value": p_more},
        {"measure": f"{name} exact friend has more", "value": r6(ex_more)},
        {"measure": f"{name} sim friend has fewer", "value": p_fewer},
        {"measure": f"{name} exact friend has fewer", "value": r6(ex_fewer)},
        {"measure": f"{name} sim tie", "value": r6(1 - p_more - p_fewer)},
        {"measure": f"{name} people whose friends average more than them", "value": feld},
        {"measure": f"{name} friend list average", "value": r6(listed)},
    ]
    return k


def buses(name, gaps, rng, res, riders):
    times, t = [], 0.0
    for g in gaps:
        t += g
        times.append(t)
    s1, s2 = sum(gaps), sum(g * g for g in gaps)
    ex_wait, ex_gap = s2 / (2 * s1), s2 / s1
    w_sum = w_sq = g_sum = g_sq = long = 0
    for r in range(1, RUNS + 1):
        a = rng.random() * times[-1]
        b = bisect.bisect_right(times, a)
        wait, gap = times[b] - a, gaps[b]
        w_sum += wait; w_sq += wait * wait; g_sum += gap; g_sq += gap * gap; long += gap == 18
        if r <= SHOW:
            riders.append({"timetable": name, "rider": r, "arrive": round(a, 4), "wait": round(wait, 4),
                           "gap": round(gap, 4)})
    wait, se_w = mean_sd(w_sum, w_sq, RUNS)
    gap, se_g = mean_sd(g_sum, g_sq, RUNS)
    gate(name + " wait", wait, ex_wait, se_w)
    gate(name + " gap", gap, ex_gap, se_g)
    res += [{"measure": f"bus {name} mean gap", "value": round(s1 / len(gaps), 6)},
            {"measure": f"bus {name} sim wait", "value": round(wait, 6)},
            {"measure": f"bus {name} exact wait", "value": round(ex_wait, 6)},
            {"measure": f"bus {name} sim gap landed in", "value": round(gap, 6)},
            {"measure": f"bus {name} exact gap landed in", "value": round(ex_gap, 6)}]
    if name == "bunched":
        gate("bunched long", long / RUNS, 0.9, math.sqrt(0.09 / RUNS))
        res.append({"measure": "bus bunched sim share in 18-minute gap", "value": long / RUNS})
    return times


def next_pass(phase):
    """Time until this elevator next passes floor 2, and whether it is going down then."""
    up, down = MY_FLOOR - 1, LOOP - (MY_FLOOR - 1)  # up at loop time 1, down at 17
    if up < phase <= down:
        return down - phase, True
    return (up - phase) % LOOP, False


def elevators(n, rng, res, runs_out, conv):
    down = 0
    for r in range(1, RUNS + 1):
        t, d = min(next_pass(rng.random() * LOOP) for _ in range(n))
        down += d
        if r <= SHOW:
            runs_out.append({"elevators": n, "run": r, "wait": round(t, 4), "down": int(d)})
        if r < 1000 or r % 1000 == 0:
            conv.append({"elevators": n, "runs": r, "down": down})
    ex = 0.5 + 0.5 * (1 - 2 * 2 * (MY_FLOOR - 1) / LOOP) ** n  # short leg below me is 2(f-1) of the loop
    gate(f"elevators {n}", down / RUNS, ex, math.sqrt(ex * (1 - ex) / RUNS))
    res += [{"measure": f"elevators {n} sim first going down", "value": down / RUNS},
            {"measure": f"elevators {n} exact first going down", "value": round(ex, 6)}]


def layout(k, rng):
    """Hubs in the middle: radius falls with friend count, angle random. Positions in a unit square."""
    top = max(k)
    pos = []
    for x in k:
        rad = 0.06 + 0.88 * (1 - x / top) ** 1.5 * 0.5 + rng.random() * 0.04
        ang = rng.random() * 2 * math.pi
        pos.append((0.5 + rad * math.cos(ang), 0.5 + rad * math.sin(ang)))
    return pos


def main():
    t0 = time.time()
    rng = random.Random(SEED)
    res, picks, conv, nodes, edges, degrees = [], [], [], [], [], []
    raw = [rng.lognormvariate(0, 1) for _ in range(PEOPLE)]
    scale = MEAN_SCORE * PEOPLE / sum(raw)
    towns = {"even": chung_lu([MEAN_SCORE] * PEOPLE, rng),
             "sociable": chung_lu([x * scale for x in raw], rng),
             "clockwork": clockwork(rng)}
    for name, friends in towns.items():
        k = town(name, friends, random.Random(f"{SEED}-{name}"), res, picks, conv)
        for i, (x, y) in enumerate(layout(k, random.Random(f"{SEED}-layout-{name}"))):
            nodes.append({"town": name, "person": i, "friends": k[i], "x": round(x, 5), "y": round(y, 5)})
        edges += [{"town": name, "a": i, "b": j} for i, f in enumerate(friends) for j in f if i < j]
        hist = {}
        for x in k:
            hist[x] = hist.get(x, 0) + 1
        degrees += [{"town": name, "friends": x, "people": c, "times_listed_as_friend": x * c}
                    for x, c in sorted(hist.items())]

    riders = []
    brng = random.Random(f"{SEED}-buses")
    timetables = {"clockwork": [10.0] * RUNS,
                  "bunched": [2.0 if i % 2 == 0 else 18.0 for i in range(RUNS)],
                  "random": [brng.expovariate(1 / 10) for _ in range(RUNS)]}
    bus_rows = []
    for name, gaps in timetables.items():
        times = buses(name, gaps, random.Random(f"{SEED}-riders-{name}"), res, riders)
        bus_rows += [{"timetable": name, "bus": i + 1, "time": round(t, 4)} for i, t in enumerate(times[:40])]
    rnd = timetables["random"]
    assert abs(sum(rnd) / len(rnd) - 10) < 0.05  # the drawn timeline really averages 10 minutes

    ele_runs, ele_conv = [], []
    for n in (1, 3, 10):
        elevators(n, random.Random(f"{SEED}-elevators-{n}"), res, ele_runs, ele_conv)

    res.append({"measure": "runs per case", "value": RUNS})
    res.append({"measure": "seed", "value": SEED})
    write("towns.csv", nodes)
    write("edges.csv", edges)
    write("picks.csv", picks)
    write("degrees.csv", degrees)
    write("town_conv.csv", conv)
    write("buses.csv", bus_rows)
    write("riders.csv", riders)
    write("elevator_runs.csv", ele_runs)
    write("elevators.csv", ele_conv)
    write("results.csv", res)
    (HERE / "run.txt").write_text(
        f"python {platform.python_version()} on {platform.platform()}\nseed {SEED}\n"
        f"seconds {time.time() - t0:.1f}\n")
    print("done", round(time.time() - t0, 1), "s")


if __name__ == "__main__":
    main()
