"""Why do traffic jams appear out of nowhere? 1,000,000 loops of the Nagel-Schreckenberg model.

The model (Nagel and Schreckenberg, 1992): one lane cut into 7.5 m cells, one-second steps, top speed 5
cells a second (135 km/h). Every second, all drivers at once: speed up by one if below the top speed; slow
to the gap ahead; with chance P ease off by one for no reason; move.

Run 1 (loops): a 133-cell loop (997.5 m), cars evenly spaced and nobody braking. Does any car stop dead in
five minutes? 50,000 loops for each count of 10 to 29 cars: 1,000,000 loops.
Run 2 (flow): cars an hour past one point, for every count of cars on the loop, with and without easing off.
Run 3 (wave): a 7.5 km loop with heavy traffic; how far back does the stopped traffic sit 20 s later?
Run 4 (showcase): ten minutes of one 24-car loop, every car every second, with and without easing off, for the
picture. The first seed is used, not a picked one.

Checks (the run stops if one fails): with no easing off and top speed 5 the flow is exactly
min(5c, 1 - c) cars per cell per second (c = cars per cell), and with top speed 1 the flow is exactly
(1 - sqrt(1 - 4(1 - P)c(1 - c))) / 2 (Schadschneider and Schreckenberg, cond-mat/9511037, eq. 2).

Run: python3 sim.py    Writes: loops.csv, convergence.csv, flow.csv, check-vmax1.csv, wave.csv,
showcase.csv, results.csv, run.txt
"""
import csv
import math
import platform
import random
import time
from pathlib import Path

SEED = 20261103
HERE = Path(__file__).parent
CELL = 7.5  # metres
VMAX = 5  # cells a second
P = 0.2  # chance a driver eases off in any one second
LOOP = 133  # cells: 997.5 m
MINUTES = 5
LOOPS_EACH = 50_000
CARS = range(10, 30)  # 20 counts x 50,000 = 1,000,000 loops
KMH = CELL * 3.6  # one cell a second in km/h
CONV_CARS = 18  # the count whose share of stopped loops is drawn settling
SHOW_CARS = 24


def step(pos, vel, length, vmax, p, rand):
    """One second for every car at once. pos is in driving order, each car behind the next, on a loop."""
    n = len(pos)
    for i in range(n):
        v = vel[i] + 1
        if v > vmax:
            v = vmax
        g = (pos[(i + 1) % n] - pos[i] - 1) % length
        if v > g:
            v = g
        if p and v and rand() < p:
            v -= 1
        vel[i] = v
    for i in range(n):
        pos[i] = (pos[i] + vel[i]) % length


def even_start(n, length, vmax):
    pos = [i * length // n for i in range(n)]
    vel = [min(vmax, (pos[(i + 1) % n] - pos[i] - 1) % length) for i in range(n)]  # nobody braking
    return pos, vel


def first_stop(n, p, rand, seconds=MINUTES * 60):
    """Second of the first dead stop in a loop of n evenly spaced cars, or 0 if nobody stops."""
    pos, vel = even_start(n, LOOP, VMAX)
    for t in range(1, seconds + 1):
        step(pos, vel, LOOP, VMAX, p, rand)
        if 0 in vel:
            return t
    return 0


def loops():
    rows, conv = [], []
    marks = {round(10 ** (i / 4)) for i in range(19)} | {LOOPS_EACH}
    for n in CARS:
        rand = random.Random(SEED + n).random
        stops = []
        for k in range(1, LOOPS_EACH + 1):
            stops.append(first_stop(n, P, rand))
            if n == CONV_CARS and k in marks:
                conv.append((k, sum(1 for s in stops if s)))
        hit = sorted(s for s in stops if s)
        calm = first_stop(n, 0.0, None)  # no easing off: no randomness, one run says it all
        rows.append({
            "cars": n, "cars_per_km": round(n / (LOOP * CELL / 1000), 2), "loops": LOOPS_EACH,
            "stopped": len(hit), "share": round(len(hit) / LOOPS_EACH, 5),
            "median_first_stop_s": hit[len(hit) // 2] if hit else "",
            "stopped_without_easing_off": int(calm > 0),
        })
    return rows, conv


def flow(n, length, vmax, p, rand, warm, measure, start=None):
    """Mean cars per cell per second past a point, and per-batch values for a standard error."""
    pos, vel = start or (sorted(random.Random(n).sample(range(length), n)), [0] * n)
    for _ in range(warm):
        step(pos, vel, length, vmax, p, rand)
    batches, moved = [], 0
    for t in range(1, measure + 1):
        step(pos, vel, length, vmax, p, rand)
        moved += sum(vel)
        if t % (measure // 20) == 0:
            batches.append(moved / length / (measure // 20))
            moved = 0
    mean = sum(batches) / len(batches)
    se = (sum((b - mean) ** 2 for b in batches) / (len(batches) - 1) / len(batches)) ** 0.5
    return mean, se


def flows():
    rows = []
    for p in (0.0, P):
        rand = random.Random(SEED + 1000).random
        for n in range(1, LOOP + 1):
            f, se = flow(n, LOOP, VMAX, p, rand, 1000, 10_000)
            exact = min(VMAX * n / LOOP, 1 - n / LOOP) if p == 0 else None
            if exact is not None:
                assert abs(f - exact) < 1e-9, ("deterministic flow", n, f, exact)
            rows.append({
                "p": p, "cars": n, "cars_per_km": round(n / (LOOP * CELL / 1000), 2),
                "cars_per_hour": round(f * 3600), "se_cars_per_hour": round(se * 3600, 1),
                "mean_speed_kmh": round(f * LOOP / n * KMH, 1),
                "exact_cars_per_hour": round(exact * 3600) if exact is not None else "",
            })
    return rows


def check_vmax1():
    # Independent rings, not batches of one ring: the whole ring's flow drifts slowly, so batch errors
    # come out too small (a first try failed at z = -4.6 on one ring, z = -0.7 across 12 rings).
    rows, length, rings = [], 2000, 10
    rand = random.Random(SEED + 2000).random
    starts = random.Random(SEED + 2001)
    for c10 in range(1, 10):
        c = c10 / 10
        n = round(c * length)
        fs = [flow(n, length, 1, P, rand, 1000, 2000, start=(sorted(starts.sample(range(length), n)), [0] * n))[0]
              for _ in range(rings)]
        f = sum(fs) / rings
        se = (sum((x - f) ** 2 for x in fs) / (rings - 1) / rings) ** 0.5
        exact = (1 - math.sqrt(1 - 4 * (1 - P) * c * (1 - c))) / 2
        z = (f - exact) / se
        assert abs(z) < 4, ("vmax 1 flow", c, f, exact, z)
        rows.append({"density": c, "flow": round(f, 5), "exact": round(exact, 5), "se": round(se, 5), "z": round(z, 2)})
    return rows


def wave(lag=20, length=1000, n=300, seconds=5000):
    """Correlate where cars are stopped now with where they are stopped `lag` seconds later."""
    rand = random.Random(SEED + 3000).random
    pos, vel = sorted(random.Random(SEED + 3001).sample(range(length), n)), [0] * n
    for _ in range(2000):
        step(pos, vel, length, VMAX, P, rand)
    stopped = []
    for _ in range(seconds + lag):
        step(pos, vel, length, VMAX, P, rand)
        stopped.append({x for x, v in zip(pos, vel) if v == 0})
    score = {dx: 0 for dx in range(-60, 21)}
    for t in range(seconds):
        later = stopped[t + lag]
        for x in stopped[t]:
            for dx in score:
                score[dx] += (x + dx) % length in later
    best = max(score, key=score.get)
    return [{"lag_s": lag, "shift_cells": dx, "matches": m} for dx, m in score.items()], best


def showcase():
    """The same 24 evenly spaced cars twice: easing off with chance P (seed SEED + 4000), and never easing off."""
    seed = SEED + 4000
    rand = random.Random(seed).random
    pos, vel = even_start(SHOW_CARS, LOOP, VMAX)
    cpos, cvel = even_start(SHOW_CARS, LOOP, VMAX)
    rows = []
    for t in range(601):
        if t:
            step(pos, vel, LOOP, VMAX, P, rand)
            step(cpos, cvel, LOOP, VMAX, 0.0, None)
        rows += [{"t": t, "car": i, "cell": pos[i], "speed": vel[i], "calm_cell": cpos[i], "calm_speed": cvel[i]}
                 for i in range(SHOW_CARS)]
    assert min(r["calm_speed"] for r in rows) > 0  # without easing off nobody ever stops
    return rows, seed


def write(name, rows):
    with open(HERE / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main():
    t0 = time.time()
    checks = check_vmax1()
    flow_rows = flows()
    loop_rows, conv = loops()
    wave_rows, shift = wave()
    show_rows, show_seed = showcase()
    write("check-vmax1.csv", checks)
    write("flow.csv", flow_rows)
    write("loops.csv", loop_rows)
    write("convergence.csv", [{"loops": k, "stopped": s} for k, s in conv])
    write("wave.csv", wave_rows)
    write("showcase.csv", show_rows)
    best = max((r for r in flow_rows if r["p"] == P), key=lambda r: r["cars_per_hour"])
    write("results.csv", [
        {"measure": "loops run", "value": sum(r["loops"] for r in loop_rows)},
        {"measure": "jam wave speed km/h (backwards)", "value": round(-shift * KMH / 20, 1)},
        {"measure": "jam wave shift cells per 20 s", "value": shift},
        {"measure": "most cars an hour with easing off", "value": best["cars_per_hour"]},
        {"measure": "at cars on the loop", "value": best["cars"]},
        {"measure": "showcase seed", "value": show_seed},
        {"measure": "showcase first stop s", "value": min(r["t"] for r in show_rows if r["speed"] == 0)},
    ])
    (HERE / "run.txt").write_text(
        f"python {platform.python_version()} on {platform.platform()}\nseed {SEED}\n"
        f"seconds {time.time() - t0:.1f}\n")
    print("done", round(time.time() - t0, 1), "s")


if __name__ == "__main__":
    main()
