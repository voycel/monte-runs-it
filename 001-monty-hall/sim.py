"""Monte Runs It, episode 001: the Monty Hall game, played 1,000,000 times per variant.

Every game is played literally: hide the car, pick a door, Monty opens doors, then compare staying
with switching to the one other closed door. The exact answer is computed separately by enumerating
every equally likely case with fractions, and the run fails if the two disagree.

Run: python3 sim.py   (standard library only; the same seed gives the same results)
"""
import csv
import math
import platform
import random
import time
from fractions import Fraction
from pathlib import Path

SEED = 20261027
GAMES = 1_000_000
VARIANTS = [  # name, doors, does Monty know where the car is?
    ("classic", 3, True),
    ("monty-forgets", 3, False),
    ("100-doors", 100, True),
]
OUT = Path(__file__).parent


def play(doors, monty_knows, rng):
    """One game: (car, pick, opened doors, the other closed door), or None if Monty revealed the car."""
    car = rng.randrange(doors)
    pick = rng.randrange(doors)
    others = [d for d in range(doors) if d != pick]
    # Monty opens all but one of the other doors; a knowing Monty only opens goat doors.
    pool = [d for d in others if d != car] if monty_knows else others
    opened = set(rng.sample(pool, doors - 2))
    if car in opened:
        return None
    (last,) = [d for d in others if d not in opened]
    assert (pick == car) != (last == car), "exactly one strategy wins each game"
    return car, pick, opened, last


def exact(doors, monty_knows):
    """Exact (P(switch wins | car not revealed), P(car not revealed)) by full enumeration."""
    switch = valid = Fraction(0)
    for car in range(doors):
        for pick in range(doors):
            others = [d for d in range(doors) if d != pick]
            # Monty leaves exactly one other door closed; a knowing Monty never opens the car.
            keep = [k for k in others if car in (k, pick)] if monty_knows else others
            for k in keep:
                if car not in (k, pick):
                    continue  # a forgetful Monty opened the car
                w = Fraction(1, doors * doors * len(keep))
                valid += w
                switch += w * (k == car)
    return switch / valid, valid


def main():
    results, convergence, first_games = [], [], []
    checkpoints = {round(10 ** (e / 4)) for e in range(25)}  # 1 … 1,000,000, log-spaced
    for i, (name, doors, knows) in enumerate(VARIANTS):
        rng = random.Random(SEED + i)
        stay = switch = spoiled = 0
        t0 = time.perf_counter()
        for g in range(1, GAMES + 1):
            game = play(doors, knows, rng)
            if game is None:
                spoiled += 1
            elif game[1] == game[0]:
                stay += 1
            else:
                switch += 1
            if name == "classic" and g <= 30 and game:
                car, pick, opened, last = game
                first_games.append({"game": g, "car": car + 1, "pick": pick + 1,
                                    "monty_opens": "+".join(str(d + 1) for d in sorted(opened)),
                                    "other_closed": last + 1,
                                    "winner": "stay" if pick == car else "switch"})
            if g in checkpoints and stay + switch:
                convergence.append({"variant": name, "games": g, "valid": stay + switch,
                                    "switch_share": f"{switch / (stay + switch):.6f}"})
        seconds = time.perf_counter() - t0

        p_exact, valid_exact = exact(doors, knows)
        n = stay + switch
        p_sim = switch / n
        se = math.sqrt(p_exact * (1 - p_exact) / n)
        z = (p_sim - p_exact) / se
        # Correctness gate: the simulated share must sit within 4 standard errors of the exact answer.
        ok = abs(z) <= 4
        if not knows:
            se_v = math.sqrt(valid_exact * (1 - valid_exact) / GAMES)
            ok = ok and abs(n / GAMES - valid_exact) / se_v <= 4
        results.append({"variant": name, "doors": doors, "monty_knows": knows, "seed": SEED + i,
                        "games": GAMES, "car_revealed": spoiled, "counted": n,
                        "stay_wins": stay, "switch_wins": switch,
                        "switch_share": f"{p_sim:.6f}", "exact_switch": str(p_exact),
                        "exact_decimal": f"{float(p_exact):.6f}", "z": f"{z:+.2f}",
                        "gate": "PASS" if ok else "FAIL", "seconds": f"{seconds:.1f}"})
        print(f"{name:14} stay {stay:>7,}  switch {switch:>7,}  car shown {spoiled:>7,}  "
              f"switch share {p_sim:.4f}  exact {p_exact} ({float(p_exact):.4f})  z {z:+.2f}  "
              f"{'PASS' if ok else 'FAIL'}  {seconds:.1f}s")

    for path, rows in [("results.csv", results), ("convergence.csv", convergence),
                       ("first-games.csv", first_games)]:
        with open(OUT / path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=rows[0].keys())
            w.writeheader()
            w.writerows(rows)
    (OUT / "run.txt").write_text(f"python {platform.python_version()} on {platform.platform()}\n"
                                 f"seed {SEED} (+ variant index), {GAMES:,} games per variant\n")
    if any(r["gate"] == "FAIL" for r in results):
        raise SystemExit("correctness gate FAILED")


if __name__ == "__main__":
    main()
