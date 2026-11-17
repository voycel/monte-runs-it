# 006: how many packets does it take to fill a sticker album?

The 2026 World Cup album has 980 stickers and a packet holds 7, so with perfect luck you'd need 140 packets. How many does it really take, where do they go, and what does swapping with a friend save?

`sim.py` runs every case 1,000,000 times from one seed, 20261117.

- **Album.** 980 different stickers, all equally likely. A packet holds 7 different stickers (no repeat inside a packet). Buy packets until every sticker is stuck in.
- **Warm-up.** A 10-sticker album bought one sticker at a time.
- **Two friends.** They pool every packet and swap perfectly, so they are done when every sticker has turned up at least twice.

Every simulated number is checked against its exact value from a Markov chain, and the run stops if one sits more than 4 standard errors away. The chain and the checks are in the docstring of `sim.py`.

Run it with `uv run -q sim.py`. Python 3.14.7 with numpy 2.5.5, 102.8 s on a Mac (`run.txt`).

## Results

| Question | Answer | File |
|---|---|---|
| Warm-up: stickers to fill a 10-sticker album | 29.27 (exact 10 × H₁₀ = 29.29) | `results.csv`, `warmup_stair.csv` |
| Packets to fill the 980-sticker album | 1,042.7 on average (exact 1,042.4), 7.4 × the perfect 140 | `results.csv`, `album_hist.csv` |
| Spread | median 1,013; luckiest 10% 846 or fewer; unluckiest 10% 1,275 or more; fewest 595, most 2,901 | `results.csv`, `album_hist.csv` |
| Where the packets go | first half 97; half to 90% 224; next 88 stickers 312; last 10 stickers 409 | `stages.csv` |
| The last 98 stickers (10%) | 721 packets, 69% of the total | `stages.csv` |
| The very last sticker | 140 packets on average (exact 140); 99.7% of albums spent at least one packet stuck at 979 | `results.csv` |
| Two friends who swap perfectly | 1,374 packets together, 687 each: 355 fewer each than alone | `results.csv`, `friends.csv` |

Before the first run I guessed about 900 packets on average, about 50 for the last sticker, a third of the packets on the last 10%, about 600 each for two friends, and 1,200 or more for the unluckiest 1 album in 10. Every guess was low: the last sticker alone takes as many packets as a perfect album.

At £1.25 a packet in the UK (about $2 in the US; Sports Illustrated, 13 May 2026), the average album costs about £1,303 or $2,085. No swapping group can average fewer than 140 packets per album, because every album still needs 980 stickers out of packets of 7.

No odds per sticker are published, so every sticker is assumed equally likely; whether a real packet can repeat a sticker isn't stated; Extra Stickers outside the 980 are ignored; and my friends swap perfectly.

## Files

`warmup.csv` and `warmup_stair.csv` (the 10-sticker warm-up), `albums.csv` and `album_hist.csv` (the first albums, and all 1,000,000 in 25-packet bins with the exact share), `album_conv.csv` (the running average), `shown_albums.csv` and `album1_packets.csv` (the albums opened sticker by sticker for the picture), `stages.csv`, `friends.csv`, `results.csv`, `run.txt`.
