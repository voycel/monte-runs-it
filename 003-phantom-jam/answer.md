# Why do traffic jams appear out of nowhere?

Because a busy road turns a small wobble into a stop. On a simulated one-lane loop about a kilometre round, 24 evenly spaced cars with nothing in the way had a car stopped dead after 17 seconds. One driver eases off a little, the driver behind slows a little more, and each car further back loses a bit more until someone reaches zero. With perfect drivers who never ease off, nobody ever stopped.

How many cars it takes: across 1,000,000 simulated loops, a kilometre with 15 cars jammed about 1 time in 40 within five minutes, and with 19 cars 97 times in 100. Half the loops jammed somewhere between 17 and 18 cars.

The jam moves backwards, at about 19 km/h in this model, because cars leave its front while new ones join its back. Real drivers on a test track in Japan in 2008 produced the same backward wave at about 20 km/h. The wobble costs the road about a quarter of its capacity: the most cars an hour fell from 2,977 with perfect drivers to 2,188.

These are simulations of the Nagel–Schreckenberg model with one lane, identical drivers and an assumed 1-in-5 chance of easing off each second, not measurements of a real road. Code and data: this folder.
