Parallel Data Processing Example
================================

This is a very simple example that shows how one could use
Python's `multiprocessing` library to do parallel analyses of
data from single participants.

Please note that this is a silly example in which we merely
compute the average response time of each participant. Due to
the additional overhead of using parallel Processes, the example
likely runs slower than a serial example. However, if the
computations per participant become heavier, you should start to
see that using more parallel Processes should decrease the time
it takes for the analysis to complete. To illustrate this, I
added a bit in the `analyse_participants` function that computes
pi. This is completely useless, but does a nice job of adding
some computational load.

Feel free to use this example as a basis to do your own parallel
analyses with. I should warn you that, for simplicity's sake,
some parts of the script aren't the optimal way of doing things.

If you have any questions, feel free to email me at
edwin.dalmaijer@psy.ox.ac.uk, or tweet at me [@esdalmaijer](https://twitter.com/esdalmaijer)