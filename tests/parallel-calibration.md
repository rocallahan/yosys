# Calibrating the number of threads in parallel algorithms

Every time we run a parallel algorithm on a design, we need to use some
heuristic to decide how many threads to use. If we use too many, then
hardware resource contention and memory overheads increase and can reduce performance.
If we use too few then we leave hardware resources unused that could have improved
performance. In general the optimal number of threads depends on the user's exact
hardware, design, and even what other work is running on the
computer, so we can't expect to get this perfectly right. However we still
need a methodology that gives generally good results without being
complicated.

## Methodology

At each site where we choose a number of threads, we count the number of "work units" that we're parallelizing over (e.g. cells in a module) and divide by some constant. The constant is
chosen by picking some (machine, design) pairs, running those combinations with various values of the constant for that site (keeping all other thread-count choices fixed), and picking the value of the
constant that seems close to optimal.
We expect this to be effective because decisions made at each site should be independent in their
performance effects. At any given point in time, the number of threads running should be governed
by exactly one site.

## Designs

Synthesis of the `fft64_width_64.il` benchmark is a good example: some passes deal with
113K cells, which is enough for parallel threads to be useful but not so large that experiments
are slow.

## Machines

For reproducibility and recalibration over time, we should use standard machines that
many people can access. We need machines that have access to many cores, but where we can
get stable performance results. Bare metal cloud machines are probably a good choice.

## Code
