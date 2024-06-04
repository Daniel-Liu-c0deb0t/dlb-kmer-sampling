# Optimal distance lower bound k-mer sampling
Guaranteed larger distance between consecutive $k$-mers for the same density.

## Motivation
Given a set of $k$-mers (substrings of length $k$) derived from a long (genomic) sequence,
it is often desirable to subsample this $k$-mer set and perform computation (eg., comparing genomes, aligning reads, etc.) on these
smaller sets for reduced memory usage and runtime.

With the reasonable assumption that all $k$-mers in the sequence is unique,
the *density* of a sampling scheme is the number of sampled $k$-mers divided by the total number of $k$-mers.
In practice, density is usually fixed to be a specific based on the desired amount of subsampling.
One simple way to subsample $k$-mers is to hash each $k$-mer and pick those with hashes
(for simplicity, assume hashes are uniform between 0 and 1) less than or equal to the density.
This is commonly called the FracMinHash scheme.

An issue with the FracMinHash is that it does not enforce an upper bound or lower bound on the
distances between $k$-mers that are sampled consecutively in a long sequence.
Windowed schemes like minimizers provide a distance *upper bound* by sampling at least one $k$-mer per sliding window
in the sequence.
In this work, we examine the problem of sampling $k$-mers with a distance *lower bound*.

Intuitively, with a higher lower bound on the distances between $k$-mers while maintaining the same density,
the $k$-mers are more evenly spaced along a sequence, resulting in less clumping of $k$-mers and hopefully
retaining more information about the sequence within the same budget of the number of $k$-mers selected.

## Results
First, we show the main results on a random sequence over the alphabet {A, C, G, T} of length $10^6$
and a fixed density of $d = \frac{1}{11}$. We plot the (relative) frequency for each possible distance
between two consecutively sampled $k$-mers.

![Plot of distances for each algorithm](distances.png)

With our algorithm, which we simply call DLB (for Distance Lower Bound), the distance lower bound increases as $k$ increases.
This is not true for the other open syncmer schemes.
Note that the upper bound on the distance lower bound is $\frac{1}{d}$ as $k \to \infty$.
In this case, DLB with $k = 51$ and $k = 61$ both achieve a distance lower bound of 9, which is closer to $\frac{1}{d} = 11$.

It is interesting to note that increasing the distance lower bound suppresses distances
greater than the lower bound (see the DLB distances in the plot). This makes sense, because
the same density has to be maintained.

## Algorithm
### Previous work
Our DLB algorithm is based on *open syncmers*:
* [**Edgar (2021)**](https://doi.org/10.7717/peerj.10805) introduced the open syncmers scheme. A $k$-mer is sampled iff
the minimum hashed $s$-mer ($s \leq k$) within the $k$-mer is at the start (position 0) of the $k$-mer.
The density of this scheme is $\frac{1}{k - s + 1}$, since there are $k - s + 1$
positions for the $s$-mer. The distance lower bound for this scheme is 1, since nothing prevents
two neighboring $k$-mers in the sequence from both being sampled.

* [**Shaw & Yu (2022)**](https://doi.org/10.1093/bioinformatics/btab790) showed that the optimal position that maximizes conservation
is at the middle position: $\left\lfloor \frac{k - s}{2} \right\rfloor$. This also gives a distance lower bound.
Assume that the current $k$-mer is sampled. Then the minimum $s$-mer is at the middle position.
Two events can happen that samples another $k$-mer, assuming $s$-mers are distinct:
either the minimum $s$-mer is shifted out the left side of the $k$-mer,
or a smaller minimum $s$-mer is shifted in on the right side of the $k$-mer.
Therefore, the distance lower bound is the minimum of the distance between the middle and the left and right sides
of the $k$-mer: $\min(\left\lfloor \frac{k - s}{2} \right\rfloor, \left\lceil \frac{k - s}{2} \right\rceil) + 1 = \left\lfloor \frac{k - s}{2} \right\rfloor + 1$.

* [**Dutta et al. (2022)**](https://doi.org/10.1371/journal.pcbi.1010638) proposed the concept of parameterized syncmer schemes, which generalizes
open syncmers and allows specifying a set of positions where the $k$-mer will be sampled if the
minimum $s$-mer occurs in any of those positions. Although our novel scheme is a parameterized syncmer scheme,
we optimize for a different metric and we provide a sequence agnostic construction that can effectively use more
than 2 or 3 positions in practice.

DLB can be seen as a generalization of Shaw & Yu (2022)'s ideas to parameterized syncmer schemes.

### DLB
We are given $k$ and density $d$. Define a parameter $t$ which we will show how to set later.
Then, compute
```math
s = k - t \cdot \frac{1}{d} + 1.
```
Notice that as $t$ increases, $s$ decreases, causing density to decrease by $t$ times for the original
open syncmer scheme.
To compensate for this and maintain the density $d$, we chose $t$ positions such that if the
minimum $s$-mer (break ties by choosing the rightmost one) in the current $k$-mer occurs at any of the $t$ positions,
then we sample the $k$-mer.
These positions are maximally spaced apart in the $k - s + 1$ $s$-mers.
First, compute
```math
\begin{align}
stride &= \left\lfloor \frac{k - s + 1 - t}{t + 1} \right\rfloor\\
&= \left\lfloor \frac{t}{t + 1} \left( \frac{1}{d} - 1 \right) \right\rfloor.
\end{align}
```
Then, the positions are
```math
\{ (i + 1) * stride + i | i = 0, 1, \ldots, t - 1 \}.
```

To find the optimal value for the parameter $t$, we loop through values of $t$ starting from 1
and find the maximum value of $stride$ such that $s \geq 7$ (to ensure that each $s$-mer is likely to be unique).
Ties are broken by choosing the minimum $t$ to maximize $s$.
Therefore, $t$ is a function of $k$ and $d$.

The distance lower bound is
```math
stride + 1 = \left\lfloor \frac{t}{t + 1} \left( \frac{1}{d} - 1 \right) \right\rfloor + 1
```

As $k$ increases and $d$ is held constant, $t$ can be larger while still giving valid values of $s$.
Thus,
```math
\lim_{k \to \infty} \left\lfloor \frac{t}{t + 1} \left( \frac{1}{d} - 1 \right) \right\rfloor + 1 = \frac{1}{d}
```
In other words, the distance lower bound reaches the upper bound of $\frac{1}{d}$ as $k \to \infty$.
Therefore, DLB is optimal for large $k$.

## How to run
Install [Rust](https://www.rust-lang.org/tools/install), clone this repo, and then run this code using
```bash
cargo run --release --quiet > results.csv
python plot.py
```
