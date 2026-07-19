# The Collatz Sessions, part 1 — The experiment: a human, an AI, and the most dangerous problem in mathematics

*This is the first post in a six-part series about a months-long research collaboration between me (Martien de Jong) and an AI research agent, attacking the Collatz conjecture. Everything described here is documented, with code, proofs and data, in the public repository [github.com/martiendejong/collatz_ai](https://github.com/martiendejong/collatz_ai).*

## The problem

Take any whole number. If it is even, halve it. If it is odd, triple it and add one. Repeat.

7 → 22 → 11 → 34 → 17 → 52 → 26 → 13 → 40 → 20 → 10 → 5 → 16 → 8 → 4 → 2 → 1.

Every number ever tried — and computers have tried every number up to 295 quintillion — ends at 1. The **Collatz conjecture** says every number does. It has been open for ninety years. Paul Erdős said mathematics is not ready for such problems, and famously offered $500 for a solution — a suspiciously low bounty for a problem this famous, which tells you what he thought the hourly rate would work out to.

## The experiment

I wanted to know what happens when you point a modern AI research agent at this problem, full time, with three rules:

1. **Everything gets verified.** Every claim is tested numerically before it is believed, ideally against tens of thousands of cases.
2. **Everything gets labeled honestly.** Proved, verified, measured, conjectured, or open — and when something we found turns out to already exist in the literature, we say so, by name.
3. **Everything gets published.** Every session, every dead end, every same-day correction of our own mistakes is in the repo. About 3,300 numbered research rounds so far.

I supplied direction and intuition — I have no formal mathematical training beyond an engineering background, but I had a framework of ideas about how the problem is organized (what I called *families* and *sequences* of numbers, and a suspicion that the whole problem is really about the conversion between binary and ternary notation). The AI supplied formalization, computation, proof-writing, and — importantly — pushback when my intuitions were wrong, which happened regularly and is documented just as carefully as when they were right.

## What came out

Over the series you will see the whole arc, but here is the honest headline version:

- **We independently rebuilt most of ninety years of Collatz understanding** — often without knowing it existed, discovering afterwards that our "new" ideas dated from 1976, 1987, 2012 or 2020. That was humbling and, in a strange way, validating: the framework I started from kept turning out to be the one the professionals use.
- **We produced three things that appear to be genuinely new** (after a targeted literature check): a computational record on the best known density-bound method (an exponent of 0.9146 where the published record is 0.84), a new analytical lens on that same method, and a new open question about where that method ultimately leads. Part 6 covers these with all the caveats.
- **We did not solve the Collatz conjecture.** Nobody has. But we can now tell you, with unusual precision, *why* it is hard — and that story, which unfolds over parts 2 to 5, is the real payoff: a tiny machine, a fuel economy, a fair coin, and a calendar of near-misses, all provable, surrounding three walls that are each one sentence long.

## Why publish this

Partly because the material is beautiful and I want people to see it. Partly as a case study in human-AI research: what such a collaboration can and cannot do, and what intellectual honesty has to look like when an AI can produce plausible-sounding mathematics at industrial speed. Every "theorem" in this series links to a numbered entry in the repository with its proof or its verification code. Where we merely rediscovered, the original authors are named.

*Next: [Part 2 — The machine: Collatz is a two-rule automaton, and you can watch it burn](02-the-machine.md).*
