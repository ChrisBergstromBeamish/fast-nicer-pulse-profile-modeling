# Roadmap

One sprint at a time. A sprint starts only when the one before it is accepted.

1. **The repository stands on its own.** A clean copy sets up by the setup page alone; one command repeats every
   check for the Python solution and for the compiled solution; the documentation shows the work step by step, so a
   reader can check each stated number.
2. **The compiled solution in parallel on a processor's cores.** Many stars worked out at once (a sampler asks for
   many at a time). Measured on the laptop's cores against the one-core time, with the same numbers as the Python
   solution. Done (sprint 2, 23 November 2026): on eight threads (lines of work of one program, each on a core), 1,711
   stars per second at the working setting against 363 on one thread (4.7 times; 4,096 stars from the public example's
   prior, mains power; `CLAIMS.md`, row I1); in the checks, each of 108 stars the same, bit for bit, in batches on 1 to 8
   threads as alone (row J1).
3. **The compiled solution on a graphics card.** The target is the graphics hardware of the author's own laptops (the one
   measured so far has built-in Intel Arc 140V graphics and no NVIDIA card), written so that others can port it to their
   own cards.

Later, in an order still to be decided: the public code's time measured on the same machine (done on 23 November 2026:
`PUBLIC_CODES.md`, section 1.3); a tabulated hot-gas (atmosphere) spectrum in place of the blackbody; the second path
of light around very compact stars; rings visible in two pieces.
