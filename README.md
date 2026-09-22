# Short terminal-path closures — verification code

Companion to “Short terminal-path closures and local distance domination in
graphs excluding K_{2,t} as a minor”.

Authors: **Marek Gałązka** and **Hanna Wdowicka**.
The companion manuscript is being prepared for submission to *Theoretical
Computer Science*; it has not been published or accepted.

The results are theoretical. These standard-library Python scripts provide
finite checks of the proofs and constructions; they do not establish the
general theorems and are not a distributed-system implementation.

## Requirements and execution

Python 3.10 or newer, with no third-party packages. Keep all scripts in the same
directory because the colouring verifier imports the independent minor oracle
and the graph/closure constructors from the earlier verifiers.

```bash
python3 verify_boundary_bounds.py
python3 verify_extension.py
python3 verify_colouring.py
python3 verify_bouquet_lower.py
```

Each script writes its associated JSON file. Deterministic random seeds are
recorded in the scripts and reports. Elapsed times depend on the machine.

## Scope

| Script | Object checked | Recorded coverage |
|---|---|---|
| `verify_boundary_bounds.py` | Pairwise interface bounds | 207,387 connected-set pairs; 622,161 inequalities |
| `verify_extension.py` | Product recurrence, global boundaries, explicit constructions | 131,068 closure inequalities; 21,930 boundary inequalities; 480 Voronoi instances; 12 ladder chains; 6 iteration examples |
| `verify_colouring.py` | Monochromatic-component contractions, witness probabilities and exponential closure bound | 35,634 exhaustive and 2,880 sampled colouring outcomes; 34 exact rational probability equalities; 8 nonvacuous constructed bounds |
| `verify_bouquet_lower.py` | Actual DomSet lower-bound construction | 128 parameter instances; direct balls, unique maximum choices, dominating sets and packing certificates; largest graph 9,911 vertices |

The minor oracle enumerates connected partitions, contracts them, and counts
common neighbours. It does not assume the leaf-contraction or random-colouring
proof. The probability checks use exact rational arithmetic. The lower-bound
verifier checks strict uniqueness of each forced maximum, so the checked
witnesses do not depend on identifier tie-breaking. The test files record no
counterexample or failure.

There are no external empirical datasets. ChatGPT/Codex (OpenAI, September
2026) assisted in developing arguments, counterexample searches and the code.
The mathematical claims and the scope of the finite checks should be assessed
from the proofs and executable sources, respectively.

## Reproducibility

The JSON files record the checks accompanying the manuscript. Running the
scripts replaces these files with new results; elapsed-time fields may differ.
When citing this code in a manuscript, identify the exact Git commit.
