## What is “strong evidence” in this bundle?

This `publication_ready/` snapshot is built to be **referee-reproducible** for the parts
that are actually deterministic and stable, and to clearly separate them from claims that
are currently *conjectural*.

### What is solid / reproducible

- **n=5 splitting field is Q(i)** (theorem-level statement; algebraic + reproducible checks).
  - Proof objects live in `paper/` and `referee_checks/verify_gram_levi.wl`,
    `referee_checks/verify_n5_discriminant.sage`.

- **n=7 inert-vanishing on dataset D₁** (computational proposition).
  - Run: `make verify-full` (or `make verify-sage`).
  - Driver: `referee_checks/verify_n7_theorem.sage`
  - Dataset: `data/D1_seeds.json`
  - Output: `logs/n7_full.json`

- **n=7 conjugation acts with cycle type 2¹²** at a fixed complex kinematic point (computational proposition).
  - Run: `julia referee_checks/verify_conjugates.jl`

- **n=7 needs more than Q(i)** (computational evidence).
  - Example: at inert primes we frequently see `N_p = 0` but `N_{p^2} > 0` with `N_{p^2} != 24`,
    which indicates the full splitting field is larger than Q(i).
  - Run: `sage referee_checks/verify_n7_np2_extension.sage --seed 0 --prime 7`

### What is conjectural (supported, but not proven here)

- **Conjectural candidate family**: a small imprimitive order-48 subgroup compatible with the 12 conjugate-pair block structure (dihedral-type families are consistent with this constraint).

The Julia monodromy scripts in `code/julia/` are retained as an *attempt* to extract candidate
rotation generators numerically, but this is not treated as “push-button proof” in the core
referee pipeline.

### Why we separate these layers

Numerical monodromy is sensitive to:
- path-tracking instabilities near discriminant loci
- cross-run solution labeling inconsistencies (can create spurious non-transitive “monodromy groups”)

The finite-field / Groebner / saturation checks are much more stable and are what the
`Makefile`-based reproduction targets.

