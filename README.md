# Certified field model of elasmobranch electroreception

Code and data for

> An error-bounded field model of elasmobranch electroreception: canal funnelling, frequency
> tuning and the effect of body plan in a hammerhead shark and a manta ray. H. Talleb.

Every figure and every number in the article is produced by the scripts here. Nothing is
transcribed by hand: the tables of the supplementary information are generated directly from
the output of the sensitivity campaign.

## Reproduce everything

```
pip install -r requirements.txt
python3 reproduce_all.py
```

About two minutes. The run reports what succeeded, what was skipped and why, and writes
`reproduce_all.log`.

## What is where

| Path | Contents |
|---|---|
| `dgm/` | frozen subset of the field engine (5 modules), see `dgm/__init__.py` |
| `scripts/` | the model and figure scripts for the article |
| `sensitivity_lorenzini.py` | the sensitivity and robustness campaign |
| `make_sensitivity_figures.py` | figures S1 to S4 |
| `build_si.py`, `build_jtb.py` | the supplementary information and the manuscript |
| `reproduce_all.py` | entry point, runs all of the above in order |
| `meshes/` | the two surface models, see `meshes/README.md` |
| `LICENSE` | MIT for the code; the surface models keep their own CC BY-NC 4.0 |

## The method, in one paragraph

The organ is treated as an electroquasistatic conduction problem. The discrete geometric
method writes the field laws on a primal Delaunay mesh and its circumcentric dual, so the
topological operators are exact and only the constitutive (Hodge) operator is approximate. A
primal solve in the potential over-estimates the dissipated energy; a complementary,
divergence-conforming RT0 flux solve under-estimates it. The exact value lies between them.
This is what lets a comparison be settled outright: when the bracket computed with the gel
canal and the bracket computed without it are disjoint, the ordering holds for the exact
solution and not merely for the computed one.

## Surface models

The great hammerhead (*Sphyrna mokarran*) and manta ray (*Mobula birostris*) surface models
are by DigitalLife3D (digitallife3d.org) under CC BY-NC 4.0. See `meshes/README.md`. The two
scripts that need them are skipped with an explicit message when they are absent, so the rest
of the pipeline still runs.

## Licence

The code is under the MIT licence, so anyone reviewing or extending the article can run,
modify and redistribute it without asking. The two surface models under `meshes/` are not
mine to relicense: they stay under CC BY-NC 4.0, their author's terms. See `LICENSE`.

## On tool use

A generative artificial-intelligence assistant was used while developing this code and
drafting the article, and this is declared in the manuscript. No figure, number or table is
produced by that assistant: each is the output of the deterministic, seeded scripts in this
repository, which is why the repository exists.
