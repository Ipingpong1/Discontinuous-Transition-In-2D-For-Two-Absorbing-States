# Manuscript v3 — 4 November 2026

**Discontinuous transition into two symmetric absorbing states in two dimensions**

This revision implements the v2 referee report using the existing simulation archive. It contains a focused main paper and a separate supplement. No new production simulations were run. The v2 manuscript, analyses, and reports are preserved.

## Read and build

`main.pdf` is the main paper; `supplement.pdf` contains detailed protocols, time-window comparisons, secondary models, and numerical checks. Author names and affiliations remain to be supplied by the authors before submission. No public repository DOI or reviewer URL is asserted in this draft.

The standalone source ZIP contains `main.tex`, `supplement.tex`, the shared preamble, all included sections, bibliography, figure PDFs, and cross-document label caches. Extract it into a directory and run:

```sh
make TECTONIC=tectonic
```

Tectonic must be installed and able to access its TeX bundle, or have that bundle cached. The four build passes resolve references between the two documents. Keep the final PDFs together so cross-document hyperlinks work. The source ZIP does not contain the numerical input archive.

## Release packages

| Archive | Contents |
| --- | --- |
| `paper_v3_sources.zip` | Standalone typesetting sources and the exact publication figure PDFs |
| `paper_v3_reproducibility.zip` | CPU/CUDA source, parameter and seed manifests, archived reduction tables, revised analysis and figure scripts, manuscript, and verification records |
| `paper_v3_inputs.tar.zst` | Locally available trajectories, profiles, final configurations, and legacy figure inputs; about 13.8 GB before compression |

Extract the reproducibility ZIP and the input archive into the same project directory. Follow `REPRODUCE.md` at its root, also available as `campaign/v3/README.md`. Each archive has a file manifest with SHA-256 hashes; the release directory also provides hashes of the three archives.

`CHANGES.md` maps the referee's requests to the revision. `evidence_index.json` maps the figures and new analyses to project-relative scripts and results. Paths in that index describe the complete project/reproducibility package, not just the standalone typesetting ZIP.

## Scope of reproducibility

The v3 driver recalculates the new bulk flip-rate, front-estimator, droplet, endpoint, wall, and five-cell analyses, then regenerates all publication figures. Retained v2 analyses use explicitly archived reduction tables; their original reducer sources are included. The driver does not rerun the HPC campaign or independently repeat CUDA validation.

Some original snapshots underlying older roughness reductions are unavailable locally. The invalid free-band comparison has therefore been removed; the revised wall-only figure is rebuilt from the archived wall-edge reductions. Redundant restart copies, recovery/bundle mirrors, and inaccessible remote archives are not represented as part of the supplied input package.

The revision intentionally leaves unresolved the infinite-time voter endpoint, an asymptotic wetting exponent, and a possible much narrower intermediate regime. The five-cell front result is a one-sided extrapolation of coherent expansion, not a measured two-sided coexistence bracket.
