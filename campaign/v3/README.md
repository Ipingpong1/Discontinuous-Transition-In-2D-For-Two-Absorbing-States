# Reproducing the v3 paper and supplement

This release regenerates the revised numerical analyses, every figure used in
the v3 paper and supplement, and both linked PDFs from the supplied project.
The procedure uses archived simulations. It launches no new simulations,
submits no cluster jobs, and requires no GPU or simulator executable.

## Release layers

| Archive | Contents and use |
| --- | --- |
| `paper_v3_sources.zip` | Standalone typesetting sources, bibliography, exact figure PDFs, cross-reference auxiliary files, and a Makefile. Extract into its own directory for typesetting only. |
| `paper_v3_reproducibility.zip` | Project-relative CPU/CUDA source, manifests, analysis and figure scripts, fixed earlier reduction tables, v3 results, manuscript sources, and verification records. This README also appears as `REPRODUCE.md` at its root. |
| `paper_v3_inputs.tar.zst` | Companion input trees: locally available trajectories, profiles, snapshots, dense final configurations, and retained legacy reductions. Extract into the **same project root** as the reproducibility ZIP. |

The full regeneration route needs the reproducibility ZIP and input archive.
The source ZIP is a separate typesetting route; do not unpack its flat contents
over the project root. Input paths and code paths are complementary, not
alternative copies of a complete project. The input manifest reports about
13.8 GB before compression; allow additional space for the archives, code, and
generated outputs. Exact file sizes are in the supplied manifests.

From a directory containing the release archives and their external
`SHA256SUMS`, verify and extract them:

```sh
sha256sum -c SHA256SUMS
mkdir dn-v3
unzip paper_v3_reproducibility.zip -d dn-v3
tar --zstd -xf paper_v3_inputs.tar.zst -C dn-v3
cd dn-v3
sha256sum -c SHA256SUMS
```

The first checksum list covers the archives. The second is the list inside the
reproducibility ZIP and covers its extracted files. Run these checks before
regeneration, which intentionally overwrites v3 results, figures, and logs.
On systems without `tar --zstd`, use
`zstd -dc paper_v3_inputs.tar.zst | tar -xf - -C dn-v3` from the download
directory. The source ZIP has its own `SHA256SUMS` and `SOURCE_MANIFEST.json`.

To check every extracted input against its recorded size and SHA-256 digest,
run this from the project root:

```sh
python -B - <<'PY'
from pathlib import Path
import hashlib
import json

manifest = json.loads(Path('releases/v3/INPUTS_MANIFEST.json').read_text())
for item in manifest['files']:
    path = Path(item['path'])
    if path.stat().st_size != item['bytes']:
        raise SystemExit(f'Size mismatch: {path}')
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    if digest.hexdigest() != item['sha256']:
        raise SystemExit(f'Checksum mismatch: {path}')
print(f"Verified {len(manifest['files'])} input files")
PY
```

## Dependencies and commands

The full analysis/figure/test run was verified in the local `SSU` environment
with Python 3.12.13, NumPy 2.4.3, SciPy 1.18.1, Matplotlib 3.11.2, and Pillow
12.1.1. These are recorded working versions, not tested minimum versions.
With Python 3.12 available, an isolated environment can be prepared using:

```sh
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install numpy==2.4.3 scipy==1.18.1 matplotlib==3.11.2 pillow==12.1.1
```

Install Tectonic separately and put it on `PATH`; the verified build used
Tectonic 0.17.0. Its first build may download TeX packages. A populated Tectonic
cache is needed for an offline build. Archive extraction also requires `unzip`,
`tar`, and `zstd`; `sha256sum` performs the shell checksum checks above.
Matplotlib uses a noninteractive backend, and figure rendering does not need
an external LaTeX installation.

Run the complete regeneration from the extracted project root:

```sh
python -B -m campaign.v3.reproduce --tests
python -B -m campaign.v3.build_paper --engine tectonic
```

The first command runs four v3 analysis modules, 18 figure scripts, and the
existing analysis unit-test suite. The archived verification run passed all
51 tests. Some scripts produce multiple figures. The second command builds
`main`, `supplement`, `main`, and `supplement` in that order to resolve the
namespaced cross-document references. It fails if either final TeX log
contains undefined references, multiply defined labels, or overflowing boxes.
Pass `--engine /path/to/tectonic` if the executable is outside `PATH`, and
optionally `--only-cached` to prohibit package downloads.

The final documents are `paper/v3/main.pdf` and `paper/v3/supplement.pdf`.
Numerical outputs are under `campaign/v3/results/`; figure PDFs and PNGs are
under `paper/v3/figs/`. The run writes `reports/v3/reproduce.log`,
`reports/v3/reproduce_verification.json`, and `reports/v3/paper_build.log`.
The reproduction driver hashes existing files in `paper/v2`, `campaign/v2`,
and `reports/v2` before and after execution and fails if any change. It uses
the interpreter that invoked it, so no local environment name or absolute
author path is required.

For rendering from the supplied reductions without rerunning the four v3
modules, use `python -B -m campaign.v3.reproduce --figures-only`. This still
needs the companion inputs because several figure scripts read trajectories
or legacy reductions directly. To typeset only from the separate source ZIP,
extract it into an empty directory, verify its `SHA256SUMS`, and run
`make TECTONIC=tectonic` there.

## What is recomputed

The four revised modules deliberately preserve the earlier analysis trees.
They use existing v2 readers and statistical helpers where appropriate.
Each writes its own tables and input/code provenance under
`campaign/v3/results/<module>/`.

| Module | Recomputed quantities and input boundary |
| --- | --- |
| `bulk_activity` | Reads the 12 largest flagship bulk time series; computes time-integrated bond density, actual flip fraction, and late-window density changes. |
| `fronts_review` | Reads declared completed GPU time series/profiles and all 32 large disk/hole trajectories. Recomputes front proxies with a common whole-cohort geometry gate, paired/window and direction contrasts, bulk means, alternative droplet areas, and the integrated curvature prediction. Retains v2 short-survey front points as the velocity input to that prediction, and verifies unchanged primary measurements against v2 tables. |
| `film_endpoint` | Reads the four declared 20-million-step endpoint cohorts from raw time series, extending absorbed runs through their known absorbing state. Recomputes unconditional means and paired changes. For wall roughness it averages the **archived wall-edge reductions**, computes finite-distance slopes, and combines thickness histograms; it does not regenerate those reductions from original snapshots. |
| `molofsky` | Reads 64 five-cell strip profiles and their metadata, four lifetime cohorts, and recorded voter-lineage histograms. Recomputes width/fit-window sensitivity, interior-decay diagnostics, the one-sided displayed onset extrapolations, all-run restricted mean absorption times, and conditional lineage fractions. |

The figure scripts also regenerate revised displays and some reductions, such
as the probability/density survey and droplet bootstrap summaries, from their
declared inputs. The main and supplement evidence index is
`paper/v3/evidence_index.json`; figure-specific provenance and method files
identify the actual inputs used by the revised scripts.

## Retained results and limits

This pipeline is not a fresh rerun of every v2 reducer or the original
Monte Carlo campaigns. The archived `campaign/v2/results/`,
`campaign/v2/theory/results/`, `campaign/theory/`, and `data/figdata/`
tables remain explicit inputs where the corresponding analyses are retained.
Examples include the existing spatial-correlation reduction, survival and
wetting summaries, closure threshold/matching tables, and early parameter
surveys. Their available code, manifests, and input files are included for
inspection. The v3 theory figure redraws recorded cutoff results; it neither
generates new coalescing walks nor extrapolates beyond the recorded cutoffs.

The original free-band snapshots needed to repair the old roughness comparator
are unavailable locally. That comparator has been removed. The retained wall
plot uses the archived reductions of the two physical film edges. Some legacy
inputs are therefore reductions rather than original configurations; the
input manifest states this limitation. The data archive retains dense final
states but omits redundant restart checkpoint copies, recovery snapshots,
bundle mirrors, and unavailable remote archives. Scientific run manifests
and explicit completion checks determine eligibility; the presence of a
legacy file alone does not make it an eligible measurement.

Checksum manifests verify the distributed bytes. Analysis provenance records
which bytes were used. Regenerated PDF bytes need not match the distributed
PDF digests because timestamps, fonts, and rendering versions can change;
compare numerical tables and scientific content when assessing a rebuild.
The recorded bootstrap and Student-t calculations retain their stated
finite-sample assumptions; successful reproduction does not turn estimator
or finite-window sensitivity into a confidence interval for an asymptotic
phase boundary.

Release packaging is defined by `campaign/v3/package_release.py`: its default
mode creates the two ZIPs and archive checksum list, while `--data` creates
the compressed input archive and `INPUTS_MANIFEST.json`. Packaging does not
run simulations or submit jobs.
