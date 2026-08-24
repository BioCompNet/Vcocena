# vcocena in Docker

Reproducible container for the Vertical CoCena² pipeline.

## Build

From the **repository root** (the directory containing `DESCRIPTION`), not from
`docker/`:

```bash
docker build -f docker/Dockerfile -t vcocena:1.0.0 .
```

The build is layered so that editing package source under `Vcocena/R/` only
re-runs the fast package install, not the ~30-minute dependency install.
`.dockerignore` keeps the build context at ~1 MB by excluding the 2.1 GB of
`.RData` workspaces and the datasets.

## Run

**Smoke test** (self-contained, seconds — uses the 100-gene example baked into
the package):

```bash
docker run --rm vcocena:1.0.0 Rscript /opt/docker/run_smoke.R
```

**Full run** on a real two-layer dataset:

```bash
mkdir -p out_docker
docker run --rm \
  -v "/path/to/your/dataset:/data:ro" \
  -v "$PWD/out_docker:/out" \
  vcocena:1.0.0 Rscript /opt/docker/run_full.R
```

Expect 10–30 min and **at least 8 GB** of container memory: `correlation_actions()`
materializes a data frame of all `choose(5000,2)` ≈ 12.5M gene pairs per layer.

The dataset mount must contain `data/` (count matrices + `reference_files/`) and
`sample_info/` (annotations). Paths are overridable with `VCOCENA_DATA` and
`VCOCENA_OUT`.

## Why the inputs are copied out of the mount

`run_vcocena()` creates its `save_folder` **inside** `working_directory` and
writes plots there, so the analysis directory can never be a read-only mount or
the installed package library. Both run scripts copy the (small) inputs into
`/work` first. The same applies to `vcocena_example()`, which returns a path
inside the installed library.

## What the full run emits

Alongside the pipeline's PDFs, `run_full.R` writes three TSVs to the output
mount:

| file | contents |
|---|---|
| `integrated_cluster_information.tsv` | integrated module assignments |
| `integrated_GFC_per_cluster.tsv` | Group Fold Change per module per condition |
| `integrated_edgelist.tsv` | the integrated network |

The pipeline has no tabular export path of its own — every native output is a
PDF — so these exist to give a future regression baseline something numeric to
compare against.

## Notes

- Base image `bioconductor/bioconductor_docker:RELEASE_3_22-R-4.5.2` (R 4.5.2),
  chosen to match the R 4.5.0 that produced the existing `out_stepwise/` output.
- `MCDA` is installed from the CRAN Archive: it was removed from CRAN on
  2025-11-20 but is still used on the main path.
- Only hard dependencies are installed. `devtools` is **not** present, so the
  `vCoCena_stepwise_*.R` drivers will not run in this image — they also depend on
  `plot_sample_distributions()`, which was never ported into the package. Use
  `run_vcocena()`.

## Reproduction harness

`run_paper.R` and `run_paper_integration.R` compare this package against the
published vCoCena analysis in Carraro et al., eLife 2022
([doi:10.7554/eLife.78012](https://doi.org/10.7554/eLife.78012)), using the
authors' own saved R environment as ground truth (Zenodo
[10.5281/zenodo.6984701](https://doi.org/10.5281/zenodo.6984701), CC-BY-4.0).

Two findings from that exercise are worth knowing before using either script:

- **`run_vcocena()` cannot reproduce that paper end to end.** Its published
  layers were hCoCena integrations across two cell lines, not count matrices.
  `run_paper_integration.R` therefore injects the published layer networks
  directly into the vertical-integration stage, which reproduces the published
  network exactly (3,405 nodes / 186,723 edges / 35,453 cross-layer edges).
- **Module counts will not match a pre-2024 analysis.** igraph 2.0 rewrote
  community detection; on the identical graph, `cluster_louvain` returns 6
  communities where igraph 1.2.6 returned 147. This affects any re-run of an
  older CoCena analysis, not just this one.
