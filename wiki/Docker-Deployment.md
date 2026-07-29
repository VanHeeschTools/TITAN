# Docker Deployment

TITAN ships as a single Docker image that contains the R runtime and all R packages. Study data and reference files are **never baked into the image** — they are provided via volume mounts at runtime.

---

## Prerequisites

- Docker Desktop (Mac/Windows) or Docker Engine (Linux)
- On Apple Silicon Macs: `--platform linux/amd64` flag required (the base image is amd64-only)
- For local testing: the reference files and at least one study RDS on your local machine

---

## Repository Layout (build-time)

The Dockerfile copies only `app/` into the image. The `.dockerignore` excludes:

| Excluded path | Reason |
|---|---|
| `app/ref/` | Large reference databases; provided by volume mount |
| `app/data/` _(except `.gitkeep`)_ | Study RDS files and catalog; provided by volume mount |
| `app/bin/ncbi-blast-2.14.0+/` | Extracted BLAST binaries; not needed (blastp.REAL is already there) |
| `app/bin/ncbi-blast-2.14.0+-x64-linux.tar.gz` | 243 MB tarball; not needed |
| `.git/`, `design/`, `documentation/`, `scripts/` | Not needed inside the container |

This keeps the build context small and means the image itself contains **no patient data or study-specific files**.

---

## Dockerfile

```dockerfile
FROM rocker/shiny:4.4.2

# System libraries for compiled R packages
RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# remotes is the pinning mechanism; install it first
RUN Rscript -e "install.packages('remotes', repos='https://cloud.r-project.org')"

# Pinned CRAN packages (match local environment)
RUN Rscript -e "\
  r  <- 'https://cloud.r-project.org'; \
  iv <- function(p, v) remotes::install_version(p, version=v, repos=r, upgrade='never'); \
  iv('shiny',        '1.14.0'); \
  iv('bslib',        '0.11.0'); \
  iv('DT',           '0.34.0'); \
  iv('plotly',       '4.12.0'); \
  iv('dplyr',        '1.2.1');  \
  iv('tidyr',        '1.3.2');  \
  iv('stringr',      '1.6.0');  \
  iv('ggplot2',      '4.0.3');  \
  iv('shinyWidgets', '0.9.1')"

# Unpinned utilities
RUN Rscript -e "install.packages(
  c('cowplot','data.table','jsonlite','yaml','processx','XML','BiocManager'),
  repos='https://cloud.r-project.org')"

# Bioconductor packages — rBLAST is only on Bioconductor, not CRAN
RUN Rscript -e "BiocManager::install(
  c('Biostrings', 'IRanges', 'rBLAST'), ask=FALSE, update=FALSE)"

# Copy app source (data/ and ref/ excluded by .dockerignore)
COPY app/ /srv/titan/

# Create the blastp symlink expected by rBLAST (only blastp.REAL is on disk)
RUN ln -sf /srv/titan/bin/blastp.REAL /srv/titan/bin/blastp \
 && chown -R shiny:shiny /srv/titan

USER shiny

# Cloud Run injects $PORT (default 8080); local docker run falls back to 3838.
CMD ["Rscript", "-e",
     "shiny::runApp('/srv/titan', host='0.0.0.0', port=as.integer(Sys.getenv('PORT', 3838)))"]
```

---

## Building the Image

```bash
# Standard (Linux x86-64)
docker build -t titan:dev .

# Apple Silicon Mac
docker build --platform linux/amd64 -t titan:dev .
```

The build takes 10–20 minutes on first run (Bioconductor packages are slow to install). Subsequent builds use Docker layer cache for all package layers as long as the Dockerfile doesn't change above the `COPY` instruction.

To verify the rBLAST installation before running:
```bash
docker run --rm [--platform linux/amd64] titan:dev \
  Rscript -e "packageVersion('rBLAST')"
```

---

## Running Locally

### Minimum — app only (no data)

```bash
docker run --rm -p 3838:3838 [--platform linux/amd64] titan:dev
```

The app starts but shows an empty Study Library and the cross-reactivity / BLAST checks will print "not loaded" messages. Useful to confirm the image builds and starts.

### With reference data and study data

```bash
docker run --rm \
  [--platform linux/amd64] \
  -v ~/titan-test/ref:/srv/titan/ref:ro \
  -v ~/titan-test/data:/srv/titan/data:ro \
  -p 3838:3838 \
  titan:dev
```

Open http://localhost:3838 in your browser.

**`~/titan-test/` layout** (sync from HPC first — see below):
```
~/titan-test/
├── ref/
│   └── ensembl114_pep/
│       ├── ensembl114_pep_index.rds
│       ├── ensembl_gene_annotation.rds
│       └── ensembl114_pep.*        (BLAST database files)
└── data/
    ├── catalog.yaml
    └── rms_organoids/
        ├── titan_rms_organoids.rds
        └── peptides_rms_organoids.csv
```

### Syncing data from HPC for local testing

```bash
HPC="user@hpc.example.org"
TITAN_HPC="/hpc/pmc_oatv/projects/tools_dev/titan/app"

# Reference data
rsync -av --progress \
  "${HPC}:${TITAN_HPC}/ref/ensembl114_pep/" \
  ~/titan-test/ref/ensembl114_pep/ \
  --exclude "*.fa" --exclude "*.fa.gz"   # skip raw FASTAs; not needed at runtime

# Catalog and study data (exclude raw RDS if too large; sync just what you need)
rsync -av "${HPC}:${TITAN_HPC}/data/catalog.yaml" ~/titan-test/data/
rsync -av "${HPC}:${TITAN_HPC}/data/rms_organoids/" ~/titan-test/data/rms_organoids/
```

---

## Volume Mount Reference

| Host path | Container path | Contents | Access |
|---|---|---|---|
| `~/titan-test/ref` | `/srv/titan/ref` | Ensembl pep index + BLAST db + GENCODE CSV | `:ro` |
| `~/titan-test/data` | `/srv/titan/data` | `catalog.yaml` + per-study RDS + peptide files | `:ro` |

Both mounts are read-only (`:ro`). The app never writes to these directories.

---

## Cloud Run Deployment (Overview)

> Full Cloud Run setup instructions are pending IT confirmation of the VPN/Interconnect network restriction approach. This section provides the planned architecture.

### Architecture

```
Browser (VPN-restricted) → Cloud Run service → GCS volumes
                                              ├── gs://titan-ref/
                                              └── gs://titan-data/
```

- The **Docker image** is pushed to Google Artifact Registry (private).
- **Reference data** lives in a private GCS bucket (`gs://titan-ref/`) mounted as `/srv/titan/ref`.
- **Study data** lives in a private GCS bucket (`gs://titan-data/`) mounted as `/srv/titan/data`.
- **Network access** is restricted to institutional VPN users only (IAP or VPC connector — TBD with IT).
- The **Cloud Run service** runs as a low-privilege service account with read-only access to the buckets.

### Planned deploy command

```bash
gcloud run deploy titan \
  --image REGION-docker.pkg.dev/PROJECT/titan/titan:VERSION \
  --region REGION \
  --port 8080 \
  --add-volume name=ref,type=cloud-storage,bucket=titan-ref \
  --add-volume-mount volume=ref,mount-path=/srv/titan/ref \
  --add-volume name=data,type=cloud-storage,bucket=titan-data \
  --add-volume-mount volume=data,mount-path=/srv/titan/data \
  --service-account titan-runner@PROJECT.iam.gserviceaccount.com \
  --no-allow-unauthenticated
```

### Security requirements

- The image must be in a **private** Artifact Registry repository (not Docker Hub).
- The GCS buckets must be **private** (no public access, no allUsers IAM binding).
- Study data (patient/cohort-derived RDS objects) must **never** appear in a public registry, public bucket, or public git history.

---

## Environment Variables

| Variable | Default | Description |
|---|---|---|
| `PORT` | `3838` | Port the Shiny app listens on. Cloud Run injects `8080`; local runs use the default. |

---

## Updating the Image

When app code changes (not reference data or study data):

```bash
# Rebuild
docker build [--platform linux/amd64] -t titan:dev .

# Tag for registry
docker tag titan:dev REGION-docker.pkg.dev/PROJECT/titan/titan:VERSION

# Push
docker push REGION-docker.pkg.dev/PROJECT/titan/titan:VERSION

# Redeploy (Cloud Run picks up the new image on next request by default)
gcloud run services update titan \
  --image REGION-docker.pkg.dev/PROJECT/titan/titan:VERSION \
  --region REGION
```

Reference data and study data updates do **not** require a rebuild — just update the GCS bucket or re-sync the local volume.

---

## Troubleshooting

### App crashes at startup with `Error in packageVersion(pkg)`
A package listed in the About tab sidebar is not installed. The `tryCatch` wrapper in `app.R` should prevent this from crashing the app — if it does crash, rebuild the image to pick up the latest Dockerfile.

### `blastp not found in PATH`
The `blastp` symlink at `/srv/titan/bin/blastp` is missing. Ensure the `RUN ln -sf …` step in the Dockerfile ran successfully. Verify:
```bash
docker run --rm titan:dev ls -la /srv/titan/bin/blastp
```

### `BLAST database not found: ref/ensembl114_pep/ensembl114_pep`
The `ref/` volume mount is missing or the BLAST files were not synced. Check that the `.phr`, `.pin`, and `.pdb` files are present in the mounted directory.

### `Ensembl 114 pep index not loaded`
`ref/ensembl114_pep/ensembl114_pep_index.rds` is absent. Re-run script 01 and re-sync.

### `No catalog found at data/catalog.yaml`
The `data/` volume mount is missing or `catalog.yaml` was not synced. The Study Library will be empty but the app will still start.
