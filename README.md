# gbkscan2024

'gbkscan2024' is a bioinformatics tool designed to scan bacterial genomes for a user-specified predicted promoter site.
It is optimized for searching for the WhiG promoter in Actinomycetota genomes to predict a phylum-wide WhiG regulon.

---

## Workflow Overview

The core script, 'gbk_scan.pl', performs the following steps:
1. **Parses GenBank files** to locate coordinates of annotated genes in the file
2. **Extracts the upstream/promoter regions** (e.g., 200 bp upstream of the start codon) based on the gene's strand direction.
3. **Scans these regions** for your specified transcription factor binding site or promoter motif.
4. **Outputs a structured report** detailing the exact coordinates, distance to the start codon, and sequence of matching motifs.

## Quick Start (Docker)

To run the tool without installing any dependencies, 'Common.pm', you can run it directly via Docker. https://hub.docker.com/r/streptomyces/gbkscan2024
