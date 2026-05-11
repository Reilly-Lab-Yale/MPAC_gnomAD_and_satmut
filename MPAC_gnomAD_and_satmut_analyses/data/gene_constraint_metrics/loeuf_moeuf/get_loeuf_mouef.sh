#!/bin/bash
wget -nc https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv
gzip gnomad.v4.1.constraint_metrics.tsv
