#!/bin/bash
set -e

echo "Step 1: Preprocessing"
bash 01_preprocess.sh

echo "Step 2: Alignment"
bash 02_align.sh

echo "Step 3: Polishing"
bash 03_polish.sh

echo "✅ Pipeline completed successfully"
