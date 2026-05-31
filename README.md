basenji_cre_state_hmm.py

# Basenji motif state annotation using HMM

This repository contains a Python script for annotating motif-level regulatory states
based on Basenji-predicted scores using a two-stage Hidden Markov Model (HMM) approach.

The script is designed for downstream analysis of chromatin accessibility–related
sequence effects and was developed for cross-species regulatory element comparison.

---

## Overview

The script performs the following steps:

1. Load base-resolution Basenji prediction scores stored in HDF5 format
2. Compute mean sequence effect scores across nucleotide substitutions
3. Standardize (z-score) the mean scores
4. Apply a Gaussian HMM to estimate posterior state probabilities
5. Use these probabilities as emissions in a Multinomial HMM
6. Decode the final state sequence using the Viterbi algorithm
7. Output per-base regulatory states as a CSV file

---

## Requirements

The script was tested with Python 3 and requires the following packages:

- numpy
- pandas
- h5py
- scikit-learn
- scipy
- hmmlearn
- matplotlib
- seaborn

## How were the HMM transition probabilities chosen?

Basenji-HMM uses a three-state Hidden Markov Model (HMM) to segment base-resolution ISM scores into:

- background
- CA-increasing elements
- CA-decreasing elements

The transition probabilities were chosen to encode a simple biological prior, rather than being arbitrary numerical settings. We expect functional sequence elements to occur as relatively short contiguous segments, whereas most bases in the genome should belong to longer background intervals.

For an HMM, if `a_ii` denotes the self-transition probability of state `i`, then the expected run length of that state is approximately:

`E[L_i] ≈ 1 / (1 - a_ii)`

Under the current parameterization:

- background state: `E[L_0] ≈ 1 / (1 - 0.98) = 50 bp`
- CA-increasing state: `E[L_1] ≈ 1 / (1 - 0.84) = 6.25 bp`
- CA-decreasing state: `E[L_2] ≈ 1 / (1 - 0.84) = 6.25 bp`

This makes the model behavior biologically interpretable. The two functional states are expected to span only a few base pairs on average, which matches the typical scale of core cis-regulatory motifs, often around 6-20 bp. By contrast, the background state is expected to persist much longer, reflecting the relative sparsity of functional bases across the genome.

In other words, we set the expected length of functional states to about 6 bp, so that the HMM preferentially identifies short high-impact segments, while the background state is allowed to run longer to represent non-functional sequence.

The initial state probabilities follow the same logic. We assign a high prior probability to the background state (`0.95`) and much smaller priors to the two functional states (`0.025` each), reflecting the expectation that most genomic positions belong to background and only a minority belong to high-impact regulatory segments.

This parameterization is intended to balance interpretability and biological realism:
- short functional runs consistent with motif-scale regulatory elements
- long background runs consistent with sparse regulatory architecture
- reduced tendency to produce unstable base-by-base switching between opposite functional states
