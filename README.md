# Living with IBD

## Overview

This dataset consists of stool samples collected longitudinally from patients with inflammatory bowel disease (IBD), including both Crohn’s disease and Ulcerative Colitis. All samples were processed for shotgun metagenomic sequencing to characterize taxonomic and functional dynamics of the gut microbiome across different disease states and time points.

The study design includes up to three sampling points per participant over the course of one year, providing a robust framework for longitudinal microbiome analysis in IBD.

## Cohort Structure

Participants: 121 individuals (patient IDs P101–P321)

Diseases represented:
Crohn’s disease
Ulcerative colitis

Sample type: Stool

Assay: Shotgun metagenomic sequencing

Clinical metadata collected per sample:
sex (Male/Female)
severity: remission or flare
diagnosis: Crohns_disease or Ulcerative_Colitis
timePoint: wk0, wk26, wk52
patient ID: unique identifier (P###)

## Sampling Time Points

Participants were sampled at:

|Time Point          |Label | Description                    |
|--------------------|------|--------------------------------|
| Baseline (week 0)  | wk0  | Initial sample upon enrollment |
| Mid-study (week 26)| wk26 | 6-month follow-up              |
| One-year (week 52) | wk52 | 12-month follow-up             |

Not all participants have samples at all three time points. Many have paired baseline–wk52 samples, while a subset includes the intermediate wk26 collection.

## Sample Naming Scheme

Two naming patterns occur:

X###_Baseline / X###_WK52
e.g., X101_Baseline, X101_WK52
The prepending "X" was added because R does not like variables that strts with a number. It will preprend the "X" anyways.

W26_P### (wk26 samples)
e.g., W26_P101 for week-26 sample from patient P101

Every row in the metadata table corresponds to a single stool sample with complete clinical annotations.

## Longitudinal Clinical Dynamics

The dataset captures transitions in disease activity over time, such as:

flare → remission
remission → flare
persistent remission
persistent flare

This enables modeling of:

intra-patient microbiome changes
disease-specific microbial signatures
flare-associated functional shifts
baseline predictors of future disease state
sex- or diagnosis-dependent trajectories

Please click on the following ling to see a descriptive list of files included in this [repo](initial_analysis/description_of_files.md)


