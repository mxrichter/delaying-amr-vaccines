# Evaluating the impact of malaria vaccination on drug-resistant <i>Plasmodium falciparum</i>: a modelling study

This code repository accompanies the manuscript "Evaluating the impact of malaria vaccination on drug-resistant <i>Plasmodium falciparum</i>: a modelling study" and can be used to reproduce simulations and publication figures.

Authors: Maximilian Richter<sup>1,2</sup>, Thiery Masserey<sup>1,2</sup>, Melissa A Penny<sup>3,4\*</sup>

<sup>1</sup> Swiss Tropical and Public Health Institute, Allschwil, Switzerland

<sup>2</sup> University of Basel, Basel, Switzerland

<sup>3</sup> The Kids Research Institute Australia, Nedlands, Western Australia, Australia

<sup>4</sup> University of Western Australia, Crawley, Western Australia, Australia

\* Corresponding author: Melissa A Penny (melissa.penny\@uwa.edu.au)

All simulations were performed using R 4.3.2 and OpenMalaria version 45 on sciCORE, the Center for Scientific Computing at the University of Basel, Switzerland (<https://scicore.unibas.ch/>).

## Overview

The scenario definitions including all parameters are defined in `config/scenarios.yaml`. The prefix of the scenario ID defines the type of each scenario:

-   *spread*: to estimate the impact of various model factors on the spread of a drug-resistant genotype (outcome metric: selection coefficient/rate of spread)

-   *morbidity*: to quantify the public health impact of vaccination in the context of drug resistance (outcome metric: reduction in cumulative uncomplicated cases in children younger than five years of age)

-   *treat*: to assess the impact of vaccination in the context of drug resistance on ACT treatment failure (outcome: reduction in proportion of failed treatments over 30-day horizon in children younger than five years of age)

Any standard options that are defined in `options.R` can be overwritten in `config/my_options.yaml`.

After specifying the scenario in `config/my_options.yaml`, sourcing the `launch.R` file in RStudio starts the simulation and analysis pipeline.

## Publication figures

All publication figures are saved in the folder `output/publication_figures/` in the working directory.

The following table outlines the process with which the publication figures can be reproduced:

| Figure     | Process (source file in RStudio)     | Required simulations |
| ---------- | ------------------------------------ | ---- |
| Figure 1   | `launch.R`                           | `spread_seasonal` |
| Figure 2   | `launch.R`                           | `morbidity_seasonal_eir`, <br /> `morbidity_seasonal_ve` |
| Figure 3   | `manuscript_treatment_failure.R`     | `treat_seasonal_eir`, <br /> `treat_seasonal_red`, <br /> `treat_seasonal_massvax_red`, <br /> `treat_seasonal_altvax_red`, <br /> `treat_seasonal_altvax_massvax_red`, <br /> `treat_perennial_eir`, <br /> `treat_perennial_red`, <br /> `treat_perennial_massvax_red`, <br /> `treat_perennial_altvax_red`, <br /> `treat_perennial_altvax_massvax_red` | 
| Figure S1  | `manuscript_further-si-plots.R`      | n/a |
| Figure S2  | `launch.R`                           | any `spread_seasonal_*` scenario |
| Figure S3  | `manuscript_further-si-plots.R`      | `spread_seasonal`, <br /> `spread_seasonal_massvax`, <br /> `spread_seasonal_altvax`, <br /> `spread_seasonal_altvax_massvax`, <br /> `spread_perennial`, <br /> `spread_perennial_massvax`, <br /> `spread_perennial_altvax`, <br /> `spread_perennial_altvax_massvax` |
| Figure S4  | `launch.R`                           | `spread_seasonal_altvax` |
| Figure S5  | `launch.R`                           | `spread_seasonal_massvax` |
| Figure S6  | `launch.R`                           | `spread_seasonal_altvax_massvax` |
| Figure S7  | `launch.R`                           | `spread_perennial` |
| Figure S8  | `launch.R`                           | `spread_perennial_altvax` |
| Figure S9  | `launch.R`                           | `spread_perennial_massvax` |
| Figure S10 | `launch.R`                           | `spread_perennial_altvax_massvax` |      
| Figure S11 | `manuscript_vaccine_effectiveness.R` | `morbidity_seasonal_eir`, <br /> `morbidity_seasonal_ve`, <br /> `morbidity_seasonal_massvax_ve`, <br /> `morbidity_seasonal_altvax_ve`, <br /> `morbidity_seasonal_altvax_massvax_ve`, <br /> `morbidity_perennial_eir`, <br /> `morbidity_perennial_ve`, <br /> `morbidity_perennial_massvax_ve`, <br /> `morbidity_perennial_altvax_ve`, <r /> `morbidity_perennial_altvax_massvax_ve` |
| Figure S12 | `manuscript_treatment_failure.R`     | `treat_seasonal_eir`, <br /> `treat_seasonal_red`, <br /> `treat_seasonal_massvax_red`, <br /> `treat_seasonal_altvax_red`, <br /> `treat_seasonal_altvax_massvax_red`, <br /> `treat_perennial_eir`, <br /> `treat_perennial_red`, <br /> `treat_perennial_massvax_red`, <br /> `treat_perennial_altvax_red`, <br /> `treat_perennial_altvax_massvax_red` |

*SI = Supplementary Information*
