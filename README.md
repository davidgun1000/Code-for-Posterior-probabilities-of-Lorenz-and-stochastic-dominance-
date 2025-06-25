# Code-for-Posterior-probabilities-of-Lorenz-and-stochastic-dominance-
This repository contains the MATLAB code for the paper:

"Posterior Probabilities of Lorenz and Stochastic Dominance for Australian Income Distributions"
David Gunawan, William Griffiths, and Duangkamon Chotikapanich (2021)
The Economic Record, 97(317), pp. 397–419.

🔧 Main Script
The main script is:

main_prog.m — Computes the posterior probabilities of Lorenz and stochastic dominance between two income distributions.

Please refer to the paper for full methodological details and interpretation of results.

📂 Input Requirements
The user must supply their own data in the form of two .mat files — one for each year (or group) being compared. Each file should contain:

A variable representing income

A variable representing the corresponding sampling weights

Example:
Example files included in this repository:

income_overall_2001_used.mat

income_overall_2010_used.mat

These files contain simulated data to illustrate the required format. You should replace them with your own income and weight data for actual analysis.

📝 Notes
The code assumes that both income and sampling weights are stored as column vectors.

The script outputs posterior probabilities for:

Lorenz dominance

First-order stochastic dominance

Second-order stochastic dominance


 

