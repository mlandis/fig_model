The `anolis` directory contains the input files used for the empirical analysis.

One file describes the state space:
- `state_labels.csv` contains a matrix that defines how range-states (0, 1, 2, ..., 511) relate to presence-absence region-sets (100000000, 010000000, 001000000, ..., 111111111).

Three files describe the regional features:
- `anolis_nr9.area.txt` provides the relative sizes of regions
- `anolis_nr9.distance.txt` provides the relative (asymmetrical) distances from one region (row) to another (column)
- `anolis_nr9.category.txt` contains a matrix whose diagonal entries define which regions are continental (0) or insular (1), and whose off-diagonal entries define which regions have connections over land (0) or over water (1)

The remaining files describe biological datasets that were derived from the Poe et al. 2017 dataset. The Supplement describes how the Poe et al. data was converted into the nine region system. The first set includes the same original 383 species (379 *Anolis*, 4 outgroup) used in Poe et al 2017.
- `anolis.383spp.mcc.tre` is a time-calibrated *Anolis* phylogeny (time in units of Myr)
- `anolis_nr9_ns383.data.nex` is a nexus file containing presence/absence data; optional use in RevBayes
- `anolis_nr9_ns383.data.tsv` is a tab-delimited file containing integer-valued range data; used in RevBayes analyses
- `anolis_nr9_ns383.data_table.tsv` is tab-delimited file containing presence/absence data for species (rows) across regions (columns); mostly used to simplify reading

The second dataset is identical, but only includes the 379 *Anolis* species (no outgroups). We found that both datasets generate similar results.
