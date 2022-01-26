The `inference` directory contains four files for running a FIG analysis:
- the `run_job.Rev` script contains all variables that most users will want to modify
- the `mcmc_FIG.Rev` script performs the MCMC analysis by loading the phylogegenetic model (`model_FIG.Rev`) and geographical model (`geo_FIG.Rev`) for FIG, creates the appropriate MCMC moves and monitors, then runs the MCMC analysis
- the `model_FIG.Rev` script instantiates the feature-informed GeoSSE model components needed for likelihood-based inference
- the `geo_FIG.Rev` script reads and parameterizes the geographical features of regions that are used to inform rates in the FIG model
