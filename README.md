# PNAS2018
## Data and analysis code deposition for Vavrek et al, PNAS 2018

This public git repo contains the data and relevant analysis code for Vavrek et al, PNAS 2018. The raw 32768-channel HPGe spectra are stored in `.TKA` format, which is just a text format giving the counts per channel. This data is kept in the `hpge/` directory (~50 MB total). The top-level directory also contains `conciseData.txt`, which summarizes the beam current and live time measurements.

The analysis code is given in the ROOT script `analysis.C`

Example usage: generate all the NRF spectrum comparison plots:

$ root

root> .L analysis.C

root> plotAll()

NOTE: this may take a few minutes to run, especially for the first time as ROOT will need to generate autodicts for vectors of TH1D's, etc.

