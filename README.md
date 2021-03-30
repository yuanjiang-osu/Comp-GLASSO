# Comp-GLASSO

This paper analyzes both simulated and real data. The simulated data are generated with random number seed. The real data are OTU abundance data downloadable from the Tara Oceans Project data repository (http://taraoceans.sb-roscoff.fr/EukDiv/), including:

- Database_W5_OTU_occurences.tsv: OTU abundance data (zipped due to storage limit on GitHub)
- taxa_mapping.csv: taxonomic hierarchy for each OTU
- TARA_truth.csv: literature validated genus interactions

The code aggregates the OTU abudance data (Database_W5_OTU_occurences.tsv) into genus abundance data based on the taxonomic information (taxa_mapping.csv), estimates the genus interaction networks, and compares with the literature validated genus interactions (TARA_truth.csv).
