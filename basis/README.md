# Preparing MS basis set

Selection of the basis can be divided into four stages.

The first stage is **to limit the range of substances that can be encountered in the particular experiment**. This step is absolutely problem-specific. For the thermal decomposition of PE in an inert atmosphere, we can limit the range of substances by carbohydrogens.

The next step is **to search for mass-spectra of these compounds**. We get the mass spectra of hydrocarbons from the NIST Chemistry WebBook database using the [NistChemPy](https://github.com/EPiCs-group/NistChemPy) package ([retrieve_carbohydrogens.py](retrieve_carbohydrogens.py), [CnHm_raw.csv](csvs/CnHm_raw.csv)). After that we repeated the first step, removing highly branched aliphatic compounds and aromatic compounds with several aromatic nuclei ([filter_carbohydrogens.py](filter_carbohydrogens.py), [CnHm_filtered.csv](csvs/CnHm_filtered.csv)).

The third step is **to filter the basis set by rejecting substances with similar spectra**. To do this, we use hierarchical clustering followed by filtering of similar compounds, i.e. various isomers of dimethylbenzene ([ms_clustering.py](ms_clustering.py), see the dendrograms in the [dendrograms](dendrograms) directory). This procedure can be automated, but in this project we limited ourselves to manual data processing.

The last step is to discard compounds that do not contribute to the selected spectra. This step is optional but provides significant acceleration of the deconvolution process. The final basis set is stored in the [basis.csv](csvs/basis.csv) file.
