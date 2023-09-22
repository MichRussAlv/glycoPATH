# glycoPATH

**Repository for correlating 3'-Tag RNAseq data with LC-MS/MS glycomics data**

The code consists of 2 main pipelines: calculating the correlations and creating an N-glycan biosynthesis network from the correlations, and predicting N-glycomic abundances from RNAseq gene counts. The first pipeline, correlations and network generation, uses the Python scripts glycogene_nglycan-correlation.py and network_generation.py, respectively.  The second pipeline, prediction, uses the predictions.py script.

# Installation

The scripts were written and tested in Python 3.9.7. The scripts use the following Python packages:

Correlations calculation
1. matplotlib.pyplot
2. numpy

Network generation
1. graphviz
2. pandas

Predictions
1. pandas
2. numpy
3. scipy.interpolate
4. pickle
5. glob
6. matplotlib.pyplot


# Correlations calculation and network generation

Both RNAseq (gene_averages.csv) and N-glycomic LC-MS/MS (n-glycan_averages.csv) data are required to calculate the correlations. Both genes and N-glycan compounds are categorized and filtered based on composition - high-mannose (HM), undecorated (U), fucosylated (F), sialylated (S), and sialofucosylated (SF). For heatmap and network visualization, the N-glycan structures in colored CFG notation are provided. All these files should be in the same working directory as the scripts glycogene_nglycan-correlation.py and network_generation.py.

# Predictions

Glycogene and N-glycan correlations are used to generate spline models saved in pickle files. The pickle files are called to calculate N-glycan values from RNAseq gene expressions.
