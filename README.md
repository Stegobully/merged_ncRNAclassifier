# merged_ncRNAclassifier: Classifying ncRNA sequences based on sequence, structure and graph encoding
Machine Learning tools for ncRNA classification 

# Introduction

Welcome to the Merged ncRNA classifier. ncRNAs are RNA sequences that do not encode for proteins.
They have been found to perform various tasks within the cell and have been linked with various medical conditions
including cancer. We use a mix of Convolutional Neural Networks (CNNs) and Artificial Neural Networks (ANNs) to classify
non-coding RNA into one of the 6 main types (lncRNA, miRNA, rRNA, snRNA, snoRNA and tRNA). The highest 
classification accuracy is achieved with a hybrid CNN/ANN model that concatenates a CNN branch using the primary
sequence as input with an ANN model using encoded graph features created by GraphProt as input.

# Installation

Create a new virtual environment with a tool of your choice (e.g. PyCharm). Clone the git repository into this folder and make sure you install the following packages:

|Package|Version|
|---|---|
|pandas   |1.5.3   |
|h5py   |3.8.0   |
|matplotlib   |3.6.3   |
|keras   |2.11   |
|tensorflow   |2.11   |
|forgi   |2.0.2   |
|numpy   |1.24.1   |
|scipy   |1.10.0   |
|scikit-learn   |1.2.1   |
|pysster   |1.2.1   |
|plotnine   |0.10.1   |
|mizani   |0.8.1   |

# How to use

We provide five different executable python files, predict_ncRNAs.py and test_grenc.py, test_merged.py, test_seqenc.py and test_strenc.py. 

## predict_ncRNAs.py
This program needs no additional parameters and should be run from console. The program will ask you, which of the four available models you want to test. Choices are:
### Merged
This classifier has the highest accuracy but requires the graph feature file created by GraphProt for the classification. For the highest accuracy you need 
