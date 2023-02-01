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

Create a new virtual environment with Python >=3.9 with a tool of your choice (e.g. PyCharm). Clone the git repository into this folder and make sure you install the following packages:

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

We provide five different executable python files, ´predict_ncRNAs.py´ and ´test_grenc.py´, ´test_merged.py´, ´test_seqenc.py´ and ´test_strenc.py´. 

## predict_ncRNAs.py
This program needs no additional parameters and should be run from console. The program will ask you, which of the four available models you want to test. Choices are:
### Merged
This classifier has the highest accuracy but requires the graph feature file created by GraphProt for the classification. For how to create graph feature files see section "GraphProt".
### StrEnc
This classifier has the second highest accuracy but requires the structure annotation file created by pysster. While creating the structure annotation is faster than retrieving the graph feature vectors, if you already have both secondary structure encodings, this model will take longer due to the additional transformation of the input by combining sequence and structure. 
### SeqEnc
This classifier is faster than Merged and StrEnc and works using just ´.fasta´ files as input. The accuracy is slightly lower due to the loss of information about secondary structure. 
### GrEnc
This classifier is the fastest at classification, but shows lower scores. If you already have the graph feature files you are better off using Merged. If you wish to test this model anyway, you may do so by choosing this option.

---

Once you have chosen a model, you will be asked whether you would like to classify new sequences or test the model with already classified sequences. 
### new
You will only need to provide a ´.fasta´ file and the secondary structure (if necessary). The output will be a file in the results folder with the same name as the fasta file, but with ´_predictions.txt´ as a suffix. One line of this file may look like this: ´URS00019662A9_9685	lncRNA	0.9983231425285339´, where the first item is the sequence ID provided by the ´.fasta´ file (everything before the first space in the header), the second item is the predicted ncRNA type (one in lncRNA, miRNA, rRNA, snRNA, snoRNA, tRNA) and the third item is the probability output by the softmax function in the output layer of the model. The order of the sequence IDs is the same as in the corresponding ´.fasta´ file.

### test
If you choose this option, the fasta file has to have headers of following type: ´>sequenceid rnatype´, where rnatype is one of "lncRNA", "miRNA", "rRNA", "snRNA", "snoRNA" or "tRNA" (case sensitive). Additionally to the afformentioned ´_predictions.txt´, this method also outputs a ´classification_scores.txt´ which displays scikit-learn's classification report (including class-wise and over all recall, precision and F1-score) and also the Matthews Correlation Coefficient for the prediction. Additionally, ´confusion_matrix.png´ shows the confusion matrix (normalized over each row) for the prediction. Note, that these files are not specific to the fasta input, meaning rerunning the model will overwrite them. 
