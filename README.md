# merged_ncRNAclassifier: Improved late merge ML-Classifier to predict ncRNA classes using sequence and graph encoded structures

# Introduction

The merged_ncRNAclassifier allows you to predict six different classes of non-coding (nc)RNAs. In comparison to protein coding messenger (m)RNAs, ncRNAs perform a variety of different “housekeeping” and “regulatory" tasks within the cell and can be further classified into ncRNA classes differing based on sequential, functional and/or structural features. Our approach compared the prediction capacity of ML classifiers using different features as input and model architectures based on Convolutional Neural Networks (CNNs) and/or Artificial Neural Networks (ANNs) ([Publication](#publication). Our provided classifiers explained below in detail can predict ncRNAs into one of the six ncRNA classes: lncRNA, miRNA, rRNA, snRNA, snoRNA and tRNA. The different ML classifier differentiate in their inputs using wether only sequential information from the Fasta Sequence (SeqEnc), structural information (StrEnc) using the tool Pysster (Ref), graph encoded structural information (GrEnc) using the tool Graphprot (Ref) or the best performing merged model (Merged) combining the primary sequence and structural graph encoded features as input.

# Installation

We recommend creating a virtual environment with Python >=3.10 for the usage of our merged_ncRNAclassifier. Clone the git repository into your folder, assign it to the virtual environment and make sure you have installed the following packages for your virtual environment:

|Package|Version|
|---|---|
|pandas   |1.5.3   |
|h5py   |3.8.0   |
|matplotlib   |3.6.3   |
|keras   |2.11   |
|tensorflow   |2.11   |
|numpy   |1.24.1   |
|scipy   |1.10.0   |
|scikit-learn   |1.2.1   |
|plotnine   |0.10.1   |

# How to use
The Git repository contains six different python scripts. To perform a step-by-step manual in the console you can use the `predict_ncRNAs.py`, which allows you to choose between the four different ML classifier. If you prefer a specific ML classifier you can also directly use one of the four provided scripts for the different models: run_strenc.py, run_grenc.py, run_seqenc.py or run_merged.py. If you want to benchmark your tool against our ML classifier you can use the python script “XX”, which adds ….

# File `predict_ncRNAs.py`: ncRNA classifier with step-by-step manual in the console
The executable python script `predict_ncRNAs.py` provides a step-by-step manual in the console and needs as input always a standard Fasta file of RNA sequences for classification in tRNA, rRNA, snoRNA, snRNA, miRNA and lncRNA. The output will be written directly in the working directly or can be defined using the standard output in the console. The output file will be named after the prefix option asked in the second step and create a folder with the prefix.  The output file for all individual ncRNA sequence results will be `prefix_[modelname]_predictions.txt`, where [modelname] is replaced by the chosen ML classifier. The output file contains for each header from your input file: ‘name\tpredicted class\tconfidence score\n’. The predicted ncRNA type can be "lncRNA", "miRNA", "rRNA", "snRNA", "snoRNA" or "tRNA" and the confidence score is between 0 (not confident at all) and 1 (confident) based on the softmax function in the output layer of the ML classifier. Beside the output .txt file 
If you want to use the StrEnc ML classifier you have in addition to prepare a Pysster format file using Pysster (ref) and for the usage of GrEnc or Merged as second input file a Graphprot (ref) output is required.  The individual steps are: 

## 1. ML classifier choice: You can choose between the following four models, each needing different input for the classification.
### Merged
This classifier had in our benchmarks the best overall accuracy and F1-score between all four ML classifiers. As input beside the Fasta file of ncRNA sequences the classifier requires the graph feature file created by GraphProt (ref) for the classification. For how to create graph feature files see section [GraphProt](#graphprot). For more detailed information about the merged model see [run_merged.py](#run_mergedpy). 
### StrEnc
This classifier performed in our benchmarks best for the classes XXX and YYy. The ML classifier needs a structure encoding file based on the Fasta file of ncRNA sequences. For how to create the structure encoding file see section [Pysster](#pysster). For more detailed information about the model see [run_strenc.py](#run_strencpy). 
### SeqEnc
This classifier performed in our benchmarks best for the classes XXX and YYy. As input the ML classifier directly uses the Fasta format file (Header starting with > and in the next line then the sequence).For more detailed information about the model see [run_seqenc.py](#run_seqencpy). 
### GrEnc
This classifier performed in our benchmarks best for the classes XXX and YYy. As input the ML classifier uses the Fasta format file transformed in the graph encoded file by Graphprot. For how to create graph feature files see section [GraphProt](#graphprot). For more detailed information about the model see [run_grenc.py](#run_grencpy).

--- 


### test
If you choose this option, the fasta file has to have headers of following type: `>sequenceid rnatype`, where rnatype is one of "lncRNA", "miRNA", "rRNA", "snRNA", "snoRNA" or "tRNA" (case sensitive). Additionally to the afformentioned `_modelname_predictions.txt`, this method also outputs a `classification_scores.txt` which displays scikit-learn's classification report (including class-wise and over all recall, precision and F1-score) and also the Matthews Correlation Coefficient for the prediction. Additionally, `confusion_matrix.png` shows the confusion matrix (normalized over each row) for the prediction. Note, that these files are not specific to the fasta input, meaning rerunning the model will overwrite them. 

---
## 2. Output file naming
Giving the name of the prefix of your result text file. If nothing here is specified, the output file will be named result as default.
### 3. Input Sequence Files


You will be asked to enter the full path and name of the  input file or input files for the merged ML classifier without the file format ending like “.fa”. The only requirement is, that it needs to be readable by BioPython's `SeqIO.parse`, additional to the ncRNA types in the header of each sequence for the test option. If you type nothing or  "default" the ML classifier will perform the output on the test set from our publication including XXX sequences (see `merged_test_file_XX `). 
If the ML classifier GrENC or merged are selected in step 1 the file names with the ending `.gspan.gz.feature` are used or for StrENCthe `.pysster.txt.



The output consists of a .txt file in the results folder. It has the same name as the fasta file (without the .fasta suffix) followed by `\_[model]\_predictions.txt`, where model is the chosen model. If you chose the "test" option, two more files will be output to the results folder, `classification_scores.txt` containing class-wise recall, precision and F1-score as well as over all recall, precision, F1-score and MCC, and `confusion_matrix.png`, which contains the normalized confusion matrix of the RNA types, created using Plotnine. 


# Standalone versions of the ML classifiers able to select in the python script ‘predict_ncRNAs.py’:


## run_merged.py
This model used a late integration of  the GrEnc and SeqEnc model by concatenating them with an additional last hidden layer before the softmax and output layer.. The SeqEnc model consists of a convolutional neural network (CNN) for the sequence input and the GrEnc a fully connected neural network (NN) for the graphprot input. The graph input file and fasta file need to have the same order of sequence ids in the input files. Both files need to be provided as a run parameter when executing the python script. If you are interested in the exact architecture of the network load `models/merged_fold7.hdf5` using keras' `load_model()` method and run `model.summary()`. Alternatively, you can find the architecture described in the [Publication](#publication). 

Example call: `python test_merged.py path/to/fasta.fasta path/to/graph_enc.gspan.gz.feature`

## run_strenc.py
This model uses the structure encoding created by [pysster](#pysster) combined with the primary sequence as input. Pysster uses RNAfold by Vienna RNA (ref) to predict the secondary structure and then annotates each nucleotide with the substructure it belongs to (F: 5'-end; I: Internal Loop; M: Multi Loop; S: Stem; H: Hairpin Loop; T: 3'end) (ref). The nucleotide and structure sequence are then combined through arbitrary encoding of each combination of nucleotide and. This combined sequence is then used as input for a convolutional neural network (CNN). If you are interested in the exact architecture of the network load `models/strenc_fold7.hdf5` using keras' `load_model()` method and run `model.summary()`. Alternatively, you can find the architecture described in the [Publication](#publication). 

Example call: `python test_strenc.py path/to/fasta.fasta path/to/structure_pysster.txt`

## run_grenc.py
This model uses the graph encoding created by [GraphProt](#graphprot) as input for an ANN model. The graph encoding comes in the form of 32,768 features in a sparse vector. This encoding is first derived using secondary structure prediction with the tool ‘RNAstruct’ (Ref). The secondary structure is then encoded into a vector using a graph kernel approach that derives long distance graphical features of the secondary structures. The ML classifier is built with a fully connected Neural Network (NN). If you are interested in the exact architecture of the network load `models/grenc_fold2.hdf5` using keras' `load_model()` method and run `model.summary()`. Alternatively, you can find the architecture described in the [Publication](#publication). 

Example call: `python test_grenc.py path/to/fasta.fasta path/to/graph_enc.gspan.gz.feature`

## run_seqenc.py
This model uses just the primary sequence of a fasta file as input. is the sequence is encoded into numerical values and padded to a length of 12,000 nt. For the classification, a convolutional neural network (CNN) is used. If you are interested in the exact architecture of the network load `models/seqenc_fold8.hdf5` using keras' `load_model()` method and run `model.summary()`. Alternatively, you can find the architecture described in the [Publication](#publication). 

Example call: `python test_seqenc.py path/to/fasta.fasta`

---

## File “XXX.py”: Benchmark our ML classifiers on given test datasets including ncRNA labels
This python script allows you to check the performance of a labelled set of ncRNA sequences. As input the script needs the corresponding input files for the different ML classifier provided in the Git repository. Also an additional text file ‘.txt’ is needed with the same order of the input sequences for the fasta file and/or additional structural files. The ‘.txt’ file has to be the same name like the prefix of the other input files and is containing the following structure: ‘name \t label \n’. The name is the name form the ncRNA sequence in the other input files and the label is lncRNA, miRNA, rRNA, snRNA, snoRNA or tRNA (Important: case sensitive). Additionally, to the aforementioned `_[modelname]_predictions.txt`, this method also creates in the output folder a text file `classification_scores.txt`, which displays a classification report from scikit-learn (including class-wise and over all recall, precision and F1-score) as well as the Matthews Correlation Coefficient (MCC) for the prediction. A confusion matrix as image will be created `confusion_matrix.png`. 


# Input files: How to create structural input files for ncRNA classification
## Fasta File
The Fasta file can have as ending whether `.fasta` or `.fa`. It does not matter if after the header line starting with ‘>’ the sequence is written in a single line or with line breaks over multiple lines(readable by Biopythons SeqIO module). 


## GraphProt
To install GraphProt, refer to https://github.com/dmaticzka/GraphProt. You will need to be able to execute `fasta2shrep_gspan.pl` and `EDeN`. To create the Graph Feature files, execute:

1. `path/to/fasta2shrep_gspan.pl -abstr -stdout -M 3 -wins '150,' -shift '25' -fasta {fasta_file} -t 3 | gzip > {gspan_file}`

2. `path/to/EDeN -a FEATURE -i {gspan_file}`, 

where fasta_file is the path to the fasta file you want to predict and gspan_file is a new file that ends in `.gspan.gz`
Confirm the output file now ends in `.gspan.gz.feature`, then you are ready to predict the sequences from the fasta file.


## Pysster
To install Pysster, refer to https://github.com/budach/pysster. You need to be able to import and run `predict_structures()`.
In Python, run:

1. `from pysster.utils import predict_structures`
2. `predict_structures(fasta_file, output_pysster.txt, annotate=True)`, 

where fasta file is the corresponding file to your sequence. Make sure the Pysster output file ends in `_pysster.txt`.


## Datasets available for retraining the models and testing
Uploaded you find the training and validation fasta and structural files for our best performing models GrEnc, StrEnc, SeqEnc and Merged. The sequences are all downloaded from the RNAcentral (ref) and contain the …. 

Training and validation data sets for model Merged and StrEnc -> `fasta_files/train_ StrEnc.fasta`; `fasta_files/val_StrEnc.fasta`

Training and validation data sets for model GrEnc -> `fasta_files/train_GrEnc.fasta`; `fasta_files/val_GrEnc.fasta`

Training and validation data sets for model SeqEnc -> `fasta_files/train_SeqEnc.fasta`; `fasta_files/val_SeqEnc.fasta`
Training and validation data sets for model Merged -> `fasta_files/train_SeqEnc.fasta`; `fasta_files/val_SeqEnc.fasta`

For our benchmarks we used the test set -> `fasta_files/rnacentraltestset.fasta`



# Publication

XXXXXXXXXXXXXX Insert Link to Publication here



