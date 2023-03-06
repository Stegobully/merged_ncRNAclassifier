# Merged_ncRNAclassifier (MncR): Improved late integrated ML-Classifier to predict ncRNA classes using sequence and graph encoded structures

# Introduction

The merged_ncRNAclassifier allows you to predict six different classes of non-coding (nc)RNAs. In comparison to protein coding messenger (m)RNAs, ncRNAs perform a variety of different “housekeeping” and “regulatory" tasks within the cell and can be further classified into ncRNA classes varying based on sequential, functional and/or structural features. Our approach compared the prediction capacity of ML classifiers using different features as input and model architectures based on Convolutional Neural Networks (CNNs) and/or Artificial Neural Networks (ANNs) ([Publication](#publication)). Our provided classifiers explained in detail below can predict ncRNAs into one of the six ncRNA classes: lncRNA, miRNA, rRNA, snRNA, snoRNA and tRNA. The different ML classifiers differentiate in their inputs using either only sequential information from the Fasta Sequence (SeqEnc), structural information (StrEnc) using the tool [Pysster](#pysster), graph encoded structural information (GrEnc) using the tool [GraphProt](#graphprot) or the best performing merged model (Merged) combining the primary sequence and structural graph encoded features as input.

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
The Git repository contains six different python scripts. To perform a step-by-step manual in the console you can use the `predict_ncRNAs.py`, which allows you to choose between the four different ML classifiers. If you prefer a specific ML classifier you can also directly use one of the four provided scripts for the different models: run_strenc.py, run_grenc.py, run_seqenc.py or run_merged.py. If you want to benchmark your tool against our ML classifier you can use the python script `benchmark_classifiers.py`, which additionally to the functions of `predict_ncRNAs.py` also adds the additional output of the prediction scores (Precision, Recall, F1, Matthews Correlation Coefficient). This script requires the labels of each sequence to be provided in the fasta file headers, so that the scores can be calculated.

# File `predict_ncRNAs.py`: ncRNA classifier with step-by-step manual in the console
The executable python script `predict_ncRNAs.py` provides a step-by-step manual in the console and always needs a standard fasta file of RNA sequences as input for classification into tRNA, rRNA, snoRNA, snRNA, miRNA and lncRNA. The user will be asked to enter the path to the output file. Alternatively, the output file for all individual ncRNA sequence results will be `[input-path]_[modelname]_predictions.txt`, where \[input-path\] is replaced by the path to the original file without the file ending and \[modelname\] is replaced by the chosen ML classifier. For each header from your input file the output file contains: ‘name\tpredicted class\tconfidence score\n’, where name corresponds to the name of the sequence in the input file. The predicted ncRNA type can be "lncRNA", "miRNA", "rRNA", "snRNA", "snoRNA" or "tRNA" and the confidence score is between 0 (not confident at all) and 1 (confident) based on the softmax function in the output layer of the ML classifier.
If you want to use the StrEnc ML classifier in addition to the fasta sequence input, you have to prepare a Pysster format file using [Pysster](#pysster) and for the usage of GrEnc or Merged a [GraphProt](#graphprot) graph feature file is required. The individual steps are: 

## 1. ML classifier choice: You can choose between the following four models, each needing different input for the classification.
### Merged
In our benchmark, this classifier had the best overall accuracy and F1-score between all four ML classifiers. Besides the fasta file of ncRNA sequences the classifier requires the graph feature file created by [GraphProt](#graphprot) as input for the classification. For how to create graph feature files see section [GraphProt](#graphprot). For more detailed information about the merged model see [run_merged.py](#run_mergedpy). 
### StrEnc
In our benchmark, this classifier performed best for the classes tRNA and miRNA with respect to the F1-score and overall second best. The ML classifier needs a structure encoding file based on the fasta file of ncRNA sequences. For how to create the structure encoding file see section [Pysster](#pysster). For more detailed information about the model see [run_strenc.py](#run_strencpy). 
### SeqEnc
This classifier performed best for the class rRNA in our benchmark. As input the ML classifier directly uses the fasta format file (Header starting with > and in the next line then the sequence). For more detailed information about the model see [run_seqenc.py](#run_seqencpy). 
### GrEnc
As input this ML classifier uses a fasta file transformed into the graph encoded file by GraphProt. Additionally, you will may also provide the original fasta file, because the GraphProt output does not contain sequence identifiers/names. If you do not provide the fasta, the output for each sample will simply be labeled `sequence_n`, where n is the row in which the feature vector is within the feature file. The sequences are not used for the classification, only the graph features. For how to create graph feature files see section [GraphProt](#graphprot). For more detailed information about the model see [run_grenc.py](#run_grencpy).

---
## 2. Output file naming

You will be asked, if you would like to provide the path to the output file. This may either be an absolute path or simply a file name that will be saved in the current working directory. Alternatively, you may also simply press 'enter', which will create `[inputname]_[modelname]_predictions.txt`, where `inputname` is simply the name of the input file (fasta input for Merged and SeqEnc; Structure input for StrEnc; Graph Feature input for GrEnc) up until the first "." and model name is the name of the model you chose in the previous step.

E.g.: If you chose the SeqEnc model and later 'sequencefile.fasta' as input, the output will be saved in 'sequencefile_seqenc_predictions.txt'.

### 3. Input Sequence Files

You will be asked to enter the full path and name of the input file or input files for the Merged classifiers. The only requirement to the sequence input is, that it needs to be readable by BioPython's `SeqIO.parse` method. If you type "default" the ML classifier will perform the output on a small subset from test set of our publication including 180 sequences (see `merged_test_file_XX `).
If the ML classifier GrEnc or Merged are selected in step 1 you will need to provide a `.feature`-file created by GraphProt. If you chose the StrEnc classifier, you will need to provide the path to the output created by pysster. 
Alternatively, you can also validate the results from out publication. For how to choose the test files used in our results, see [Datasets](#datasetsavailableforretrainingthemodelsandtesting)


# Standalone versions of the ML classifiers able to select in the python script ‘predict_ncRNAs.py’:


## run_merged.py
This model uses a late integration of the GrEnc and SeqEnc model by concatenating them within the last hidden layer before the softmax output layer. The SeqEnc model consists of a convolutional neural network (CNN) for the sequence input and the GrEnc a fully connected neural network (NN) for the GraphProt input. The graph input file and fasta file need to have the same order of sequences in the input files. Both files need to be provided as a run parameter when executing the python script. If you are interested in the exact architecture of the network load `models/merged_fold7.hdf5` in a python script using keras' `load_model()` method and run `model.summary()`. Alternatively, you can find the architecture described in the [Publication](#publication). 

Example call: `python test_merged.py path/to/fasta.fasta path/to/graph_enc.gspan.gz.feature`

## run_strenc.py
This model uses the structure encoding created by [pysster](#pysster) combined with the primary sequence as input. Pysster uses RNAfold by [Vienna RNA](http://rna.tbi.univie.ac.at/) to predict the secondary structure and then annotates each nucleotide with the substructure it belongs to (F: 5'-end; I: Internal Loop; M: Multi Loop; S: Stem; H: Hairpin Loop; T: 3'end). This encoding is derived from the dot-bracket notation by Pysster. 

Our model reads in the original nucleotide sequence as well as the structure sequence from the file created by Pysster. These two sequences are then combined into one sequence using an arbitrary encoding, for instance if at position i of the original nucleotide sequence an "A" nucleotide is present, and position i of the structure sequence created by Pysster is part of an inner loop (encoded as "I" by Pysster), then the i-th position of our sequence is encoded as an "E". The codes for all combinations of nucleotides and structure codes can be found in Table A1 of the [Publication](#publication). This combined sequence is then used as input for a convolutional neural network (CNN). If you are interested in the exact architecture of the network load `models/strenc_fold7.hdf5` in a python script using keras' `load_model()` method and run `model.summary()`. Alternatively, you can find the architecture described in the [Publication](#publication). 

Example call: `python test_strenc.py path/to/pysster_output.txt`


## run_grenc.py
This model uses the graph encoding created by [GraphProt](#graphprot) as input for an ANN model. The graph encoding comes in the form of 32,768 features in a feature vector. This encoding is first derived within GraphProt using secondary structure prediction with the tool ‘RNAshapes’ ([PubMed](https://pubmed.ncbi.nlm.nih.gov/16357029/)). The secondary structure is then encoded into a vector in GraphProt using a graph kernel approach that derives long distance graphical features of the secondary structures. The ML classifier is built with a fully connected Neural Network (NN). If you are interested in the exact architecture of the network load `models/grenc_fold2.hdf5` in a python script using keras' `load_model()` method and run `model.summary()`. Alternatively, you can find the architecture described in the [Publication](#publication). For the model to run, you will need to provide the feature file as a run parameter. Since the feature files do not contain sequence identifiers, the output identifiers will simply be labeled "sequence_n" where n is the row in which the feature vector lies within the file. If you wish to keep the original identifiers, you may additionally provide the fasta file that was used for the creation of the GraphProt feature file. The sequences will not be used for the classification, but the identifiers will be written to the output, meaning the sequences have to be in the same order as the graph feature files.

Example call without fasta: `python test_grenc.py path/to/graph_enc.gspan.gz.feature`
Example call with fasta: `python test_grenc.py path/to/graph_enc.gspan.gz.feature path/to/fasta.fasta`


## run_seqenc.py
This model uses just the primary sequence of a fasta file as input. The sequence is encoded into numerical values and padded to a length of 12,000 nt. For the classification, a convolutional neural network (CNN) is used. If you are interested in the exact architecture of the network load `models/seqenc_fold8.hdf5` in a python script using keras' `load_model()` method and run `model.summary()`. Alternatively, you can find the architecture described in the [Publication](#publication). 

Example call: `python test_seqenc.py path/to/fasta.fasta`

---

## File “benchmark_classifiers.py”: Benchmark our ML classifiers on given test datasets including ncRNA labels
This python script allows you to check the performance of our classifiers on a labelled set of ncRNA sequences. As input the script needs the corresponding input files for the different ML classifiers provided in the Git repository. Also, an additional text file is needed with the same order of the input sequences for the fasta file and/or additional structural files. The file has to have one line for each sequence in the same order as the input. Each line needs to contain one of lncRNA\n, miRNA\n, rRNA\n, snRNA\n, snoRNA\n or tRNA\n (Important: case sensitive). Additionally to the aforementioned `_[modelname]_predictions.txt`, this method also creates a text file `classification_scores.txt` in the working directory, which displays a classification report from scikit-learn (including class-wise and overall recall, precision and F1-score) as well as the Matthews Correlation Coefficient (MCC) for the prediction. A confusion matrix as a figure will be created using plotnine in `confusion_matrix.png`. 

# Input files: How to create structural input files for ncRNA classification
## Fasta File
The Fasta file can have a `.fasta` or `.fa` ending. It does not matter if after the header line starting with `>` the sequence is written in a single line or with line breaks over multiple lines, as both of these methods are readable by BioPythons SeqIO module. 

## GraphProt
To install GraphProt, refer to https://github.com/dmaticzka/GraphProt. You will need to be able to execute `fasta2shrep_gspan.pl` and `EDeN`. To create the Graph Feature files, execute:

1. `path/to/fasta2shrep_gspan.pl -abstr -stdout -M 3 -wins '150,' -shift '25' -fasta path/to/fasta -t 3 | gzip > path/to/gspan`

2. `path/to/EDeN -a FEATURE -i path/to/gspan`, 

where path/to/fasta is the path to the fasta file you want to predict and path/to/gspan is a new file that ends in `.gspan.gz`. This will create a new file in the same directory as your gspan file that ends in `.gspan.gz.feature`. This is the graph features file you will need to provide to the Merged and GrEnc models.
If you successfully created the `.feature` file, then you are ready to predict the sequences from the fasta file.


## Pysster
If you want to use the StrEnc model for classification, you will need to provide a structure file created by Pysster.
To install Pysster, refer to https://github.com/budach/pysster. You need to be able to import and run `predict_structures()`.
In Python, run:

1. `from pysster.utils import predict_structures`
2. `predict_structures(path/to/fasta, path/to/pysster.txt, annotate=True)`, 

where path/to/fasta is the corresponding file to the sequences you want to predict and path/to/pysster is the path to the file you will have to provide when choosing the StrEnc model. Before using the model, confirm that the pysster output is of the following format:

`>sequence name`
`nucleotide sequence`
`structure sequence`

Example: 

`>URS0000CB8740_9606`
`ATCTGCTCGCCGGAGCTCACTCT`
`FFFFSSSSHHHHSSSSTTTTTTT`

Every third line needs to contain the structure sequence made up of one of the six structure codes (5'-End: "F"; Stem: "S"; Inner Loop: "I"; Multi Loop: "M"; Hairpin Loop: "H"; 3'-End: "T"). The python script then creates the input for the StrEnc model by combining the second and third line into an arbitrary encoding, e.g. an A nucleotide as part of an inner loop ("I") is encoded as the letter "E". The full table for each combination of nucleotide and structure can be found in the [Publication](#publication). For the above example, the created input sequence would be

`QLULFIYIJAAJFWFIBSZSBSB`

## Datasets available for retraining the models and testing
Uploaded you find the training and validation fasta and structural files for our best performing models GrEnc, StrEnc, SeqEnc and Merged. The sequences are all downloaded from the RNAcentral (ref) and contain the …. 

Training and validation data sets for model Merged and StrEnc -> `fasta_files/train_ StrEnc.fasta`; `fasta_files/val_StrEnc.fasta`

Training and validation data sets for model GrEnc -> `fasta_files/train_GrEnc.fasta`; `fasta_files/val_GrEnc.fasta`

Training and validation data sets for model SeqEnc -> `fasta_files/train_SeqEnc.fasta`; `fasta_files/val_SeqEnc.fasta`
Training and validation data sets for model Merged -> `fasta_files/train_SeqEnc.fasta`; `fasta_files/val_SeqEnc.fasta`

For our benchmarks we used the test set -> `fasta_files/rnacentraltestset.fasta`



# Publication

XXXXXXXXXXXXXX Insert Link to Publication here



