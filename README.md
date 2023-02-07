# merged_ncRNAclassifier: Classifying ncRNA sequences based on sequence, structure and graph encoding
Machine Learning (ML) classifiers for ncRNA prediction

# Introduction

Welcome to the Merged ncRNA classifier. ncRNAs are non-coding RNA sequences that do not encode for proteins like messenger RNAs (mRNAs).
ncRNAs perform a variety of different “housekeeping” and “regulatory" tasks within the cell and can be further classified in RNA families based on
functional, structural and sequential features. Our ML classifiers are based on Convolutional Neural Networks (CNNs) and/or Artificial Neural Networks (ANNs)
to predict non-coding RNA sequences into one of six main ncRNA-types (lncRNA, miRNA, rRNA, snRNA, snoRNA and tRNA). The different ML classifier differentiate in their inputs using structural or sequential information, of which the merged model combining the primary sequence and structural graph features as input achieved the highest accuracy in our benchmark.

# Installation

Create a new virtual environment with Python >=3.10 with a tool of your choice (e.g. PyCharm). Clone the git repository into this folder and make sure you have installed the following packages:

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

## predict_ncRNAs.py
The executable python file `predict_ncRNAs.py` provides a step-by-step manual in the console to classify a fasta file of ncRNA sequences using one of our four implemented ML classifiers. Additionally, this script gives the user the option to validate the prediction if the true labels are known by outputting a classification report as well as the confusion matrix for the different ncRNA classes. The individual steps are: 

## 1. Choice of Model: You can choose between the following four models, each needing different input for the classification
### Merged
This classifier has the highest accuracy but requires the graph feature file created by GraphProt for the classification. For how to create graph feature files see section [GraphProt](#graphprot). For more information about the model see [run_merged.py](#run_mergedpy). 
### StrEnc
This classifier has the second highest accuracy but requires the structure annotation file created by pysster. While creating the structure annotation is faster than retrieving the graph feature vectors, if you already have both secondary structure encodings, this model will take longer due to the additional transformation of the input by combining sequence and structure. For more information about the model see [run_strenc.py](#run_strencpy). 
### SeqEnc
This classifier is faster than Merged and StrEnc and works using just `.fasta` files as input. The accuracy is slightly lower due to the loss of information about secondary structure. For more information about the model see [run_seqenc.py](#run_seqencpy). 
### GrEnc
This classifier is the fastest at classification, but shows lower scores. If you already have the graph feature files you are better off using Merged. If you wish to test this model anyway, you may do so by choosing this option. For how to create graph feature files see section [GraphProt](#graphprot). For more information about the model see [run_grenc.py](#run_grencpy).

--- 

## 2. Do you want to "test" the model with sequences whose RNA type is known or do you simply want to predict "new" sequences
### new
You will only need to provide a `.fasta` file and the secondary structure (if necessary). The output will be a file in the results folder with the same name as the fasta file, but with `_modelname_predictions.txt` as a suffix, where modelname is the one provided in the previous question. One line of this file may look like this: `URS00019662A9_9685	lncRNA	0.9983231425285339`, where the first item is the sequence ID provided by the `.fasta` file (everything before the first space in the header), the second item is the predicted ncRNA type (one of "lncRNA", "miRNA", "rRNA", "snRNA", "snoRNA", "tRNA") and the third item is the probability output by the softmax function in the output layer of the model. The order of the sequence IDs is the same as in the corresponding `.fasta` file.

### test
If you choose this option, the fasta file has to have headers of following type: `>sequenceid rnatype`, where rnatype is one of "lncRNA", "miRNA", "rRNA", "snRNA", "snoRNA" or "tRNA" (case sensitive). Additionally to the afformentioned `_modelname_predictions.txt`, this method also outputs a `classification_scores.txt` which displays scikit-learn's classification report (including class-wise and over all recall, precision and F1-score) and also the Matthews Correlation Coefficient for the prediction. Additionally, `confusion_matrix.png` shows the confusion matrix (normalized over each row) for the prediction. Note, that these files are not specific to the fasta input, meaning rerunning the model will overwrite them. 

---

## 3. Input Sequence File in Fasta Format


After choosing whether to test or predict new sequences, you will be asked to enter the path to a `.fasta` file. The entered file may also end in `.fa`. The only requirement is, that it needs to be readable by Biopython's `SeqIO.parse`, additional to the rna types in the header of each sequence for the test option. You may also choose the default option by typing "default". This will test the model on `merged_test_file_30.fasta`, which contains 30 sequences of each ncRNA type. 

## 4. (Optional) for model GrEnc/Merged and StrEnc: Input Graph Encoding or Structure Encoding files

Next, you will be asked to provide the link to the structure file, unless you chose "SeqEnc" for the model. The graphprot graph encoding file will need to end in `.gspan.gz.feature`, the pysster structure encoding will need to end in `pysster.txt`, otherwise the program will repeat the prompt to enter the file. If you chose the default option, you will not be asked to enter a structure/graph encoding file and the program will choose the corresponding file automatically. 

Lastly, the program will read in the `.fasta` file (and the structure file if needed) and predict with the chosen model. If the number of graph feature vectors in the provided `.feature` file and sequences in the `.fasta` file do not match, the program will throw an error before prediction. If prediction does not fail, the output is written to the corresponding file and if the option "test" was chosen, the plots are created and saved in the results folder.

## 5. Output will be created automatically

The output consists of a .txt file in the results folder. It has the same name as the fasta file (without the .fasta suffix) followed by "\_[model]\_predictions.txt", where model is the chosen model. If you chose the "test" option, two more files will be output to the results folder, "classification_scores.txt" containing class-wise recall, precision and F1-score as well as over all recall, precision, F1-score and MCC, and "confusion_matrix.png", which contains the normalized confusion matrix of the RNA types, created using Plotnine. 


Different models for selection in predict_ncRNAs.py or as standalone version:

## run_[model].py
The `test_[model].py` programs are standalone versions of each of the models. They require the fasta file and (if needed) the structure file as run parameters. The  `_[model]_predictions.txt` file. The order of the entered files is not flexible. 

### run_merged.py
This model combines the GrEnc and SeqEnc model by concatenating the last layer before the output. It has the highest accuracy of the four models and uses a mix of primary sequence and secondary structure graph encoding created by [GraphProt](#graphprot). The model consists of a convolutional neural network for the sequence input and a fully connected neural network for the graph input. The graph input file needs to contain the corresponding feature vectors in the same order as the sequence input. Both files need to be provided as a runtime parameter when executing the python script.

Example call: `python test_merged.py path/to/fasta.fasta path/to/graph_enc.gspan.gz.feature`

### run_strenc.py
Example call: `python test_strenc.py path/to/fasta.fasta path/to/structure_pysster.txt`

### run_grenc.py
Example call: `python test_grenc.py path/to/fasta.fasta path/to/graph_enc.gspan.gz.feature`

### run_seqenc.py
Example call: `python test_seqenc.py path/to/fasta.fasta`



---
Required Input files:

Primary Sequence 



IF you want to usethe model XXX oder YYY for predciting ncRNAs based on graph encoded you have beforehand to use the ncRNA sequence FASTA file as input for Graphprot to create the needed input file:

## GraphProt
To install GraphProt, refer to https://github.com/dmaticzka/GraphProt. You will need to be able to execute `fasta2shrep_gspan.pl` and `EDeN`. To create the Graph Feature files execute:

1. `path/to/fasta2shrep_gspan.pl -abstr -stdout -M 3 -wins '150,' -shift '25' -fasta {fasta_file} -t 3 | gzip > {gspan_file}`

2. `path/to/EDeN -a FEATURE -i {gspan_file}`, 

where fasta_file is the path to the fasta file you want to predict and gspan_file is a new file that ends in `.gspan.gz`
Confirm the output file now ends in `.gspan.gz.feature`, then you are ready to predict the sequences from the fasta file.


IF you want to usethe model XXX oder YYY for predciting ncRNAs based on graph encoded you have beforehand to use the ncRNA sequence FASTA file as input for Graphprot to create the needed input file:

## Pysster
To install Pysster, refer to https://github.com/budach/pysster. You need to be able to import and run `predict_structures()`.

`predict_structures(fasta_file, output_pysster.txt, annotate=True)`, 

where fasta file is the corresponding file to your sequence. Make sure the Pysster output file ends in `_pysster.txt`.


## Datasets available for retraining the models and testing

Training data set model Merged and StrEnc -> `fasta_files/train_fold_7.fasta`

Training data set model GrEnc -> `fasta_files/train_fold_2.fasta`

Training data set model SeqEnc -> `fasta_files/train_fold_8.fasta`

Validation files used have the same number as the training files

RNAcentral test set -> `fasta_files/rnacentraltestset.fasta`

To retrieve the exact model parameters for recreating the classifiers load the model hdf5 files (`model_files/[model]_foldx.hdf5`) using keras. You may then use `model.summary()` for the parameters. For further parameters see the publication. 
