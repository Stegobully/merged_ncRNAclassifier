# merged_ncRNAclassifier: Classifying ncRNA sequences based on sequence, structure and graph encoding
Machine Learning tools for ncRNA classification 

# Introduction

Welcome to the Merged ncRNA classifier. ncRNAs are RNA sequences that do not encode for proteins.
They have been found to perform various tasks within the cell and have been linked with various medical conditions
including cancer. We use a mix of Convolutional Neural Networks (CNNs) and Artificial Neural Networks (ANNs) to classify
non-coding RNA into one of the 6 main types (lncRNA, miRNA, rRNA, snRNA, snoRNA and tRNA). The highest 
classification accuracy is achieved with a hybrid CNN/ANN model that concatenates a CNN branch using the primary
sequence as input with an ANN model using encoded graph features created by GraphProt as input.


