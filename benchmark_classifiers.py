import os
import output_analysis
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

import predict_ncRNAs

if __name__ == '__main__':

    # The user will first be asked a couple of question to determine what model they want to use,
    # whether they want to test the model with annotated sequences or predict new ones,
    # what sequences (in the form of a fasta file) they want to predict,
    # and lastly, if they chose Merged or GrEnc, they will have to provide a graph features file.
    # If the user wishes to test sequences, the fasta files will need to have the label for each sequence
    # after each identifier (e.g.: >URS00009C41D9_10090 lncRNA) in the fasta file.
    # The model outputs will be saved in a txt file with each sequence identifier and the according probability.
    # In the case of testing, classification report and mcc are saved in results/classification_scores.txt
    # and a normalized confusion matrix is saved in results/confusion_matrix.png

    print("\nWelcome to the benchmark tool for ncRNA classification\n"
          "This program has the same functionality as predict_ncRNAs.py, but you will need to provide the classifier\n"
          "with the labels of each sequence in the following way: \n"
          "Provide the path to a .txt file where each line consists of the RNA type for the sequence.\n"
          "The n-th line in the labels file has to correspond to the n-th sequence in the input file.\n\n"
          "Additionally to the results, the program will provide two files in the current working directory,\n"
          "classification_scores.txt and confusion_matrix.png, in which scikit-learn's classification report and\n"
          "a normalized confusion matrix are found.\n\n"
          "Which model would you like to benchmark?\n"
          "Merged provides the best results, but requires Graph Features predicted by GraphProt\n"
          "StrEnc provides the second best result, but requires Structural Encoding created by Pysster\n"
          "SeqEnc predicts using only the sequence, but does not require additional files\n"
          "GrEnc predicts using only Graph Encoding, but requires Graph Features predicted by GraphProt\n")

    # Run the predictions
    results = predict_ncRNAs.predict_ncRNAs()

    # Prompt the user to enter the file with the labels for each sequence
    print("Now you will need to provide the labels for each sequence in a txt file.\n"
          "This file needs to only contain the RNA types, which for each line needs to be one of\n"
          "lncRNA\\n, miRNA\\n, rRNA\\n, snRNA\\n, snoRNA\\n or tRNA\\n (case sensitive)\n"
          "If you are using the default test files, you may also enter 'default'.")
    labels_file = "None"
    while not (os.path.isfile(labels_file) or labels_file.upper() == "DEFAULT"):
        labels_file = input("Please enter a valid path to the labels file or enter 'default'\n")

    # Set the standard labels file
    if labels_file.upper() == "DEFAULT":
        labels_file = "testing_datasets/small_testset_30_labels.txt"

    # Read the labels into a list
    labels_file = open(labels_file, "r")
    labels = labels_file.readlines()
    labels = list(map(lambda x: x.strip("\n"), labels))
    labels_file.close()

    # Print classification report, MCC and confusion matrix to files
    output_analysis.return_output_analysis(labels, results)
    output_analysis.plot_confusion_matrix(labels, results)
