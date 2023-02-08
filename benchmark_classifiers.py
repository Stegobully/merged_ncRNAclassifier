import os
import sys
import data_processing
import output_analysis
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

# Needed for testing the models
import run_merged
import run_grenc
import run_seqenc
import run_strenc

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
          "This program has the same functionality as predict_ncRNAs.py, but you will need to provide the classifier"
          "with the labels of each sequence in the following way: For each sequence in the fasta file the header"
          "needs to be in the following format:\n"
          "'>SeqID rna_type', where rna_type is one of 'lncRNA', 'miRNA', 'rRNA', 'snRNA', 'snoRNA' or 'tRNA' for"
          "each sequence (case sensitive)\n"
          "Additionally to the results the program will provide two files in the results folder, "
          "classification score.txt and confusion_matrix.png, in which scikit-learn's classification report and "
          "a normalized confusion matrix are found.\n"
          "Which model would you like to benchmark?\n"
          "Merged provides the best results, but requires Graph Features predicted by GraphProt\n"
          "StrEnc provides the second best result, but requires Structural Encoding created by Pysster\n"
          "SeqEnc predicts using only the sequence, but does not require additional files\n"
          "GrEnc predicts using only Graph Encoding, but requires Graph Features predicted by GraphProt\n")

    # The following code asks the user what model they want to use to predict
    model = "None"
    while model.upper() not in {"MERGED", "STRENC", "SEQENC", "GRENC"}:
        model = input("Type one of 'Merged', 'StrEnc', 'SeqEnc' or 'GrEnc'\n")
        if model.upper() not in {"MERGED", "STRENC", "SEQENC", "GRENC"}:
            print("Please write one of the model names")
    print(f"You chose {model}\n\n")

    # Here the user has to enter the fasta file path
    fasta_file_input = ""
    while not ((os.path.isfile(fasta_file_input) & fasta_file_input.endswith((".fa", ".fasta")))
               or fasta_file_input.upper() == "DEFAULT"):
        fasta_file_input = input("Please enter a valid fasta path or 'default'\n")

    # Here the user has to enter the graph feature file path created by GraphProt if the model requires it
    if model.upper() in {"MERGED", "GRENC"} and fasta_file_input.upper() != "DEFAULT":
        graph_input = ""
        print("You will now need to enter the path to the corresponding graph feature file created by GraphProt\n"
              "This file needs to have the ending '.gspan.gz.feature'\n"
              "For how to create new gspan feature files see section XYZ of the readme\n")
        while not os.path.isfile(graph_input) & graph_input.endswith(".gspan.gz.feature"):
            graph_input = input("Please enter a valid path to a feature file\n")

    # Here the user has to provide the structure file path created by Pysster
    if model.upper() == "STRENC" and fasta_file_input.upper() != "DEFAULT":
        struct_input = ""
        print("You will now need to enter the path to the corresponding structure file created by Pysster\n"
              "This file needs to have the ending 'pysster.txt'\n"
              "For how to create new structure files see section 'Pysster' of the readme\n")
        while not (os.path.isfile(struct_input) & struct_input.endswith("pysster.txt")):
            struct_input = input("Please enter a valid path to a structure file\n")

    # Print output path. If the user wishes to test with provided sequences the necessary files are set to the default
    if fasta_file_input.upper() != "DEFAULT":
        print(f"The output will be provided in {fasta_file_input.split('.')[0]}_predictions.txt\n")
    else:
        fasta_file_input = "merged_test_file_30.fasta"
        graph_input = "graphprot_output/merged_test_file_30.gspan.gz.feature"
        struct_input = "merged_test_file_30_pysster.txt"

    # Read in fasta file
    sequence_df = data_processing.read_fasta_file(fasta_file_input, labels=True)

    # Run the Merged model
    if model.upper() == "MERGED":
        # Test if the length of the sequence df matches the number of graph feature vectors
        if not data_processing.test_graphfeat_seq_match(sequence_df, graph_input):
            sys.exit()
        results, pred_probabilities = run_merged.test_merged(sequence_df, graph_input, fasta_file_input)

    # Run the GrEnc model
    elif model.upper() == "GRENC":
        # Test if the length of the sequence df matches the number of graph feature vectors
        if not data_processing.test_graphfeat_seq_match(sequence_df, graph_input):
            sys.exit()
        results, pred_probabilities = run_grenc.test_grenc(sequence_df, graph_input, fasta_file_input)

    # Run the SeqEnc model
    elif model.upper() == "SEQENC":
        results, pred_probabilities = run_seqenc.test_seqenc(sequence_df)

    # Run the StrEnc model
    elif model.upper() == "STRENC":
        results, pred_probabilities = run_strenc.test_strenc(sequence_df, struct_input)

    # Save seq IDs for the result output
    sequence_ids = sequence_df.index

    # Write output to file
    if not os.path.isdir("results"):
        os.mkdir("results")
    output_file = f"results/{fasta_file_input.split('.')[0]}_{model}_predictions.txt"
    output = open(output_file, "w")
    for id, pred, pred_probability in zip(sequence_ids, results, pred_probabilities):
        output.write(f"{id}\t{pred}\t{pred_probability}\n")
    output.close()
    print(f"Results are saved in {output_file}")

    # Print classification report, MCC and confusion matrix to files
    labels = sequence_df.label
    output_analysis.return_output_analysis(labels, results)
    output_analysis.plot_confusion_matrix(labels, results)