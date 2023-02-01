import os
import sys
import data_processing
import output_analysis
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

# Needed for testing the models
import test_merged
import test_grenc
import test_seqenc
import test_strenc

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

    print("\nWelcome to ncRNA classification\n"
          "Which model would you like to predict sequences saved in fasta file format?\n"
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

    # The following prompt asks the user whether they want to predict new sequences or test labeled ones
    print("Do you wish to test labeled sequences or annotate new ones?"
          "\n'new' will provide a results.txt with the labels for each sequence.\n"
          "'test' will additionally provide a set of plots for you to analyze the test results.\n"
          "Note, that testing requires the header of each sequence to be '>SEQID rnatype'.\n"
          "rna type will need to be one of 'lncRNA', 'miRNA', 'rRNA', 'snRNA', 'snoRNA' or 'tRNA'\n")
    new_or_test = ""
    while new_or_test.upper() not in {"NEW", "TEST"}:
        new_or_test = input("Type 'new' for unlabeled sequences or 'test' for testing labeled sequences\n")
        if new_or_test.upper() not in {"NEW", "TEST"}:
            print("Please write either 'new' or 'test'")

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
              "For how to create new structure files see section XYZ of the readme\n")
        while not (os.path.isfile(struct_input) & fasta_file_input.endswith("pysster.txt")):
            struct_input = input("Please enter a valid path to a structure file\n")

    # Print output path. If the user wishes to test with provided sequences the necessary files are set to the default
    if fasta_file_input.upper() != "DEFAULT":
        print(f"The output will be provided in {fasta_file_input.split('.')[0]}_predictions.txt\n")
    else:
        fasta_file_input = "merged_test_file_30.fasta"
        graph_input = "graphprot_output/merged_test_file_30.gspan.gz.feature"
        struct_input = "merged_test_file_30_pysster.txt"

    # Option whether to read in the labels from the fasta file for testing
    if new_or_test.upper() == "TEST":
        sequence_df = data_processing.read_fasta_file(fasta_file_input, labels=True)
    else:
        sequence_df = data_processing.read_fasta_file(fasta_file_input, labels=False)

    # Run the Merged model
    if model.upper() == "MERGED":
        # Test if the length of the sequence df matches the number of graph feature vectors
        if not data_processing.test_graphfeat_seq_match(sequence_df, graph_input):
            sys.exit()
        results, pred_probabilities = test_merged.test_merged(sequence_df, graph_input, fasta_file_input)

    # Run the GrEnc model
    elif model.upper() == "GRENC":
        # Test if the length of the sequence df matches the number of graph feature vectors
        if not data_processing.test_graphfeat_seq_match(sequence_df, graph_input):
            sys.exit()
        results, pred_probabilities = test_grenc.test_grenc(sequence_df, graph_input, fasta_file_input)

    # Run the SeqEnc model
    elif model.upper() == "SEQENC":
        results, pred_probabilities = test_seqenc.test_seqenc(sequence_df)

    # Run the StrEnc model
    elif model.upper() == "STRENC":
        results, pred_probabilities = test_strenc.test_strenc(sequence_df, struct_input)

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
    if new_or_test.upper() == "TEST":
        labels = sequence_df.label
        output_analysis.return_output_analysis(labels, results)
        output_analysis.plot_confusion_matrix(labels, results)

"""if sys.argv == None:
    print("please enter fasta file")
elif not (sys.argv[1].endswith((".fasta", ".fa"))):
    print("Please provide a valid fasta file with .fasta or .fa ending")
elif not os.path.isfile(sys.argv[1]):
    print(f"{sys.argv[1]} does not exist.")
else:
    fasta_file = sys.argv[1]
    df = dp.read_fasta_file(fasta_file)
    dp.transform_seq_into_graphfeatures(df, fasta_file)
"""