import os
import sys
import data_processing
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

# Needed for testing the models
import run_merged
import run_grenc
import run_seqenc
import run_strenc

if __name__ == '__main__':

    # The user will first be asked a couple of question to determine what model they want to use
    # what sequences (in the form of a fasta file) they want to predict,
    # and lastly, if they chose Merged or GrEnc, they will have to provide a graph features file.
    # The model outputs will be saved in a txt file with each sequence identifier and the according probability.

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
              "For how to create new gspan feature files see section 'GraphProt' of the readme\n")
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
        print(f"The output will be provided in results/{fasta_file_input.split('.')[0]}_{model}_predictions.txt\n")
    else:
        fasta_file_input = "merged_test_file_30.fasta"
        graph_input = "graphprot_output/merged_test_file_30.gspan.gz.feature"
        struct_input = "merged_test_file_30_pysster.txt"

    # Reading in the fasta
    sequence_df = data_processing.read_fasta_file(fasta_file_input, labels=False)

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
