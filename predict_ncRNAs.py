import os
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

# Needed for testing the models
import run_merged
import run_grenc
import run_seqenc
import run_strenc

def predict_ncRNAs():

    # The user will first be asked a couple of question to determine what model they want to use
    # what sequences (in the form of a fasta file) they want to predict,
    # and lastly, if they chose Merged or GrEnc, they will have to provide a graph features file.
    # The model outputs will be saved in a txt file with each sequence identifier and the according probability.

    # The following code asks the user what model they want to use to predict
    model = "None"
    while model.upper() not in {"MERGED", "STRENC", "SEQENC", "GRENC"}:
        model = input("Type one of 'Merged', 'StrEnc', 'SeqEnc' or 'GrEnc'\n")
        if model.upper() not in {"MERGED", "STRENC", "SEQENC", "GRENC"}:
            print("\nPlease write one of the model names\n")
    print(f"You chose {model}\n")

    # Here the user has to enter the fasta file path if the model requires it
    fasta_file_input = ""
    if model.upper() in {"MERGED", "SEQENC"}:
        print("\nThe model you entered requires a fasta file as input for the prediction.\n\n")
        while not ((os.path.isfile(fasta_file_input) & fasta_file_input.endswith((".fa", ".fasta")))
                   or fasta_file_input.upper() == "DEFAULT"):
            fasta_file_input = input("Please enter a valid fasta path or 'default'\n")


    # Here the user has to enter the graph feature file path created by GraphProt if the model requires it
    graph_input = ""
    if model.upper() in {"MERGED", "GRENC"} and fasta_file_input.upper() != "DEFAULT":
        print("\nThe model you selected requires a graph features input file.\n"
              "You will now need to enter the path to the corresponding graph feature file created by GraphProt.\n"
              "This file needs to have the ending '.feature'.\n"
              "For how to create new gspan feature files see section 'GraphProt' of the readme.\n"
              "Alternatively, you can enter 'default' to choose the default option.\n")
        while not ((os.path.isfile(graph_input) & graph_input.endswith(".feature"))
                or graph_input.upper() == "DEFAULT"):
            graph_input = input("Please enter a valid path to a feature file or 'default'\n")

    # GraphProt feature files do not contain sequence identifiers
    # If the user wishes to provide sequence identifiers in the output, they will need to provide a path
    # to the corresponding fasta file used for the creation of the graph features file
    if model.upper() == "GRENC" and graph_input.upper() != "DEFAULT":
        fasta_file_input = "tmp"
        print("\nThe graph feature files for the GrEnc model do not contain sequence identifiers \n" 
              "If you want the output to contain sequence identifiers, you will need to provide the fasta file \n"
              "used for the creation of the graph feature file. Alternatively, the results will be in the same \n"
              "order as the feature vectors, but simply names 'sequence_n', \n"
               " where n is the row number of the feature vector.\n"
              "If you want to choose this option simply press enter.\n")
        while not ((os.path.isfile(fasta_file_input) & fasta_file_input.endswith((".fa", ".fasta"))
                   or fasta_file_input == "")):
            fasta_file_input = input("Please enter a valid fasta path or simply press enter\n")

    # Here the user has to provide the structure file path created by Pysster
    struct_input = ""
    if model.upper() == "STRENC":
        print("\nThe model you entered requires a structure file created by Pysster. \n"
              "You will now need to enter the path to the structure file of the sequences you wish to predict. \n"
              "For how to create new structure files see section 'Pysster' of the readme. \n"
              "You may also enter 'default' to choose the default option.\n")
        while not (os.path.isfile(struct_input) or struct_input.upper() == "DEFAULT"):
            struct_input = input("Please enter a valid path to a structure file or type 'default'\n")

    # If the user wishes to test with provided sequences the necessary files are set to the default
    if fasta_file_input.upper() == "DEFAULT"\
            or graph_input.upper() == "DEFAULT"\
            or struct_input.upper() == "DEFAULT":
        fasta_file_input = "testing_datasets/small_testset_30.fasta"
        graph_input = "graphprot_output/small_testset_30_graphprot.feature"
        struct_input = "pysster_output/small_testset_30_pysster.txt"

    # Ask user to provide a path to the output or write to the standard output name
    print("\nPlease input the path to the file, to which you want the predictions to be written.\n"
          "If you do not want to specify a file, simply press enter: The predictions will be written to\n"
          "[input_name]_[model]_predictions.txt in the current working directory,\n"
          "where 'input_name' is the name of the input file and model is the name of the chosen model.\n")
    output_file = input("Enter the path to the output file or press enter\n")

    # Run the Merged model
    if model.upper() == "MERGED":
        ids, results, pred_probabilities = run_merged.test_merged(fasta_file_input, graph_input)
        if output_file == "":
            output_file = f"{fasta_file_input.split('.')[0]}_merged_prediction.txt"

    # Run the GrEnc model
    elif model.upper() == "GRENC":
        ids, results, pred_probabilities = run_grenc.test_grenc(graph_input, fasta_file_input)
        if output_file == "":
            output_file = f"{graph_input.split('.')[0]}_grenc_prediction.txt"

    # Run the SeqEnc model
    elif model.upper() == "SEQENC":
        ids, results, pred_probabilities = run_seqenc.test_seqenc(fasta_file_input)
        if output_file == "":
            output_file = f"{fasta_file_input.split('.')[0]}_seqenc_prediction.txt"

    # Run the StrEnc model
    elif model.upper() == "STRENC":
        ids, results, pred_probabilities = run_strenc.test_strenc(struct_input)
        if output_file == "":
            output_file = f"{struct_input.split('.')[0]}_strenc_prediction.txt"

    # Write output to specified path
    output = open(output_file, "w")
    for id, pred, pred_probability in zip(ids, results, pred_probabilities):
        output.write(f"{id}\t{pred}\t{pred_probability}\n")
    output.close()
    print(f"\nResults are saved in {output_file}")

    return results

if __name__ == '__main__':

    print("\nWelcome to ncRNA classification\n"
          "Which model would you like to predict sequences saved in fasta file format?\n"
          "Merged provides the best results, but requires Graph Features predicted by GraphProt\n"
          "StrEnc provides the second best result, but requires Structural Encoding created by Pysster\n"
          "SeqEnc predicts using only the sequence, but does not require additional files\n"
          "GrEnc predicts using only Graph Encoding, but requires Graph Features predicted by GraphProt\n")

    predict_ncRNAs()
