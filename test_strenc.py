import os
import data_processing
import sys
import pandas as pd
import numpy as np
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
import keras as ks
from sklearn.preprocessing import OneHotEncoder

pd.options.mode.chained_assignment = None  # default='warn'


def test_strenc(sequence_df, structure_input):

    # This method loads the trained Sequence model and tests it on sequences
    # sequence_df is a dataframe that has the sequence IDs as row names and at least one column named "Seq"
    # in which Sequences are stored as uppercase Strings including only IUPAC codes for nucleotides

    pysster_file = open(structure_input)
    lines = pysster_file.readlines()

    seq_id = ""
    nt_seq = ""
    i = 1
    sequence_df["structure_sequence"] = ""

    for line in lines:
        if line.startswith(">"):
            seq_id = line.split(" ")[0].strip(">")
        elif i % 3 == 2:
            nt_seq = line.strip("\n")
        elif i % 3 == 0:
            struct_seq = line.strip("\n")
            sequence_df["structure_sequence"][seq_id] = data_processing.struc_annotator(nt_seq, struct_seq)

        i = i + 1

    # Pad the structure sequences to the fixed length of 12,000
    structure_sequences = data_processing.pad_sequences(sequence_df["structure_sequence"], 12000)
    # Encode the structure sequence into integers
    structure_sequences = data_processing.struct_list_annotator(structure_sequences)
    # Transform the list of lists of integers into one matrix used as input for the model
    structure_sequences = data_processing.transform_seq_into_ml_input(structure_sequences)

    # Drop "Seq" and "structure_sequence" from data frame to save on memory
    sequence_df.drop(["Seq", "structure_sequence"], axis=1, inplace=True)

    # Create a one hot encoder with the 6 possible RNA types to return the output as plain text
    rna_types = ["lncRNA", "miRNA", "rRNA", "snRNA", "snoRNA", "tRNA"]
    ohe = OneHotEncoder(sparse_output=False)
    ohe.fit(np.array(rna_types).reshape(-1, 1))

    # Load the model
    model = ks.models.load_model("model_files/strenc_fold7.hdf5")

    # Predict the ncRNA types from the two inputs
    prediction = model.predict(structure_sequences)
    # Return probability for each prediction
    pred_probabilities = [np.max(x) for x in prediction]
    # Convert results to ncRNA types
    results = ohe.inverse_transform(prediction)
    results = [x[0] for x in results]

    return results, pred_probabilities


if __name__ == '__main__':
    # Exception for when the command is not properly executed with fasta and feature file
    if len(sys.argv) != 3:
        print("Please enter a fasta and a structure file when running the file\n"
              "Example:\n"
              "python test_strenc.py merged_test_file_30.fasta pysster_output/merged_test_file_30_pysster.txt")
    else:
        fasta_file_input = sys.argv[1]
        structure_input = sys.argv[2]
        # Read in sequences as dataframe
        sequence_df = data_processing.read_fasta_file(fasta_file_input, labels=False)
        sequence_ids = sequence_df.index

        # Load and test the model
        results, pred_probabilities = test_strenc(sequence_df, structure_input)

        # Save output to file
        output_file = f"results/{fasta_file_input.split('.')[0]}_strenc_predictions.txt"
        output = open(output_file, "w")
        for id, pred, pred_probability in zip(sequence_ids, results, pred_probabilities):
            output.write(f"{id}\t{pred}\t{pred_probability}\n")
        output.close()
        print(f"Results are saved in {output_file}")

    # change "from collections import Mapping" to "from collections.abc import Mapping" in linecloud.py
    # predict_structures("merged_test_file_30.fasta", "merged_test_file_30_pysster.txt", annotate=True)
