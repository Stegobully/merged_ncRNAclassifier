import os
import data_processing
import sys
import pandas as pd
import numpy as np
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
import keras as ks
from sklearn.preprocessing import OneHotEncoder

pd.options.mode.chained_assignment = None  # default='warn'


def test_seqenc(fasta_file_input):

    # This method loads the trained Sequence model and tests it on sequences
    # sequence_df is a dataframe that has the sequence IDs as row names and at least one column named "Seq"
    # in which Sequences are stored as uppercase Strings including only IUPAC codes for nucleotides

    # Read in the sequences from the fasta file input
    sequence_df = data_processing.read_fasta_file(fasta_file_input)

    # Pad sequences to the fixed length of 12,000 nt by appending "_"
    sequence_input = data_processing.pad_sequences(sequence_df["Seq"], 12000)
    # Encode the nucleotides into integers
    sequence_input = data_processing.encode_nucleotides(sequence_input)
    # Transform the list of lists of integers into one matrix used as input for the model
    sequence_input = data_processing.transform_seq_into_ml_input(sequence_input)

    # Drop "Seq" and "feature_vector" from data frame to save on memory
    sequence_df.drop(["Seq"], axis=1, inplace=True)

    # Create a one hot encoder with the 6 possible RNA types to return the output as plain text
    rna_types = ["lncRNA", "miRNA", "rRNA", "snRNA", "snoRNA", "tRNA"]
    ohe = OneHotEncoder(sparse_output=False)
    ohe.fit(np.array(rna_types).reshape(-1, 1))

    # Load the model
    model = ks.models.load_model("model_files/seqenc_fold8.hdf5")

    # Predict the ncRNA types from the two inputs
    prediction = model.predict(sequence_input)
    # Return probability for each prediction
    pred_probabilities = [np.max(x) for x in prediction]
    # Convert results to ncRNA types
    results = ohe.inverse_transform(prediction)
    results = [x[0] for x in results]

    ids = sequence_df.index

    return ids, results, pred_probabilities


if __name__ == '__main__':
    # Exception for when the command is not properly executed with fasta and feature file
    if len(sys.argv) != 2:
        print("Please enter a fasta file when running the file\n"
              "Example:\n"
              "python run_seqenc.py small_testset_30.fasta")
    else:
        # Identify fasta file from run parameter
        fasta_file_input = sys.argv[1]

        # Load and test the model
        ids, results, pred_probabilities = test_seqenc(fasta_file_input)

        # Save output to file
        output_file = f"{fasta_file_input.split('.')[0]}_seqenc_predictions.txt"
        output = open(output_file, "w")
        for id, pred, pred_probability in zip(ids, results, pred_probabilities):
            output.write(f"{id}\t{pred}\t{pred_probability}\n")
        output.close()
        print(f"Results are saved in {output_file}")
