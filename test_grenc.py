import os
import data_processing
import sys
import pandas as pd
import numpy as np
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
import keras as ks
from sklearn.preprocessing import OneHotEncoder

pd.options.mode.chained_assignment = None  # default='warn'


def test_grenc(sequence_df, graph_input, fasta_file_input):

    # This method loads the trained GrEnc model and tests it on graph features
    # sequence_df is a dataframe that has the sequence IDs as row names and at least one column named "Seq"
    # in which Sequences are stored as uppercase Strings including only IUPAC codes for nucleotides
    # graph_input is the path to the file in which the corresponding graph features are saved
    # These graph features must match the sequences in sequence_df

    # Test if number of graph feature vector is the same as number of sequences
    if not data_processing.test_graphfeat_seq_match(sequence_df, graph_input):
        print("Graph Feature file does not match sequence file. Exiting.")
        sys.exit()

    # Read graph feature vector into the dataframe in a new column called "feature_vector"
    sequence_df = data_processing.read_graphprot_vectors(sequence_df,
                                                         graph_input,
                                                         fasta_file_input,
                                                         verbose=True)
    graph_input = np.array(sequence_df.feature_vector.to_list()).reshape(len(sequence_df.feature_vector),
                                                                         len(sequence_df.feature_vector[0]))
    # Drop "Seq" and "feature_vector" from data frame to save on memory
    sequence_df.drop(["Seq", "feature_vector"], axis=1, inplace=True)

    # Create a one hot encoder with the 6 possible RNA types to return the output as plain text
    rna_types = ["lncRNA", "miRNA", "rRNA", "snRNA", "snoRNA", "tRNA"]
    ohe = OneHotEncoder(sparse_output=False)
    ohe.fit(np.array(rna_types).reshape(-1, 1))

    # Load the model
    model = ks.models.load_model("model_files/grenc_fold2.hdf5")

    # Predict the ncRNA types from the two inputs
    prediction = model.predict(graph_input)
    # Return probability for each prediction
    pred_probabilities = [np.max(x) for x in prediction]
    # Convert results to ncRNA types
    results = ohe.inverse_transform(prediction)
    results = [x[0] for x in results]

    return results, pred_probabilities


if __name__ == '__main__':
    # Exception for when the command is not properly executed with fasta and feature file
    if len(sys.argv) != 3:
        print("Please enter a fasta file and a graph features file when running the file\n"
              "Example:\n"
              "python test_grenc.py merged_test_file_30.fasta graphprot_output/merged_test_file_30.gspan.gz.feature")
    else:
        fasta_file_input = sys.argv[1]
        graph_input = sys.argv[2]

        # Read in sequences as dataframe
        sequence_df = data_processing.read_fasta_file(fasta_file_input, labels=False)
        sequence_ids = sequence_df.index

        # Load and test the model
        results, pred_probabilities = test_grenc(sequence_df, graph_input, fasta_file_input)

        # Save output to file
        output_file = f"results/{fasta_file_input.split('.')[0]}_grenc_predictions.txt"
        output = open(output_file, "w")
        for id, pred, pred_probability in zip(sequence_ids, results, pred_probabilities):
            output.write(f"{id}\t{pred}\t{pred_probability}\n")
        output.close()
        print(f"Results are saved in {output_file}")
