import os
import data_processing
import sys
import pandas as pd
import numpy as np
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
import keras as ks
from sklearn.preprocessing import OneHotEncoder

pd.options.mode.chained_assignment = None  # default='warn'


def test_grenc(graph_input, fasta_file_input = ""):

    # This method loads the trained GrEnc model and tests it on graph features
    # sequence_df is a dataframe that has the sequence IDs as row names and at least one column named "Seq"
    # in which Sequences are stored as uppercase Strings including only IUPAC codes for nucleotides
    # graph_input is the path to the file in which the corresponding graph features are saved
    # These graph features must match the sequences in sequence_df

    if fasta_file_input != "":
        sequence_df = data_processing.read_fasta_file(fasta_file_input)
        if not data_processing.test_graphfeat_seq_match(sequence_df, graph_input):
            sys.exit()
            # Read graph feature vector into the dataframe in a new column called "feature_vector"
        sequence_df.drop(["Seq"], axis=1, inplace=True)
        sequence_df = data_processing.read_graphprot_vectors(sequence_df, graph_input)
    else:
        i = 0
        index_list = []
        for line in open(graph_input):
            index_list.append(f"sequence_{i}")
            i = i+1
        sequence_df = pd.DataFrame.from_dict({"id": index_list})
        sequence_df.index = index_list
        sequence_df = data_processing.read_graphprot_vectors(sequence_df,
                                                             graph_input,
                                                             verbose=False)



    # Test if number of graph feature vector is the same as number of sequences

    graph_matrix = np.array(sequence_df.feature_vectors.to_list()).reshape(len(sequence_df.feature_vectors),
                                                                         len(sequence_df.feature_vectors[0]))
    # Drop "Seq" and "feature_vector" from data frame to save on memory
    sequence_df.drop(["feature_vectors"], axis=1, inplace=True)

    # Create a one hot encoder with the 6 possible RNA types to return the output as plain text
    rna_types = ["lncRNA", "miRNA", "rRNA", "snRNA", "snoRNA", "tRNA"]
    ohe = OneHotEncoder(sparse_output=False)
    ohe.fit(np.array(rna_types).reshape(-1, 1))

    # Load the model
    model = ks.models.load_model("model_files/grenc_fold2.hdf5")

    # Predict the ncRNA types from the two inputs
    prediction = model.predict(graph_matrix)
    # Return probability for each prediction
    pred_probabilities = [np.max(x) for x in prediction]
    # Convert results to ncRNA types
    results = ohe.inverse_transform(prediction)
    results = [x[0] for x in results]

    ids = sequence_df.index

    return ids, results, pred_probabilities


if __name__ == '__main__':
    # Exception for when the command is not properly executed with fasta and feature file
    if len(sys.argv) not in (2, 3) :
        print("Please enter a graph features file when running the file\n"
              "Optionally, you may also enter the fasta file used for the creation of the feature file.\n"
              "If you do not provide the fasta file, the output rows will simply have 'sequence_n' as identifier.\n"
              "Run the file by providing the path(s) to the input(s) as run parameter(s)"
              "Example:\n"
              "python run_grenc.py graphprot_output/small_testset_30_graphprot.feature\n"
              "or\n"
              "python run_grenc.py graphprot_output/small_testset_30_graphprot.feature small_testset_30.fasta\n")
    else:
        graph_input = sys.argv[1]
        if len(sys.argv) == 3:
            fasta_file_input = sys.argv[2]
        else:
            fasta_file_input = ""

        # Read in sequences as dataframe

        # Load and test the model
        ids, results, pred_probabilities = test_grenc(graph_input, fasta_file_input)

        # Save output to file
        output_file = f"results/{fasta_file_input.split('.')[0]}_grenc_predictions.txt"
        output = open(output_file, "w")
        for id, pred, pred_probability in zip(ids, results, pred_probabilities):
            output.write(f"{id}\t{pred}\t{pred_probability}\n")
        output.close()
        print(f"Results are saved in {output_file}")
