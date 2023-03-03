import os
import data_processing
import sys
import pandas as pd
import numpy as np
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
import keras as ks
from sklearn.preprocessing import OneHotEncoder

pd.options.mode.chained_assignment = None  # default='warn'


def test_strenc(structure_input):

    # This method loads the trained Sequence model and tests it on sequences
    # sequence_df is a dataframe that has the sequence IDs as row names and at least one column named "Seq"
    # in which Sequences are stored as uppercase Strings including only IUPAC codes for nucleotides

    # Open the input file
    pysster_file = open(structure_input)
    lines = pysster_file.readlines()

    id_list = []
    structure_list = []

    # Variable i is used to determine what line we are in
    # using the modulo operation, we can determine, if the current line contains identifier, sequence or structure
    i = 0
    for line in lines:
        # Read in identifiers
        if i % 3 == 0:
            id_list.append(line.strip(">").split(" ")[0])
        # Save the most recent nucleotide sequence
        elif i % 3 == 1:
            sequence = line.strip("\n")
        # Read in structure sequences
        elif i % 3 == 2:
            struct_seq = (line.strip("\n"))
            # Annotate the structure sequence using struc_annotator()
            structure_list.append(data_processing.struc_annotator(sequence, struct_seq))
        i = i + 1

    pysster_file.close()

    structure_df = pd.DataFrame.from_dict({"id": id_list, "structure": structure_list})


    # Pad the structure sequences to the fixed length of 12,000
    structure_sequences = data_processing.pad_sequences(structure_df.structure, 12000)
    # Encode the structure sequence into integers
    structure_sequences = data_processing.struct_list_annotator(structure_sequences)
    # Transform the list of lists of integers into one matrix used as input for the model
    structure_sequences = data_processing.transform_seq_into_ml_input(structure_sequences)

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

    return id_list, results, pred_probabilities


if __name__ == '__main__':
    # Exception for when the command is not properly executed with structure file
    if len(sys.argv) != 2:
        print("Please enter a structure file when running the file\n"
              "Example:\n"
              "python run_strenc.py pysster_output/merged_test_file_30_pysster.txt")
    else:
        structure_input = sys.argv[1]

        # Load and test the model
        ids, results, pred_probabilities = test_strenc(structure_input)

        # Save output to file
        output_file = f"{structure_input.split('.')[0]}_strenc_predictions.txt"
        output = open(output_file, "w")
        for id, pred, pred_probability in zip(ids, results, pred_probabilities):
            output.write(f"{id}\t{pred}\t{pred_probability}\n")
        output.close()
        print(f"Results are saved in {output_file}")

    # change "from collections import Mapping" to "from collections.abc import Mapping" in linecloud.py
    # predict_structures("merged_test_file_30.fasta", "merged_test_file_30_pysster.txt", annotate=True)
