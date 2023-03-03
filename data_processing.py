import sys
import os
import time
import numpy as np
import pandas as pd
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
import tensorflow as tf
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from scipy import sparse
from sklearn.preprocessing import OrdinalEncoder

pd.options.mode.chained_assignment = None  # default='warn'

########################################################################################################################

def read_fasta_file(filename):
    # This function reads in a single fasta file into a dataframe
    # The dataframe will contain one column "Seq" for the sequence
    # as well as columns for the rna type, rna subtype and sequence length
    # To read in all fastas in a folder refer to read_fastas
    # "labels" allows the user to add the rna type
    # to each fasta header and read them into the df as a column


    sequence_dict = {rec.id: [str(rec.seq).upper()] for rec in SeqIO.parse(filename, "fasta")}
    df = pd.DataFrame.from_dict(sequence_dict, orient="index", columns=["Seq"])

    df["length"] = df["Seq"].map(len)

    return df

########################################################################################################################


def transform_seq_into_graphfeatures(fasta_file,
                                     path=os.getcwd(),
                                     fasta2shrep_path="~/tools/pkgs/graphprot-1.1.7-2/libexec/graphprot/fasta2shrep_gspan.pl",
                                     eden_path="~/tools/envs/graphprot/libexec/graphprot/EDeN"):
    
    # This function applies the graphprot transformation to a fasta file
    # The first step is to create a graph representation of the secondary structure using fasta2shrep
    # The second step is then to create the feature vectors using EDeN
    # This function takes a very long time for longer sequences
    # The user needs to provide the path to the fasta2shrep.pl executable
    # as well as the path to EDeN
    
    if not os.path.isdir(path):
        print(path)
        print("Please provide a valid directory for the output")
    else:
        # Create path for graphprot outputs
        if not os.path.isdir(f"{path}/graphprot_output"):
            os.mkdir(f"{path}/graphprot_output")

        start_time = time.time()

        gspan_file = f"{path}/graphprot_output/{fasta_file.split('.')[0]}.gspan.gz"
        # Create graph representation of the secondary structure
        os.system(f"""{fasta2shrep_path}
         -abstr -stdout -M 3 -wins '150,' -shift '25' -fasta {fasta_file} -t 3 | gzip > {gspan_file}""")
        # Use the Graph Kernel method to create the feature vectors
        os.system(f"{eden_path} -a FEATURE -i {gspan_file}")

        end_time = time.time()
        print(f"Time elapsed: {(end_time-start_time)/(60*60)}h")
                    
########################################################################################################################


def read_graphprot_vectors(df, feature_file):
    
    # This function reads in the feature vectors created by transform_seq_into_graphfeatures
    # You have to provide a dataframe with the corresponding sequence identifiers
    # The order within the dataframe has to be the same as the rows in the feature file
    # The graphprot feature vectors are read in as sparse np vectors
    # verbose indicates, whether the name of the currently processed file is printed out

    # Create empty column in which the feature vectors are written in
    df["feature_vector"] = None
    if feature_file.endswith(".feature"):
        try:
            feature_file = open(feature_file, "r")
        except IOError:
            print(f"Could not open {feature_file}")
            sys.exit()

        lines = feature_file.readlines()
        feature_file.close()

        df["feature_vectors"] = ""
        # Every line is saved as a vector where each entry is one feature
        for line, i in zip(lines, range(len(lines))):
            # Create an empty vector to fill with the features
            unsparse_vector = np.zeros(32768)
            for feature in line.split(" "):
                if feature != "":
                    # Every feature is of the form "position:value", where position is the index in the vector
                    # and value is the value for the feature
                    unsparse_vector[int(feature.split(":")[0])] = feature.split(":")[1].strip("\n")
            # Save the vector in the dataframe
            df.feature_vectors[i] = unsparse_vector

    return df

########################################################################################################################


def pad_sequences(seq_list, length, char="_"):
    
    # This function pads a list of sequences of nucleotides to a specified length
    # seq_list is a column of a dataframe, length is the desired padding length
    # left indicates, on which side the sequences are padded (True => padded on the left)
    # char is the padding character
    # For sequences that are longer than the specified length,

    # Pad sequences on the right with the provided character
    pad_list = seq_list.str.ljust(length, char)
    # Reduce length of sequences that are to long
    pad_list = pad_list.map(lambda x: x[0:length])

    return pad_list

########################################################################################################################


def test_graphfeat_seq_match(seq_df, graph_file):

    # This function tests if the sequence df has the same amount of sequences as the gspan.gz.feature file has
    # graph feature vectors
    # It returns True if the numbers match and False else

    num_lines = sum(1 for line in open(graph_file))
    if num_lines != len(seq_df):
        print(f"Number of lines in {graph_file} does not match number of sequences")
        return False
    else:
        return True

########################################################################################################################


def encode_nucleotides(seq_list):
    
    # This function encodes nucleotides into integer values
    # Allowed characters include the IUPAC ambiguity codes, as well as the padding character "_"
    # Unknown characters will be assigned the value 16
    # Input is a column of a dataframe of padded sequences
    # Output is a list of lists of integers encoding the sequence
    
    categories = ["A", "C", "G", "T", "N", "R", "K", "S", "Y", "M", "W", "B", "H", "D", "V", "_"]
    # create scikit-learns ordinal encoder
    ordi = OrdinalEncoder(handle_unknown="use_encoded_value", unknown_value=16)
    ordi.fit(np.array(list(categories)).reshape(-1, 1))
    
    return seq_list.map(lambda seq: ordi.transform(np.array(list(seq)).reshape(-1, 1)))
    
##########################################################################################################


def transform_seq_into_ml_input(seq_list):
    
    # This function takes the encoded sequences created by encode_nucleotides
    # and transforms them into an array used as input for the sequence CNN model
    
    return np.array(seq_list.to_list()).reshape((len(seq_list), len(seq_list[0])))

##########################################################################################################


def write_seq_df_to_fasta(df,
                          seq_col="Seq",
                          type_col="rna_type",
                          filename="Sequence_file.fasta",
                          outpath=os.getcwd()):
    
    # This function takes a dataframe of the format produced by read_fastas 
    # and writes the sequences into new fasta files
    # Fasta files are saved to the file provided in outpath

    seq_list = []
    for i in range(0, len(df)):
        record = SeqRecord(Seq(df[seq_col][i]), df.index[i], "", df[type_col][i])
        seq_list.append(record)
    fasta_file = f"{outpath}/{filename}"
    if not os.path.isfile(fasta_file):
        SeqIO.write(seq_list, fasta_file, "fasta")
    else:
        print(f"{fasta_file} already exists. Please choose a different outpath")

###########################################################################################################


def read_ncr_results(path, df, col_name, tabs=True):

    # This method can be used to read in the results created by predicting with ncRDense or ncRDeep
    # The results are read into a column of a data frame, the name of the column is provided by col_name
    # Keys of the data frame must match the sequence IDs
    # tabs is a bool that indicates whether the result file is separated by tabs or commas
    
    df[col_name] = ""
    for filename in os.listdir(path):
        print(filename)
        if filename.endswith(".txt"):
            with open(path + "/" + filename) as file:
                for line in file:
                    if tabs:
                        seq_id = line.split("\t")[0].strip()
                        pred = line.split("\t")[1].strip()
                    else:
                        seq_id = line.split(",")[0].strip()
                        pred = line.split(",")[1].strip()
                    try:
                        df[col_name][seq_id] = pred
                    except:
                        print(f"{seq_id} not found in df")

########################################################################################################################


def struc_annotator(sequence, structure):

    # This method gets a nucleotide sequence and a structure sequence created by pysster
    # The structure sequence can only contain "F", "S", "I", "M", "H", "T"
    # while the nucleotide sequence may include all IUPAC characters
    # Annotation is done according to supplementary table 3

    if len(sequence) != len(structure):
        print("Sequence and structure not of equal length")
    else:
        annotated_struct = ""

        # Assignment of letters for combinations
        dic = {"A": ["Q", "W", "E", "R", "T", "Z"],
               "C": ["U", "I", "O", "P", "A", "S"],
               "G": ["D", "F", "G", "H", "J", "K"],
               "T": ["L", "Y", "X", "C", "V", "B"],
               "other": ["N", "N", "N", "N", "N", "N"]}
        df = pd.DataFrame(dic, index=["F", "S", "I", "M", "H", "T"])
        for i in range(0, len(sequence)):
            if sequence[i] in ["A", "C", "G", "T"]:
                # Prolong structure sequence by new combination
                annotated_struct = annotated_struct + df[sequence[i]][structure[i]]
            else:
                # If any IUPAC code other than the four nts is used, the structure sequence is prolonged with "N"
                annotated_struct = annotated_struct + df["other"][structure[i]]
        return annotated_struct

########################################################################################################################


def struct_list_annotator(seq_list):

    # This method encodes the possible characters in the pysster structure file using sklearns ordinal encoder

    categories = ["Q", "W", "E", "R", "T", "Z", "U", "I", "O", "P",
                  "A", "S", "D", "F", "G", "H", "J", "K", "L", "Y",
                  "X", "C", "V", "B", "N", "_"]
    ordi = OrdinalEncoder(handle_unknown="use_encoded_value", unknown_value=len(categories))
    ordi.fit(np.array(list(categories)).reshape(-1, 1))
    return seq_list.map(lambda seq: ordi.transform(np.array(list(seq)).reshape(-1, 1)))
