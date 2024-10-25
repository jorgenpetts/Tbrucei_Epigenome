import pandas as pd
from Bio import SeqIO

# # Path to the .bed file
# file_path = '/Users/jorgenrpettenburger/Desktop/Thesis research/SSRs_predicted.bed'

# # read in the .bed file
# data = pd.read_csv(file_path, sep='\t', header=None)

# print(data.head())

# # Set up a dictionary to store the data
# bed_dict = {
#     'chromosome': [],
#     'start': [],
#     'end': [],
#     'type': [],
#     'score': [],
#     'strand': [],
# }

# # populate the dictionary
# for index, row in data.iterrows():
#     bed_dict['chromosome'].append(row[0])
#     bed_dict['start'].append(row[1])
#     bed_dict['end'].append(row[2])
#     bed_dict['type'].append(row[3])
#     bed_dict['score'].append(row[4])
#     bed_dict['strand'].append(row[5])

# # Create a DataFrame from the dictionary
# bed_df = pd.DataFrame(bed_dict)

# bed_dict_new = {
#     'chromosome': [],
#     'start': [],
#     'end': [],
#     'type': [],
#     'score': [],
#     'strand': [],
# }


# # iterate over the rows in the DataFrame per chromosome
# for chromosome, group in bed_df.groupby('chromosome'):

#     # Get the chromosome data
#     chromosome_data = group

#     # Get the start and end positions
#     start_positions = chromosome_data['start'].values
#     end_positions = chromosome_data['end'].values

#     # Get the type
#     types = chromosome_data['type'].values

#     # Get the score
#     scores = chromosome_data['score'].values

#     # Get the strand
#     strands = chromosome_data['strand'].values

#     # Iterate over the rows in the chromosome data using index
#     i = 0
#     while i < (len(chromosome_data) - 1):
#         start = chromosome_data.iloc[i]['start']
#         end = chromosome_data.iloc[i]['end']
#         type = chromosome_data.iloc[i]['type']
#         score = chromosome_data.iloc[i]['score']
#         strand = chromosome_data.iloc[i]['strand']

#         start_next = chromosome_data.iloc[i + 1]['start']
#         end_next = chromosome_data.iloc[i + 1]['end']
#         type_next = chromosome_data.iloc[i + 1]['type']
#         score_next = chromosome_data.iloc[i + 1]['score']
#         strand_next = chromosome_data.iloc[i + 1]['strand']

#         if type == type_next and strand == '-' and strand_next == '+':
#             # Update the dictionary
#             bed_dict_new['chromosome'].append(chromosome)
#             bed_dict_new['start'].append(start)
#             bed_dict_new['end'].append(start_next)
#             bed_dict_new['type'].append(type)
#             bed_dict_new['score'].append(score)
#             bed_dict_new['strand'].append('.')
#             i += 2
#         else:
#             i += 1
            

# # Create a DataFrame from the new dictionary
# bed_df_new = pd.DataFrame(bed_dict_new)

# write new bed file
# bed_df_new.to_csv('SSRs_predicted_new.bed', sep='\t', header=False, index=False)

# TODO: cut off 10% of sequence from the start and end of the chromosome
# TODO: Fix anomalies in the .bed file by hand (e.g. overlapping SSRs)




fasta_file = '/Users/jorgenrpettenburger/Desktop/Thesis research/Code/T.brucei_427_ref.fasta'

# Read the FASTA file using Biopython
fasta_sequences = SeqIO.parse(open(fasta_file), 'fasta')

# get the lengths of each chromosome
chromosome_lengths = {}
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    chromosome_lengths[name] = len(sequence)

# read SSRs_predicted_new.bed file
ssr_file = '/Users/jorgenrpettenburger/Desktop/Thesis research/Code/SSRs_predicted_new.bed'

ssr_data = pd.read_csv(ssr_file, sep='\t', header=None)

# Set up a dictionary to store the data
ssr_dict = {
    'chromosome': [],
    'start': [],
    'end': [],
    'type': [],
    'score': [],
    'strand': [],
}

# populate the dictionary
for index, row in ssr_data.iterrows():
    ssr_dict['chromosome'].append(row[0])
    ssr_dict['start'].append(row[1])
    ssr_dict['end'].append(row[2])
    ssr_dict['type'].append(row[3])
    ssr_dict['score'].append(row[4])
    ssr_dict['strand'].append(row[5])

# Create a DataFrame from the dictionary
ssr_df = pd.DataFrame(ssr_dict)

ssr_dict_new = {
    'chromosome': [],
    'start': [],
    'end': [],
    'type': [],
    'score': [],
    'strand': [],
}

# iterate over the rows in the DataFrame per chromosome
for chromosome, group in ssr_df.groupby('chromosome'):

    # Get the chromosome data
    chromosome_data = group

    # Get the start and end positions
    start_positions = chromosome_data['start'].values
    end_positions = chromosome_data['end'].values

    # Get the type
    types = chromosome_data['type'].values

    # Get the score
    scores = chromosome_data['score'].values

    # Get the strand
    strands = chromosome_data['strand'].values

    # Iterate over the rows in the chromosome data using index
    i = 0
    while i < (len(chromosome_data)):
        start = chromosome_data.iloc[i]['start']
        end = chromosome_data.iloc[i]['end']
        type = chromosome_data.iloc[i]['type']
        score = chromosome_data.iloc[i]['score']
        strand = chromosome_data.iloc[i]['strand']

        # cut off 10% of sequence from the start and end of the chromosome

        if start >= 0.1 * chromosome_lengths[chromosome] and end <= 0.9 * chromosome_lengths[chromosome]:
            ssr_dict_new['chromosome'].append(chromosome)
            ssr_dict_new['start'].append(start)
            ssr_dict_new['end'].append(end)
            ssr_dict_new['type'].append(type)
            ssr_dict_new['score'].append(score)
            ssr_dict_new['strand'].append(strand)
            
        i += 1


# Create a DataFrame from the new dictionary
ssr_df_new = pd.DataFrame(ssr_dict_new)

# write new bed file
ssr_df_new.to_csv('SSRs_predicted_new_telomeric_cutoff.bed', sep='\t', header=False, index=False)