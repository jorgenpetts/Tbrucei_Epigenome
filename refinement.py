import pandas as pd
#import biopython

# Path to the .xls file
file_path = '/Users/jorgenrpettenburger/Desktop/Thesis research/SSRs.xls'

# Read the .xls file
data = pd.read_excel(file_path)

# Initialize an empty dictionary with specified keys
data_dict = {
    'chromosome': [],
    'First ORF': [],
    'Last ORF': [],
    'Strand': [],
    'Type': []
}


# Create a DataFrame from the dictionary


# Loop through the data and assign values to the dictionary
for index, row in data.iterrows():
    data_dict['chromosome'].append(row['Chromosome'])
    data_dict['First ORF'].append(row['First annotated ORF in unit'])
    data_dict['Last ORF'].append(row['Last annotated ORF in unit'])
    data_dict['Strand'].append(row['Strand'])
    data_dict['Type'].append(row['Type'])

df = pd.DataFrame(data_dict)

    # Print the dictionary




# Change the value in the dataframe under 'chromosome' from 927 to 427


# Get a list of all unique names under 'chromosome'
unique_chromosomes = df['chromosome'].unique().tolist()

print(unique_chromosomes)


for i in range(len(unique_chromosomes)):
    temp = unique_chromosomes[i].split('_') 
    name = "Tb427"
    for j in range(1, len(temp)):
        name += "_" + temp[j]
    df['chromosome'] = df['chromosome'].replace(unique_chromosomes[i], name)

unique_chromosomes = df['chromosome'].unique().tolist()
print(unique_chromosomes)

# Print the updated dataframe
print(df)

unique_first_orf = df['First ORF'].unique().tolist()

for i in range(len(unique_first_orf)):
    temp = unique_first_orf[i].split(".")
    name = f"Tb427.0{temp[1]}.{temp[2]}"
    df['First ORF'] = df['First ORF'].replace(unique_first_orf[i], name)
    
unique_last_orf = df['Last ORF'].unique().tolist()

for i in range(len(unique_last_orf)):
    temp = unique_last_orf[i].split(".")
    name = f"Tb427.0{temp[1]}.{temp[2]}"
    df['Last ORF'] = df['Last ORF'].replace(unique_last_orf[i], name)

for index, row in df.iterrows():
    chromosome = row['chromosome']
    first_orf = row['First ORF']
    last_orf = row['Last ORF']
    if chromosome == "Tb427_09_v4":
        temp = first_orf.split(".")
        name = f"Tb427tmp.{(temp[1])[1:]}.{temp[2]}" 
        df.at[index, 'First ORF'] = name
        temp = last_orf.split(".")
        name = f"Tb427tmp.{(temp[1])[1:]}.{temp[2]}"
        df.at[index, 'Last ORF'] = name
    if chromosome == "Tb427_06_v4":
        if first_orf == "Tb427.03A7.960":
            df.at[index, 'First ORF'] = "Tb427tmp.3A7.960"
            df.at[index, 'Last ORF'] = "Tb427tmp.3A7.960"
    if chromosome == "Tb427_07_v4":
        if last_orf == "Tb427.030D13.80":
            df.at[index, 'Last ORF'] = "Tb427tmp.30D13.80"
    if chromosome == "Tb427_10_v4":
        temp = first_orf.split(".")
        name = f"Tb427.10.{temp[2]}"
        df.at[index, 'First ORF'] = name
        temp = last_orf.split(".")
        name = f"Tb427.10.{temp[2]}"
        df.at[index, 'Last ORF'] = name

    if chromosome == "Tb427_11_01_v4":
        temp = first_orf.split(".")
        name = f"Tb427tmp.{(temp[1])[1:]}.{temp[2]}"
        df.at[index, 'First ORF'] = name
        temp = last_orf.split(".")
        name = f"Tb427tmp.{(temp[1])[1:]}.{temp[2]}"
        df.at[index, 'Last ORF'] = name


#df.to_excel('/Users/jorgenrpettenburger/Desktop/Thesis research/SSRsUpdated.xlsx', sheet_name='Sheet1',index = False)

# read from a gff file
gff_file = '/Users/jorgenrpettenburger/Desktop/Thesis research/TriTrypDB-68_TbruceiLister427.gff'
gff_data = pd.read_csv(gff_file, sep = '\t', comment = '#', header = None)

# Set up a dictionary to store the data
gff_dict = {
    'chromosome': [],
    'source': [],
    'type': [],
    'start': [],
    'end': [],
    'score': [],
    'strand': [],
    'phase': [],
    'attributes': []
}

# populate the dictionary
for index, row in gff_data.iterrows():
    gff_dict['chromosome'].append(row[0])
    gff_dict['source'].append(row[1])
    gff_dict['type'].append(row[2])
    gff_dict['start'].append(row[3])
    gff_dict['end'].append(row[4])
    gff_dict['score'].append(row[5])
    gff_dict['strand'].append(row[6])
    gff_dict['phase'].append(row[7])
    gff_dict['attributes'].append(row[8])

# Create a DataFrame from the dictionary
gff_df = pd.DataFrame(gff_dict)

bed_dict = {
    'chromosome': [],
    'start': [],
    'end': [],
    'name': [],
    'score': [],
    'strand': [],

}

print(gff_df)

print(df)

# Step 1: Pre-process 'gff_df' to extract 'id' and 'parent' columns from 'attributes'
gff_df[['id', 'parent']] = gff_df['attributes'].str.split(';', expand=True)[[0, 1]]
gff_df['id'] = gff_df['id'].str.split('=').str[1]
gff_df['parent'] = gff_df['parent'].str.split('=').str[1]

# replace chromosome Tb427_10_v4 with Tb427_10_v5 in df
df['chromosome'] = df['chromosome'].replace('Tb427_10_v4', 'Tb427_10_v5')

for index, row in df.iterrows():
    chromosome = row['chromosome']
    first_orf = row['First ORF']
    last_orf = row['Last ORF']
    strand = row['Strand']
    type = row['Type']
    filtered_gff = gff_df[(gff_df['chromosome'] == chromosome) & 
                          ((gff_df['id'] == first_orf) | (gff_df['parent'] == first_orf))]

        

    if not filtered_gff.empty:

        if type == 'Divergent':
            name = "dSSR"
        elif type == 'Internal':
            name = "iSSR"
        else:
            name = type

        row2 = filtered_gff.iloc[0]
        bed_dict['chromosome'].append(row2[0])
        bed_dict['start'].append(row2[3])
        bed_dict['end'].append(row2[4])
        bed_dict['name'].append(name)
        bed_dict['score'].append(row2[5])
        bed_dict['strand'].append(row2[6])
        
        


bed_df = pd.DataFrame(bed_dict)

bed_df.to_csv('/Users/jorgenrpettenburger/Desktop/Thesis research/SSRs_predicted.bed', sep = '\t', index = False, header = False)

print(bed_df)