import os
import json
import numpy as np

# Set the new working directory for UniProt data
new_dir_uniprot = 'C:/Users/riocx/Documents/Masters new/Summer 2023/Dr. Zhao/Mapping sequences/For annotations'
os.chdir(new_dir_uniprot)  # Change the current working directory to new_dir_uniprot

# Load the JSON data from the UniProt file
with open('uniprot.json', 'r') as file:  # Open the file in read mode
    data_uniprot = json.load(file)  # Load the JSON data from the file into the variable data_uniprot

# Create a list to store the extracted UniProt information
output_data_uniprot = []  # Initialize an empty list to store processed UniProt data

# Iterate over the proteins in the UniProt dataset
for resultsEntry in data_uniprot['results']:  # Loop through each entry in the 'results' list of the JSON data
    # Initialize variables for each entry
    protein_regions = {}  # Dictionary to store protein region information
    geneNames = []  # List to store gene names
    type = ''  # Variable to store the feature type
    description = ''  # Variable to store the feature description
    uniProtkb_Id = ''  # Variable to store the UniProtKB ID
    uniprotID = ''  # Variable to store the UniProt ID
    scientificName = ''  # Variable to store the scientific name of the organism
    taxonId = ''  # Variable to store the taxon ID
    sequence = ''  # Variable to store the sequence

    # Iterate over the features of the protein
    for featureEntry in resultsEntry.get('features', []):  # Loop through each feature in the 'features' list
        feature_type = featureEntry.get('type')  # Get the type of the feature
        # Check if the description contains specific terms
        if any(term in featureEntry.get('description', '') for term in ["linker", "binding", "Disordered"]):
            description = featureEntry.get('description')  # Get the description of the feature
            start = featureEntry['location']['start']['value']  # Get the start position of the feature
            end = featureEntry['location']['end']['value']  # Get the end position of the feature
            type = featureEntry['type']  # Get the type of the feature

            # Add feature information to the protein_regions dictionary
            if feature_type not in protein_regions:  # Check if the feature type is already in the dictionary
                protein_regions[feature_type] = []  # If not, add it with an empty list
            protein_regions[feature_type].append({'start': start, 'end': end})  # Append the start and end positions

            # Extract general protein information
            uniProtkb_Id = resultsEntry['uniProtkbId']  # Get the UniProtKB ID
            uniprotID = resultsEntry['primaryAccession']  # Get the primary accession ID
            scientificName = resultsEntry['organism']['scientificName']  # Get the scientific name of the organism
            taxonId = resultsEntry['organism']['taxonId']  # Get the taxon ID
            sequence = resultsEntry['sequence']['value']  # Get the sequence value

    # Extract recommended name and alternative names
    proteinDescription = resultsEntry.get('proteinDescription', {})  # Get the protein description dictionary
    recommendedName = proteinDescription.get('recommendedName', {}).get('fullName', {}).get('value', '')  # Get the recommended name
    alternativeNames = [name.get('fullName', {}).get('value', '') for name in proteinDescription.get('alternativeNames', [])]  # Get the alternative names

    # Extract gene names
    genes = resultsEntry.get('genes', [])  # Get the list of genes
    geneNames = [gene.get('geneName', {}).get('value', '') for gene in genes]  # Get the gene names

    # Create a dictionary with the extracted data
    extracted_data_uniprot = {
        "uniProtkb Id": uniProtkb_Id,  # Add the UniProtKB ID to the dictionary
        "uniprot ID": uniprotID,  # Add the UniProt ID to the dictionary
        "scientificName": scientificName,  # Add the scientific name to the dictionary
        "taxonId": taxonId,  # Add the taxon ID to the dictionary
        "Regions": protein_regions,  # Add the protein regions to the dictionary
        "Type": type,  # Add the feature type to the dictionary
        "Description": description,  # Add the description to the dictionary
        "Recommended Names": recommendedName,  # Add the recommended names to the dictionary
        "Alternative Names": alternativeNames,  # Add the alternative names to the dictionary
        "Gene Names": geneNames,  # Add the gene names to the dictionary
        "Sequence": sequence  # Add the sequence to the dictionary
    }

    # Append the extracted data to the output list
    output_data_uniprot.append(extracted_data_uniprot)  # Add the extracted data dictionary to the output list

# Write the output data to a JSON file
output_path_uniprot = "output_with_sequences.json"  # Define the output file path
with open(output_path_uniprot, 'w') as file:  # Open the file in write mode
    json.dump(output_data_uniprot, file, indent=4)  # Write the output data to the file with indentation for readability

print("UniProt data is saved to:", output_path_uniprot)  # Print a message indicating the output file path for UniProt data

# Set the new working directory for DisProt data
new_dir_disprot = 'C:/Users/riocx/Documents/Masters new/Summer 2023/Dr. Zhao/Mapping sequences/For annotations/Disprot'
os.chdir(new_dir_disprot)  # Change the current working directory to new_dir_disprot

# Specify the path to the DisProt JSON file and the encoding
json_file_path = "DisProt release_2023_06 with_ambiguous_evidences.json"  # Define the JSON file path
encoding = 'utf-8'  # Use 'utf-8' encoding for the JSON file

# Read data from the JSON file with the specified encoding
with open(json_file_path, "rt", encoding=encoding) as file:  # Open the file in read mode with specified encoding
    data_disprot = json.load(file)  # Load the JSON data from the file into the variable data_disprot

# Extract the relevant data from the JSON
b = data_disprot["data"]  # Extract the 'data' list from the JSON data

# Open a text file for writing the DisProt disorder annotations
with open("DisProt_disorder_annotation_ambiguous.txt", "wt") as y:  # Open the file in write mode
    lst = []  # Initialize an empty list to store the annotations

    # Iterate over each entry in the DisProt data
    for j in b:  # Loop through each entry in the 'data' list
        uniprot = j["acc"]  # Get the UniProt accession
        disid = j["disprot_id"]  # Get the DisProt ID
        seq = j["sequence"]  # Get the amino acid sequence
        dis = np.zeros(len(seq), dtype=np.uint8)  # Initialize a numpy array for disorder annotations
        disorder_consensus = j["disprot_consensus"]  # Get the disorder consensus information
        disorder = disorder_consensus["Structural state"]  # Get the structural state information

        # Annotate the sequence with disorder information
        for n in disorder:  # Loop through each disorder annotation
            if n["type"] == "D":  # Check if the type is disordered
                dis[n["start"]-1:n["end"]] = 1  # Mark disordered regions with 1
            elif n["type"] == "S":  # Check if the type is structured
                dis[n["start"]-1:n["end"]] = 2  # Mark structured regions with 2

        # Append the formatted annotation to the list
        lst.append(uniprot + "_" + disid + '\n' + seq.strip() + "\n" + "".join(list(map(str, dis))) + "\n")

    # Write the annotations to the text file
    for n in lst:  # Loop through each annotation in the list
        y.write(">" + n.strip() + "\n")  # Write the formatted annotation to the file

print("DisProt data is saved to DisProt_disorder_annotation_ambiguous.txt")  # Print a message indicating the output file path for DisProt data

