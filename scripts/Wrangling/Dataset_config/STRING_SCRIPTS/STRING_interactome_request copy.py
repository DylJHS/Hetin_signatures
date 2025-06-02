
import requests
import json
import re
import pandas as pd

pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)


string_api_url = "https://version-11-5.string-db.org/api"
output_format = "tsv-no-header"

##LEM Domain proteins, LMNA, LMNB1 & LBR
protein_genes_of_interest = [
    "MAN1","TMPO", "EMD", "TOR1AIP1", 
    'IFFO1', "LMNA", "LMNB1", "LMNB2", 
    "LBR", "LEMD2", "PRR14", "CBX5",
    "CBX1","CBX3", "LEMD3", "EMD",
	"BANF1","BAF","TERB1","TERB2",
	"MAJIN"
    ]
    

## Set parameters

params = {

    "identifiers" : "\r".join(protein_genes_of_interest), # your protein list
    "species" : 9606, # species NCBI identifier 
    "limit" : 1, # only one (best) identifier per input protein
    "echo_query" : 1, # see your input identifiers in the output
    "caller_identity" : "www.awesome_app.org" }

## Construct URL
request_url = "/".join([string_api_url, output_format, 'get_string_ids'])

## Call STRING
results = requests.post(request_url, data=params)


## Read and parse the results
identifiers = []
for line in results.text.strip().split("\n"):
    l = line.split("\t")
    input_identifier, string_identifier = l[0], l[2]
    identifiers.append(string_identifier)
# print(identifiers)


## Set new parameters
new_params = {

    "identifiers" : "%0d".join(identifiers), # your protein list
    "species" : 9606, # species NCBI identifier 
    "limit": 1500,
    "network_type": "physical",
    "caller_identity" : "www.awesome_app.org" }


request_url2 = "/".join([string_api_url, output_format, "interaction_partners"])

response2 = requests.post(request_url2, data=new_params)


df = pd.DataFrame(columns = [
    'Lamin_Protein', 'Interactor',
    'Lamin_Protein_ENSP', 'Interacting_Protein_ENSP',
    'Combined_Score', 'Experimental_Score'
    ]
)

for line in response2.text.strip().split("\n"):


    l = line.strip().split("\t")
    query_ensp = l[0]
    query_name = l[2]
    partner_ensp = l[1]
    partner_name = l[3]
    combined_score = l[5]
    Experimental_score = l[10]

    protein_genes_of_interest.append(partner_name)

    row = [
        [
        query_name, partner_name, 
        query_ensp, partner_ensp, 
        combined_score, Experimental_score
        ]
    ]

    rowed = pd.DataFrame(
        data = row, 
        columns = [
            'Lamin_Protein', 'Interactor','Lamin_Protein_ENSP', 
            'Interacting_Protein_ENSP','Combined_Score', 'Experimental_Score'
        ]
    )

    df = pd.concat([df, rowed]).reset_index(drop = True)


full_list = list(set(protein_genes_of_interest))
print(' '.join(full_list))

print(len(full_list))





