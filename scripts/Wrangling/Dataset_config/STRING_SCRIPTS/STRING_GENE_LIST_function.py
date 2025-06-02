
import requests, json, re
import pandas as pd

pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)


def Gene_list(protein_genes):
    """
    Queries the STRING PPI database to retrieve a list of genes interacting with a given list of proteins.

    This function takes a list of proteins of interest, converts them to their corresponding STRING database identifiers,
    and then uses these identifiers to query the database for interacting genes. The output is a string containing
    the full list of unique interacting proteins.

    The function performs two main steps:
    1. Conversion of input protein names to STRING identifiers.
    2. Retrieval of physically interacting genes based on these identifiers.

    Parameters:
    protein_genes_of_interest (list of str): A list of protein names to query in the STRING database.

    Returns:
    list: A list containing unique names of physically interacting genes, derived from the initial list of proteins of interest.

    Note:
    The function assumes access to the STRING API and requires an internet connection to fetch data.
    The returned list of interacting genes is in a comma-separated list format.
    """

    # Base URL and output format for STRING API

    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "tsv-no-header"
        

    # Parameters for the first API request to convert protein names to STRING IDs

    params = {
        "identifiers" : "\r".join(protein_genes), # list of initial genes
        "species" : 9606, # Human NCBI identifier 
        "limit" : 1, # only one (best) identifier per input protein
        "echo_query" : 1, # include input identifiers in the output
        "caller_identity" : "www.awesome_app.org" # identity of the caller
    }


    # Construct the request URL for getting STRING IDs

    request_url = "/".join([string_api_url, output_format, 'get_string_ids'])


    # Call STRING for the list of IDs corresonding to those
    # of the genes in the proteins of interest list

    protein_id_results = requests.post(request_url, data=params)
    # print(protein_id_results)


    # Call STRING for the list of IDs corresonding to those
    # of the genes in the proteins of interest list

    if protein_id_results.status_code == 521 or protein_id_results.status_code == 400:
        print('STRING WebServer is down')
        return


    # Read and parse the results
    identifiers = []
    for line in protein_id_results.text.strip().split("\n"):
        new_line = line.split("\t")
        input_identifier, string_identifier = new_line[0], new_line[2]
        identifiers.append(string_identifier)


    # Parameters for the second API request to get interaction partners

    new_params = {
        "identifiers" : "%0d".join(identifiers), # list of STRING IDs
        "species" : 9606, # Human NCBI identifier 
        "limit": 1500, # limit for the number of interaction partners
        "network_type": "physical", # type of interaction network
        "caller_identity" : "www.awesome_app.org" # identity of the caller
        }


    # Construct the request URL for getting interaction partners 

    request_url2 = "/".join([string_api_url, output_format, "interaction_partners"])


    # Make the API request to get interaction partners for the STRING IDs

    interactome_results = requests.post(request_url2, data=new_params)


    # Collect all unique gene names including the input proteins and their interactors

    all_genes = set(protein_genes)
    for line in interactome_results.text.strip().split("\n"):
        l = line.strip().split("\t")
        interactor = l[3]
        all_genes.add(interactor)


    # Convert the set of unique genes to a list and return it

    all_unqiue_genes = list(all_genes)
    return all_unqiue_genes





