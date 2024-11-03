import os
import json
import urllib.request
import pandas as pd


def fetch_taxid_and_name(accession):
    try:
        # Search for the assembly ID
        search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&term={accession}&retmode=json"
        with urllib.request.urlopen(search_url) as response:
            search_response = json.loads(response.read().decode())
        assembly_id_list = search_response['esearchresult']['idlist']

        if not assembly_id_list:
            return "Not Found", "Not Found"

        assembly_id = assembly_id_list[0]

        # Link the assembly ID to the taxonomy database to get the TaxID
        link_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=assembly&db=taxonomy&id={assembly_id}&retmode=json"
        with urllib.request.urlopen(link_url) as response:
            link_response = json.loads(response.read().decode())
        tax_id = link_response['linksets'][0]['linksetdbs'][0]['links'][0]

        # Fetch the scientific name using the TaxID
        summary_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=taxonomy&id={tax_id}&retmode=json"
        with urllib.request.urlopen(summary_url) as response:
            summary_response = json.loads(response.read().decode())
        scientific_name = summary_response['result'][str(tax_id)]['scientificname']

        return tax_id, scientific_name
    except Exception as e:
        print(f"Error fetching TaxID and scientific name for {accession}: {e}")
        return "Error", "Error"


# Function to load the cache
def load_cache(cache_file: str = '/tmp/acc_to_taxid_cache.json'):
    if os.path.isfile(cache_file):
        with open(cache_file, 'r') as f:
            return json.load(f)
    return {}


# Function to save the cache
def save_cache(cache, cache_file: str = '/tmp/acc_to_taxid_cache.json'):
    with open(cache_file, 'w') as f:
        json.dump(cache, f)


def load_sample_sheet(
        sample_sheet: str,
        identifier_col: str = 'Identifier',
        expected_species_col: str = 'Species',
        expected_taxid_col: str = 'TaxID'
):
    data = pd.read_csv(sample_sheet, sep='\t')
    for col in [identifier_col, expected_species_col, expected_taxid_col]:
        assert col in data.columns, f"Column {col} not found in {data.columns}"
    # Create dictionary: {identifier: (species, taxid)}
    sample_sheet_dict = dict(zip(data[identifier_col], zip(data[expected_species_col], data[expected_taxid_col])))
    return sample_sheet_dict


def load_gtdbtk(
        gtdbtk_summary: str,
        sample_sheet: str = None,
        identifier_col: str = 'Identifier',
        expected_species_col: str = 'Species',
        expected_taxid_col: str = 'TaxID',
        out_file: str = None,
):
    if out_file is None:
        out_file = gtdbtk_summary.removesuffix('.tsv') + '.with-ncbi.tsv'

    if sample_sheet:
        sample_sheet_dict = load_sample_sheet(sample_sheet, identifier_col, expected_species_col, expected_taxid_col)

    data = pd.read_csv(gtdbtk_summary, sep='\t')

    # Extract the relevant accession numbers from your data
    accessions = data['closest_genome_reference'].unique()

    # Load the cache
    accession_taxid_name_map = load_cache()

    # Fetch missing accessions
    for acc in accessions:
        if acc not in accession_taxid_name_map:
            accession_taxid_name_map[acc] = fetch_taxid_and_name(acc)

    # Save the updated cache
    save_cache(accession_taxid_name_map)

    # Separate the TaxID and scientific name into two columns
    data['ncbi_taxid'] = data['closest_genome_reference'].map(lambda x: accession_taxid_name_map[x][0])
    data['scientific_name'] = data['closest_genome_reference'].map(lambda x: accession_taxid_name_map[x][1])

    if sample_sheet:
        data['expected_species'] = data['user_genome'].map(lambda x: sample_sheet_dict.get(x, ('-', '-'))[0])
        data['expected_taxid'] = data['user_genome'].map(lambda x: sample_sheet_dict.get(x, ('-', '-'))[1])

    # Save the updated DataFrame
    data.to_csv(out_file, sep='\t', index=False)
    print(f"Created {out_file}.")


if __name__ == '__main__':
    from fire import Fire

    Fire(load_gtdbtk)
