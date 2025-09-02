# Xavier Castellanos-Girouard
# Date first created: June 28 2023
# Date last modified: April 23 2024

## Import libraries
from bs4 import BeautifulSoup
import requests
import re
import pandas as pd


## This function extracts uniprot IDs from PDB page
def get_UniprotID(PDB_id):
	global PDB2Uniprot_list
	
	PPI_UniprotIDs = []
	PPI_UniprotIDs.append(PDB_id)
	
	html_text = requests.get(str('https://www.rcsb.org/structure/' + PDB_id)).text
	soup = BeautifulSoup(html_text, 'lxml')
	
	
	table_macromol1 = soup.find_all('table', id = "table_macromolecule-protein-entityId-1") # Find all instances of row cells in a table
	
	for table_row in table_macromol1: # For every macromolecule table
		links = table_row.find_all('a', class_ = "querySearchLink") # Find <a tag with class
		for link in links: # Find rows with UniprotID regex
			UniprotID_match = re.search("^[A-Z][0-9][A-Z0-9][A-Z0-9][A-Z0-9][0-9]$", link.text)
			if UniprotID_match != None: # If rows match UniprotID regex
				PPI_UniprotIDs.append(UniprotID_match.group(0)) # Append to ppi list
	
	
	
	table_macromol2 = soup.find_all('table', id = "table_macromolecule-protein-entityId-2") # Find all instances of row cells in a table
	
	for table_row in table_macromol2: # For every macromolecule table
		links = table_row.find_all('a', class_ = "querySearchLink") # Find <a tag with class
		for link in links: # Find rows with UniprotID regex
			UniprotID_match = re.search("^[A-Z][0-9][A-Z0-9][A-Z0-9][A-Z0-9][0-9]$", link.text)
			if UniprotID_match != None: # If rows match UniprotID regex
				PPI_UniprotIDs.append(UniprotID_match.group(0))
	
	PDB2Uniprot_list.append(PPI_UniprotIDs)


def import_PDB_id_list():
	
	# Dataframe containing PDBs and other info
	PDB_PPI_DF = pd.read_csv("../results/IntermediateFiles/Kd_PDB_list.csv", index_col = 0)
	
	PDB_id_list = list(PDB_PPI_DF["PDB_ID"].values)
	
	return(PDB_id_list)
	
	
	
def main():
	
	global PDB2Uniprot_list
	
	# Import list of PDB file names 
	PDB_id_list = import_PDB_id_list()
	for PDB_id in PDB_id_list: # For PDB file name in lab
		get_UniprotID(PDB_id) # Get Uniprot IDs
	
	#print(PDB2Uniprot_list)
	
	# Only retain instances where both Uniprot IDs were found (3 entries in output)
	PDB2Uniprot_list = [x for x in PDB2Uniprot_list.copy() if len(x) == 3]
	
	# Construct PDB to Uniprot conversion table
	PDB2Uniprot_DF = pd.DataFrame(PDB2Uniprot_list, columns = ['PDB_id', 'protein1_UniprotID', 'protein2_UniprotID'])
	
	## Output table
	PDB2Uniprot_DF.to_csv("../results/IntermediateFiles/PDB2Uniprot.csv")

PDB2Uniprot_list = []
main()

