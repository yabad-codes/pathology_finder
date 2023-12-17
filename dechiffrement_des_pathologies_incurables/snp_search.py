from .__var__ import *
from Bio import Entrez
import xml.etree.ElementTree as ET
import os, ssl
from time import sleep
import textwrap

def get_ncbi_data(gene):
	Entrez.email = "yabadwashere@gmail.com"

	handle = Entrez.esearch(db="gene", term=f"{gene}[Gene Name] AND Homo sapiens[Organism]")
	record = Entrez.read(handle)
	handle.close()

	gene_ids = record["IdList"]

	if not gene_ids:
		print(f"The gene \'{gene}\' was not found in the database. Please check the spelling or try a different gene name.")
		return
	
	handle = Entrez.esummary(db="gene", id=",".join(gene_ids), rettype="gb", retmode="text")
	gene_records = handle.read()
	handle.close()
	return (gene_records)

def extract_summary(gene_info):
	root = ET.fromstring(gene_info)
	summary = root.find('.//Summary')
	return summary.text

def display_summary(gene, summary):
	print(f"\n╔════════════════╗")
	print(f"║    Function    ║  {RED}{gene}{RESET}")
	print(f"╚════════════════╝")
	print("╔═══════════════════════════════════════════════════════════════════════════════════╗")
	try:
		wrapped_text = textwrap.wrap(summary, 80)
		for text in wrapped_text:
			print(f"{WHITE}║ {text}{RESET}")
	except:
		print("║", "  Information for this gene is currently unavailable. Please try again later.")
	print("╚═══════════════════════════════════════════════════════════════════════════════════╝")

def process_input(num, dict):
	matched_genes = []
	for key, value in dict.items():
		if value[-1] == num:
			print(f"\n\n{RESET} * This SNP ({RED}{key}{RESET}) influences the function of the following gene(s):")
			matched_genes = value[:-1]
			for gene in matched_genes:
				print(f"{GREEN}    - {gene}")
			break
	while (True):
		user_input = input(f"\n{RESET} * Please enter the name of the gene you want to search for >{WHITE} ")
		if user_input.upper() in matched_genes:
			gene_info = get_ncbi_data(user_input.upper())
			summary = extract_summary(gene_info)
			display_summary(user_input.upper(), summary)
			sleep(2)
			more_input = input(f"\n{RESET} * Do you wish to continue searching for additional genes within the same {RED}SNP{RESET}? (yes/no) > {WHITE}")
			if (more_input == "yes"):
				print(f"{WHITE} Please note that the only genes available for search are {GREEN}{matched_genes}{RESET}.")
				continue
			else:
				break
		else:
			print(f"{RED} Error: The valid gene names to select from are {GREEN}{matched_genes}{RESET}")
			continue

def valid_number(num, dict):
	return True if num in range(1, len(dict) + 1) else False

def prompt_from_snp(dict):
	num = 0
	while True:
		user_input = input(f"{RESET}\n\t* Please enter the number corresponding to the {RED}SNP{RESET} you are interested in > {WHITE}")
		print(f"{RESET}", end="")
		try:
			num = int(user_input)
		except:
			print(f"{RED}\tError: Please enter a number.{RESET}")
			continue
		if not valid_number(num, dict):
			print(f"{RED}\tError: Please enter a number between {1, len(dict)} for this pathology.{RESET}")
			continue
		process_input(num, dict)
		break

def create_context():
	if (not os.environ.get('PYTHONHTTPSVERIFY', '') and
	 	getattr(ssl, '_create_unverified_context', None)):
		ssl._create_default_https_context = ssl._create_unverified_context

def get_snp_info(dict):
	create_context()
	print("\n")
	sleep(1)
	while True:
		user_input = input(f"{RESET} * Do you wish to obtain more information on a specific {RED}SNP{RESET} (yes/no) > {WHITE}")
		if user_input.lower() == "no":
			break
		elif user_input.lower() == "yes":
			prompt_from_snp(dict)
			continue
		else:
			continue