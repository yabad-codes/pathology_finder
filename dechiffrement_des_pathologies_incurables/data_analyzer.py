import csv
from decimal import Decimal
from prettytable import PrettyTable
from .snp_search import *
from .__var__ import *

def get_data_file(csv_reader):
	data = []
	for row in csv_reader:
		data.append(row)
	return (data)

def get_min(data):
	min = Decimal(2)
	next_min = Decimal(2)

	for row in data:
		try:
			num = Decimal(row[PVALUE])
		except Exception:
			continue
		if num < min:
			next_min = min
			min = num
		elif num < next_min and num != min:
			next_min = num
		
	return min, next_min

def display_row(list, isMax, is_snp):
	indent = '\t\t' if is_snp else '    '
	color = YELLOW if not isMax else GREEN
	print(f"{color}", end="")
	for i in range(0, len(list), 3):
		row_elements = list[i:i + 3]
		print(indent.join(row_elements))
	print(f"{RESET}")

def display_reference(reference):
	print(f"\n╔═════════════════╗")
	print(f"║    Reference    ║  ")
	print(f"╚═════════════════╝\n")
	print(f"{WHITE}References of primary molecular data: {RESET}\n")
	unique_list = list(set(reference))
	for i in range(0, len(unique_list), 5):
		row_elements = unique_list[i:i+5]
		print('\t'.join(row_elements))
	print(f"\n{RED}NB: the references are presented by their \"Study accession\" codes in \"The NHGRI_EBI GWAS Catalog\".\n{WHITE}")

def my_table_setup(table):
	table.field_names = [f"{WHITE}column 1{RESET}", f"{WHITE}column 2{RESET}", f"{WHITE}column 3{RESET}"]
	table.align = 'l'
	table.padding_width = 2
	table.junction_char = '╬'
	table.top_junction_char = '╦'
	table.bottom_junction_char = '╩'
	table.left_junction_char = '╠'
	table.right_junction_char = '╣'
	table.top_left_junction_char = '╔'
	table.top_right_junction_char = '╗'
	table.bottom_left_junction_char = '╚'
	table.bottom_right_junction_char = '╝'
	table.horizontal_char = '═'
	table.vertical_char = '║'

def get_enum_dict(sorted_dict):
	indexed_dict = {}
	index = 1
	for key, value_list in sorted_dict.items():
		new_value_list = value_list[:-1] + [index]
		indexed_dict[key] = new_value_list
		index += 1
	return indexed_dict

def display_snps(dict):
	keys = list(dict.keys())
	table = PrettyTable()
	index = 1

	print(f"{WHITE}The risk SNPs are: {RESET}")
	for i in range(0, len(keys), 3):
		row = []
		for j in range(3):
			k_index = i + j
			if k_index < len(keys):
				key = keys[k_index]
				value = dict[key]
				if value[-1]:
					color_key = f"{WHITE}({index}){RESET}  {GREEN}{key}{RESET}"
				else:
					color_key = f"{WHITE}({index}){RESET}  {YELLOW}{key}{RESET}"
				row.append(color_key)
				index += 1
			else:
				row.append('')
		table.add_row(row)
	my_table_setup(table)
	print(table)


def display_genes(dict):
	table = PrettyTable()

	print(f"\n{WHITE}The genes involved in cell signaling are: {RESET}")
	rows = []
	for genes_info in dict.values():
		genes, isMax = genes_info[:-1], genes_info[-1]
		color = GREEN if isMax else YELLOW
		colored_genes = [f"{color}{gene}{RESET}" for gene in genes]
		rows.extend(colored_genes)
	
	for i in range(0, len(rows), 3):
		row = rows[i:i+3]
		while len(row) < 3:
			row.append('')
		table.add_row(row)
	my_table_setup(table)
	print(table)

def struct_dict(data, min, next_min):
	if (min == -1):
		print(f"{RED}Erreur :")
		print(f"        Nous avons pu trouver le fichier lié à la pathologie, mais il semble que le fichier soit corrompu.{RESET}")
		return
	
	allele_to_genes = {}
	reference = []
	for row in data:
		pvalue = Decimal(row[PVALUE])
		if min == pvalue or next_min == pvalue:
			isMin = True if (min == pvalue) else False
			allele = row[RISKALLELE]
			genes = row[GENES].split(',')
			genes.append(isMin)
			allele_to_genes[allele] = genes
			reference.append(row[ACCESSIONID])
	sorted_allele_to_genes = dict(sorted(allele_to_genes.items(), key=lambda item: (item[1][-1] == False, item[0])))
	display_snps(sorted_allele_to_genes)
	display_genes(sorted_allele_to_genes)
	get_snp_info(get_enum_dict(sorted_allele_to_genes))
	display_reference(reference)

def analyze_file(_file):
	with open(_file, 'r', encoding='utf-8') as file:
		csv_reader = csv.reader(file, delimiter='\t')
		next(csv_reader, None)
		data = get_data_file(csv_reader)
		[min, next_min] = get_min(data)
		struct_dict(data, min, next_min)
