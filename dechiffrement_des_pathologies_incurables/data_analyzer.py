import csv
import bisect
from decimal import Decimal
from prettytable import PrettyTable
from .snp_search import *
from .__var__ import *

def get_data_file(csv_reader):
	data = []
	for row in csv_reader:
		data.append(row)
	return (data)

def get_min(data, informationLevel):
	sorted_values = []
	min_list = []

	for row in data:
		try:
			num = Decimal(row[PVALUE])
		except Exception:
			continue
		index = bisect.bisect_left(sorted_values, num)
		if index == len(sorted_values) or sorted_values[index] != num:
			bisect.insort_left(sorted_values, num)
	for i in range(0, informationLevel):
		min_list.append(sorted_values[i])
	return min_list

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

def table_edges(table):
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

def my_table_setup(table):
	table.field_names = [f"{WHITE}column 1{RESET}", f"{WHITE}column 2{RESET}", f"{WHITE}column 3{RESET}"]
	table.align = 'l'
	table.padding_width = 2
	table_edges(table)

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

	print(f"\n{WHITE}The risk SNPs are: {RESET}")
	for i in range(0, len(keys), 3):
		row = []
		for j in range(3):
			k_index = i + j
			if k_index < len(keys):
				key = keys[k_index]
				value = dict[key]
				if value[-1] <= 5:
					color_key = f"{WHITE}({index}){RESET}  {GREEN}{key}{RESET}"
				else:
					color_key = f"{WHITE}({index}){RESET}  {RED}{key}{RESET}"
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
		genes, thres = genes_info[:-1], genes_info[-1]
		color = GREEN if thres <= 5 else RED
		colored_genes = [f"{color}{gene}{RESET}" for gene in genes]
		rows.extend(colored_genes)
	
	for i in range(0, len(rows), 3):
		row = rows[i:i+3]
		while len(row) < 3:
			row.append('')
		table.add_row(row)
	my_table_setup(table)
	print(table)

def struct_dict(data, min_list):
	if (not len(min_list)):
		print(f"{RED}Error :")
		print(f"        We were able to find the file linked to the pathology, but it seems that the file is corrupted.{RESET}")
		return
	
	allele_to_genes = {}
	reference = []
	for row in data:
		try:
			pvalue = Decimal(row[PVALUE])
		except:
			continue
		if pvalue in min_list:
			try:
				index = min_list.index(pvalue) + 1
			except:
				continue
			allele = row[RISKALLELE]
			genes = row[GENES].split(',')
			genes.append(index)
			allele_to_genes[allele] = genes
			reference.append(row[ACCESSIONID])
	sorted_allele_to_genes = dict(sorted(allele_to_genes.items(), key=lambda item: (item[1][-1], item[0])))
	display_snps(sorted_allele_to_genes)
	display_genes(sorted_allele_to_genes)
	get_snp_info(get_enum_dict(sorted_allele_to_genes))
	display_reference(reference)

def draw_box(message, field):
	table = PrettyTable()
	
	table.field_names = [field]
	table.add_row([message])
	table.align[field] = "l"
	table._max_width = {field: 50}
	table_edges(table)
	print(table)

def get_information_level():
	message = f"{WHITE}The threshold of significance is a scale from 1 to 10 that determines the number of SNPs (Single Nucleotide Polymorphisms) displayed. A setting of 1 shows the most significant SNP, while 10 displays the least significant ones.{RESET}"
	field = f"{WHITE}Threshold Of Significance{RESET}"
	draw_box(message, field)
	user_input = input(f"\n* Please enter threshold of significance {RED}(1 - 10){RESET} > {WHITE}")
	try:
		thsi = int(user_input)
		if thsi >= 1 and thsi <= 10:
			return thsi
		raise "Number has to be valid"
	except:
		print(f"{RED}Error : The number has to be a digit between (1 - 10), please try again.")
		exit(1)

def analyze_file(_file):
	informationLevel = get_information_level()
	with open(_file, 'r', encoding='utf-8') as file:
		csv_reader = csv.reader(file, delimiter='\t')
		next(csv_reader, None)
		data = get_data_file(csv_reader)
		min_list = get_min(data, informationLevel)
		struct_dict(data, min_list)
