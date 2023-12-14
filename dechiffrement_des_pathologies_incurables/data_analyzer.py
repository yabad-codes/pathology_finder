import csv
from decimal import Decimal
from .__var__ import *

def get_data_file(csv_reader):
	data = []
	for row in csv_reader:
		data.append(row)
	return (data)

def get_max(data):
	max = Decimal(-1)
	next_max = Decimal(-1)

	for row in data:
		try:
			num = Decimal(row[PVALUE])
		except Exception:
			continue
		if num > max:
			next_max = max
			max = num
		elif num > next_max and num != max:
			next_max = num
		
	return max, next_max

def display_row(list, isMax, is_snp):
	indent = '\t\t' if is_snp else '    '
	color = YELLOW if not isMax else GREEN
	print(f"{color}", end="")
	for i in range(0, len(list), 3):
		row_elements = list[i:i + 3]
		print(indent.join(row_elements))
	print(f"{RESET}")


def display_data(data, message, isSnp):
	print(f"{WHITE}{message}{RESET}")
	max_list = []
	next_max_list = []
	for element in data:
		if element[1] == True:
			max_list.append(element[0])
		else:
			next_max_list.append(element[0])
	display_row(max_list, True, isSnp)
	display_row(next_max_list, False, isSnp)

def display_reference(reference):
	print(f"\n╔═════════════════╗")
	print(f"║    Référence    ║  ")
	print(f"╚═════════════════╝\n")
	print(f"{WHITE}Références des données moléculaires primaires : {RESET}\n")
	unique_list = list(set(reference))
	for i in range(0, len(unique_list), 5):
		row_elements = unique_list[i:i+5]
		print('\t'.join(row_elements))
	print(f"\n{RED}NB: les références sont présentées par leurs codes \"Study accession\" sur \"The NHGRI_EBI GWAS Catalog\".\n{WHITE}")
	

def print_cols_formatted(data, max, next_max):
	risk_allele = []
	genes = []
	reference = []
	for row in data:
		if (max == -1):
			print(f"{RED}Erreur :")
			print(f"        Nous avons pu trouver le fichier lié à la pathologie, mais il semble que le fichier soit corrompu.{RESET}")
			return
		if max == Decimal(row[PVALUE]) or next_max == Decimal(row[PVALUE]):
			isMax = False
			if (max == Decimal(row[PVALUE])):
				isMax = True
			risk_allele.append((row[RISKALLELE], isMax))
			genes.extend([(f"{gene:<20}", isMax) for gene in row[GENES].split(',')])
			reference.append(row[ACCESSIONID])
	display_data(risk_allele, "Les SNPs de risque sont :", True)
	display_data(genes, "Les gènes impliqués dans la signalisation cellulaire de cette pathologie sont :", False)
	display_reference(reference)

def analyze_file(_file):
	with open(_file, 'r', encoding='utf-8') as file:
		csv_reader = csv.reader(file, delimiter='\t')
		next(csv_reader, None)
		data = get_data_file(csv_reader)
		[max, next_max] = get_max(data)
		print_cols_formatted(data, max, next_max)
