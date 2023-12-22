import os
from .data_analyzer import analyze_file
from .__var__ import *
from fuzzywuzzy import fuzz

# this is the minimum percentage allowed for similarity between pathology and the input
# you can increase it for more precision
TRESHOLD = 80

def is_tsv_file(file_path):
	_, file_extension = os.path.splitext(file_path)
	return file_extension.lower() == ".tsv"

def user_prompt(pathology_name):
	inp = input(f" * Did you mean {RED}{pathology_name}{RESET}? (yes/no) > ")
	return TRUE if inp == "yes" else CONTINUE

def correct_file(pathology_name, user_input):
	ratio = fuzz.ratio(user_input.lower(), pathology_name.lower())
	if (ratio == 100):
		return TRUE
	elif (ratio >= TRESHOLD):
		return user_prompt(pathology_name)
	else:
		if user_input.lower() in pathology_name.lower() and len(user_input) > 3:
			return user_prompt(pathology_name)
		else:
			return FALSE

def	match_file(user_input, database_folder):
	all_files = os.listdir(database_folder)
	matching_file = None

	for file_name in all_files:
		if (not is_tsv_file(file_name)):
			continue
		pathology_name, file_extension = os.path.splitext(file_name)
		res = correct_file(pathology_name.replace('_', ' '), user_input)
		if res == TRUE:
			display_pathology(pathology_name.replace('_', ' ').capitalize())
			matching_file = os.path.join(database_folder, pathology_name) + ".tsv"
			break
		elif res == CONTINUE:
			continue
		
	return (matching_file)

def display_pathology(input):
	print(f"\n╔═════════════════╗")
	print(f"║    Pathology    ║  {WHITE}{input}{RESET}")
	print(f"╚═════════════════╝\n")

def start_search(input, database_folder):
	_file = match_file(input, database_folder)
	if _file == None:
		print(f"{RED}Error : Please ensure that the pathology name is correct.")
		print(f"        * Work on this pathology is ongoing, please try again later.{RESET}")
		return
	analyze_file(_file)