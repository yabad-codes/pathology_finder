from dechiffrement_des_pathologies_incurables import start_search
from dechiffrement_des_pathologies_incurables.__var__ import GREEN, RESET, WHITE

def display_header():
	print(f"{GREEN}╔══════════════════════════════════════╗")
	print(f"║          Welcome to LIB FSAC         ║")
	print(f"╚══════════════════════════════════════╝{RESET}")

def get_input():
	display_header()
	user_input = input(f"\n * Please enter a pathology name > {WHITE}")
	print(f"{RESET}")
	return user_input

file = start_search(get_input(), "./database")
