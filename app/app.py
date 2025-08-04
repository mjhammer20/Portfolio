# Import Standard Libraries
import sys
import os
from shiny.express import input, ui, render

# Add src to system path to import custom modules
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

# Import Custom Modules 
from src.rdkit_module import rdkit_functions

# Initialize RDKit Module
rdkit_module = rdkit_functions()

# Define UI
app_ui = ui.page_fluid(
    ui.input_text("smiles", "Enter a SMILES string to add to list:"),
    ui.input_action_button("add_smiles", "Add SMILES"),
)

# Define Server Logic
def server(input, output, session):
    # Initialize a list to store SMILES strings
    smiles_list = []

    @input.add_smiles
    def add_smiles():
        # Get the SMILES string from input
        smiles = input.smiles()
        if smiles:
            smiles_list.append(smiles)
            rdkit_module.mol_from_smiles(smiles_list)
            print(f"Added SMILES: {smiles}")