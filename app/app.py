# Import Standard Libraries
from rdkit_package.essentials import rdkit_essentials
from shiny.express import input, ui, render

# Initialize RDKit Module
essential_methods = rdkit_essentials()

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
            essential_methods.mol_from_smiles(smiles_list)
            print(f"Added SMILES: {smiles}")