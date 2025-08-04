# Import Standard Libraries
from rdkit_package.essentials import rdkit_essentials
from shiny import App, ui, render, reactive

# Initialize RDKit Module
essential_methods = rdkit_essentials()

# Define UI
app_ui = ui.page_fluid(
    ui.input_text("smiles", "Enter a SMILES string to add to list:"),
    ui.input_action_button("add_smiles", "Add SMILES"),
    ui.output_text("output_smiles_list"),
)

# Define Server Logic
def server(input, output, session):
    
    # Initialize a reactive value to store list of SMILES strings
    smiles_list = reactive.Value([])

    @reactive.effect
    @reactive.event(input.add_smiles)
    def _():
        
        if input.add_smiles():
            
            # Fetch smiles string from input
            smiles = input.smiles()
            
            if smiles:

                # Add the new SMILES to the list
                updated_list = smiles_list.get() + [smiles]
                
                # Update reactive value to ensure output reflects the change
                smiles_list.set(updated_list)

                # Reset the input field to accept new SMILES
                session.send_input_message("smiles", {"value": ""})

    @output
    @render.text
    def output_smiles_list():
        lst = smiles_list.get()
        if not lst:
            return "No SMILES strings added yet."
        return "Current SMILES list: " + ", ".join(lst)
    
    @output
    @render.image
    def output_cpd__structures():
        lst = smiles_list.get()
        if not lst:
            return None
        
        # Generate RDKit Mol objects from SMILES strings
        mols = essential_methods.mol_from_smiles(lst)
        
        # Visualize the first molecule as an image
        if mols:
            return essential_methods.visualize_molecules(0)
        else:
            return None


# Lauunch App
app = App(app_ui, server)