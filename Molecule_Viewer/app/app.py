# Import Standard Libraries
from rdkit_package.essentials import rdkit_essentials
from shiny import App, ui, render, reactive
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw
import tempfile

# Initialize RDKit Module
essential_methods = rdkit_essentials()

# Define UI Layout
app_ui = ui.page_fluid(
    ui.input_text("smiles", "Enter a SMILES string to add to list:"),
    ui.input_action_button("add_smiles", "Add SMILES"),
    ui.input_action_button("visualize_molecules", "Visualize Molecules"),
    ui.output_text("output_smiles_list"),
    ui.output_image("output_mol_structures"),
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
                # Validate and import SMILES string using RDKit
                try: 
                    validated_smiles = essential_methods.import_smiles(smiles)
                    smiles_list.set(validated_smiles)
                except ValueError as e:
                    # If invalid, show error message
                    session.show_modal(
                        ui.modal(
                            f'Error: {str(e)}',
                            title="Invalid SMILES",
                            easy_close=True,
                            footer=ui.modal_button("Close")
                        )
                    )

                # Reset the input field to accept new SMILES
                session.send_input_message("smiles", {"value": ""})

    @output
    @render.text
    def output_smiles_list():
        
        # Fetch SMILES list stored in reactive value
        lst = smiles_list.get()
        
        # If the list is empty, return a message
        if not lst:
            return "No SMILES strings added yet." 
        
        # Else, return joined list as string for display
        return "Current SMILES list: " + ", ".join(lst) 
    
    @output
    @render.image
    @reactive.event(input.visualize_molecules)
    def output_mol_structures():
        
        # Fetch SMILES list stored in reactive value
        lst = smiles_list.get()
        
        # If list is empty, return None
        if not lst:
            return None
        
        # Generate RDKit Mol objects from SMILES strings
        essential_methods.mol_from_smiles(lst)

        # Generate grid image of molecule structures
        img = essential_methods.visualize_molecules(lst)

        # If there is no image, return None
        if img is None:
            return None
        
        # Store image in temporary file and return its path for rendering
        # This is necessary because Shiny expects a file path for image output
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            img.save(tmp.name)
            tmp.flush()
            return {
                "src": tmp.name,
                "alt": "Molecule Grid",
                "width": img.width,
                "height": img.height,
                "type": "image/png"
            }

# Lauunch App
app = App(app_ui, server)