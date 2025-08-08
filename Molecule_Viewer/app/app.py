# Import Standard Libraries
from rdkit_package.essentials import rdkit_essentials
from shiny import App, ui, render, reactive
import tempfile

# Define UI Layout
app_ui = ui.page_fluid(
    ui.tags.script("Shiny.setInputValue('page_loaded', True);"),
    ui.tags.title("Molecule Viewer"),
    ui.tags.div(
        ui.tags.h1(
            "Molecule Viewer",
            style="font-size: 36px; margin-bottom: 20px; margin-top: 20px;"),
        ui.tags.p(
            "This app allows you to input SMILES strings, visualize molecules, and calculate their similarities."
            ),
        ui.tags.p(
            "Enter a SMILES string to add it to the list."
            ),
        ui.tags.p(
            "You can visualize and calculate Tanimoto similarity between the molecules once they are added."
            ),
        style="text-align:center; font-size: 12px"
    ),
    ui.tags.hr(),
    ui.navset_card_tab(
        ui.nav_panel(
            "Import",
            ui.tags.div(
                ui.input_text("smiles", "Enter a SMILES string to add to list:"),
                ui.input_action_button("add_smiles", "Add SMILES"),
                ui.tags.br(),
                ui.output_text("output_smiles_list"),
                style="display: flex; flex-direction: column; align-items: center; text-align: center; margin-bottom: 20px;"
                ),
            ),
        ui.nav_panel(
            "Visualize",
            ui.tags.div(
                ui.output_image("output_mol_structures"),
                style="display: flex; flex-direction: column; align-items: center; text-align: center; margin-bottom: 20px;"
                ),
            ),
        ui.nav_panel(
            "Calculate Similarity",
            ui.tags.div(
                ui.tags.div(
                    ui.input_select("cpd1", "Select first molecule:", choices={}),
                    ui.input_select("cpd2", "Select second molecule:", choices={}),
                    ),
                ui.input_action_button('calculate_similarity', 'Calculate Similarity'),
                ui.tags.br(),
                ui.output_text("output_similarity"),
                style="display: flex; flex-direction: column; align-items: center; text-align: center; margin-bottom: 20px;"
                ),
            ),
        ui.nav_panel("Substructure Search", 
                    ui.tags.div(    
                        ui.tags.div(
                            ui.input_select("mol", "Select the molecule you want to search in:", choices={}),
                            ui.input_text("substructure", "Enter a substructure SMILES string to search for:"),
                        ),
                        ui.input_action_button("search_substructure", "Search Substructure"),
                        ui.tags.br(),
                        ui.output_image("output_substructure_search"),
                        style="display: flex; flex-direction: column; align-items: center; text-align: center; margin-bottom: 20px;"
                    )
                ),
        id="nav_tab"
        ),
    )

# Define Server Logic
def server(input, output, session):  
    
    # Initialize RDKit Module
    essential_methods = rdkit_essentials()

    # Initialize a reactive value to store list of SMILES strings
    smiles_list = reactive.Value([])

    # Initialize a reactive value to store Tanimoto similarity results
    similarity_result = reactive.Value("")

    @reactive.effect
    @reactive.event(input.page_loaded)
    def _():

        # Reset stored values when the page is loaded
        if input.page_loaded():
            essential_methods.smiles = []
            essential_methods.ms = []
            essential_methods.fps = []

    @reactive.effect
    @reactive.event(input.add_smiles)
    def _():
        if input.add_smiles():
            # Fetch smiles string from input
            smiles = input.smiles()

            if smiles:
                # Validate and import SMILES string using RDKit and update the list
                try:
                    essential_methods.import_smiles(smiles)
                    updated_smiles_list = smiles_list.get() + [smiles]
                    smiles_list.set(updated_smiles_list)

                # If invalid, show error message
                except ValueError as e:
                    ui.modal_show(
                        ui.modal(
                            f'Error: {str(e)}',
                            title="Invalid SMILES",
                            easy_close=True,
                            footer=ui.modal_button("Close")
                        )
                    )

                # Reset the input field to accept new SMILES
                session.send_input_message("smiles", {"value": ""})

    @reactive.effect
    @reactive.event(input.nav_tab)
    def _():
        if input.nav_tab() == "Calculate Similarity":
            ui.update_select('cpd1', choices={k:v for k, v in enumerate(smiles_list.get())})
            ui.update_select('cpd2', choices={k:v for k, v in enumerate(smiles_list.get())})

    @reactive.effect
    @reactive.event(input.nav_tab)
    def _():
        if input.nav_tab() == "Substructure Search":
            ui.update_select('mol', choices={k:v for k, v in enumerate(smiles_list.get())})

    @reactive.effect
    @reactive.event(input.calculate_similarity)
    def _():
        if input.calculate_similarity():
            if not input.cpd1() != input.cpd2():
                ui.modal_show(
                    ui.modal(
                        'Please select two different molecules to calculate similarity.',
                        title="Selection Error",
                        easy_close=True,
                        footer=ui.modal_button("Close")
                    ),
                )
            try:
                essential_methods.generate_fingerprints()
                similarity = essential_methods.calculate_similarity(
                    int(input.cpd1()),
                    int(input.cpd2()),
                )
                similarity_result.set(f"Tanimoto Similarity: {similarity:.4f}")
            except Exception as e:
                ui.modal_show(
                    ui.modal(
                        f'Error calculating similarity: {str(e)}',
                        title="Calculation Error",
                        easy_close=True,
                        footer=ui.modal_button("Close")
                    ),
                )
            

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
    @reactive.event(input.nav_tab)
    def output_mol_structures():

        if not input.nav_tab() == "Visualize":
            return None
        
        # If no molecules are available, show a message
        if not len(essential_methods.ms) > 0:
            ui.modal_show(
                ui.modal(
                    f'Import some SMILES strings to visualize.',
                    title="No Molecules to Visualize",
                    easy_close=True,
                    footer=ui.modal_button("Close")
                )
            )
            return None

        # Generate grid image of molecule structures
        img = essential_methods.visualize_molecules()

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
        
    @output
    @render.text
    def output_similarity():
        result = similarity_result.get()
        if result:
            return result
        
    @output
    @render.image
    @reactive.event(input.search_substructure)
    def output_substructure_search():

        # Generate substructure search result image
        if not input.mol():
            ui.modal_show(
                ui.modal(
                    'Please select a molecule to search.',
                    title="Input Error",
                    easy_close=True,
                    footer=ui.modal_button("Close")
                )
            )
        if not input.substructure():
            ui.modal_show(
                ui.modal(
                    'Please enter a substructure SMILES string to search for.',
                    title="Input Error",
                    easy_close=True,
                    footer=ui.modal_button("Close")
                )
            )

        try:   
            img = essential_methods.substructure_search(int(input.mol()), input.substructure())

            # Store image in temporary file and return its path for rendering
            # This is necessary because Shiny expects a file path for image output
            with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
                img.save(tmp.name)
                tmp.flush()
                return {
                    "src": tmp.name,
                    "alt": "Substructure Search Result",
                    "width": img.width,
                    "height": img.height,
                    "type": "image/png"
                }
        except ValueError as e:
            ui.modal_show(
                ui.modal(
                    f'Error Generating Image: {str(e)}',
                    title="Substructure Search Error",
                    easy_close=True,
                    footer=ui.modal_button("Close")
                )
            )
            return None

# Lauunch App
app = App(app_ui, server)