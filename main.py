from lib import preprocessing, agora_checking, pipeline, general_stats, stratification, normalisation
import os
import json

def main(*args, relative=False, path_to_stratification_file=None, **kwargs):

    """ 

    Parameters 
    ---------- 


    Returns 
    ------- 
    None 

    Authors
    -------
    Tim Hulshof
    Bram Nap
    """

    # Read in Agora2 unique taxa on all taxonomic levels
    with open('mars.json', 'r') as fp:
        agora2 = json.load(fp)

    agora2_phyla = set(agora2["Phylum"])
    agora2_classes = set(agora2["Class"])
    agora2_orders = set(agora2["Order"])
    agora2_families = set(agora2["Family"])
    agora2_genera = set(agora2["Genus"])
    agora2_species = set(agora2["Species"])
    agora2_strains = set(agora2["Strain"])

    # Total reads
    df, kingdom_df, phylum_df, class_df, order_df, family_df, genus_df, species_df, strain_df = preprocessing.preprocessing(**kwargs, relative=relative)

    # Retrieve present species and genus df for later use
    _, present_genus_df = agora_checking.agora_checking(genus_df, agora2_genera)
    _, present_species_df = agora_checking.agora_checking(species_df, agora2_species)

    levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Strain']

    # TODO: expand saving capabilities
    if not os.path.isdir("MARS_output"):
        os.mkdir("MARS_output")

    # Save requested files
    for arg in args:
        try:
            arg = arg.lower()
            if arg == "class":
                _, _ = pipeline.pipeline(df, class_df, levels, "Class", agora2_classes, agora2_species, agora2_genera)
            elif arg == "order":
                _, _ = pipeline.pipeline(df, order_df, levels, "Order", agora2_orders, agora2_species, agora2_genera)
            elif arg == "family":
                _, _ = pipeline.pipeline(df, family_df, levels, "Family", agora2_families, agora2_species, agora2_genera)
            elif arg == "strain":
                print("Strain pipeline still under construction...")
                pass
            else:
                print(f"\"{arg}\" did not match any of the optional taxonomic levels. Please check spelling")
        except SyntaxError:
            print("A syntax error was found in your arguments. Please check that you inputted a string.")
            pass
    
    # Phylum data for general stats
    species_phylum_list, genus_phylum_list = pipeline.pipeline(df, phylum_df, levels, "Phylum", agora2_phyla, agora2_species, agora2_genera)

    # agora_sepecies_normed - just saved?

    agora_species_normed, agora_species_normed_cut, agora_species_renormed = normalisation.normalise_and_cut(present_species_df)
    agora_genus_normed, agora_genus_normed_cut, agora_genus_renormed = normalisation.normalise_and_cut(present_genus_df)

    species_df_list = [present_species_df, species_df, agora_species_normed_cut, agora_species_renormed]
    genus_df_list = [present_genus_df, genus_df, agora_genus_normed_cut, agora_genus_renormed]

    # Get stats
    species_stats = general_stats.general_stats(df, species_phylum_list, species_df_list)
    genus_stats = general_stats.general_stats(df, genus_phylum_list, genus_df_list)

    if path_to_stratification_file is not None:
        species_strat = stratification.split_df(species_stats, path_to_stratification_file)
        genus_strat = stratification.split_df(genus_stats, path_to_stratification_file)

    # return present_genus_df, present_species_df


if __name__ == "__main__":

    genus, species = main(taxonomy_table=r"C:\Users\MSPG\Desktop\Mars_test\taxonomyWoL.tsv",
                          feature_table=r"C:\Users\MSPG\Desktop\Mars_test\feature-tableWoLgenome.txt",
                          path_to_stratification_file=r"C:\Users\MSPG\OneDrive - National University of Ireland, Galway (1)\MARS\Test_files\Strat_file.xlsx")