from lib import preprocessing, agora_checking, pipeline, general_overview, stratification, normalisation, panname_2_phylum, visuals
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

    # TODO: sort args into taxonomic levels and phyla
    # TODO: Add sanity checker -> e.g., create a txt file with sample names that exhibit strange values and the user should take a closer look at.

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
    absent_genus_df, present_genus_df = agora_checking.agora_checking(genus_df, agora2_genera)
    absent_species_df, present_species_df = agora_checking.agora_checking(species_df, agora2_species)

    # absent make up
    def absent_make_up(df, total_df):
        total=total_df.sum()
        df["loss"] = df.loc[:].div(total).sum(axis=1).div(len(df.columns)-1)
        return df
        # print(total.sum())
        # print(df.sum(axis=1).div(total.sum()))

    absent_species_df = absent_make_up(absent_species_df, species_df)
    absent_genus_df = absent_make_up(absent_genus_df, genus_df)

    levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Strain']

    # arg sorting
    taxonomic_levels = []
    extra_phyla = []

    lowered_levels = list(map(lambda x: x.lower(), levels))

    for arg in args:
        if arg.lower() in lowered_levels:
            taxonomic_levels.append(arg.lower())
        else:
            extra_phyla.append(arg.lower())

    # TODO: expand saving capabilities
    if not os.path.isdir("MARS_output"):
        os.mkdir("MARS_output")
        os.mkdir('MARS_output/Figs')
        os.mkdir('MARS_output/Species')
        os.mkdir('MARS_output/Genus')
        os.mkdir('MARS_output/Phylum')

    absent_species_df.to_csv("MARS_output/Species/absent_species.csv")
    present_species_df.to_csv("MARS_output/Species/present_species.csv")
    absent_genus_df.to_csv("MARS_output/Genus/absent_genus.csv")
    present_genus_df.to_csv("MARS_output/Genus/present_genus.csv")

    # Save requested files
    for taxa in taxonomic_levels:
        try:
            if taxa == "class":
                _, _ = pipeline.pipeline(df, class_df, levels, "Class", agora2_classes, agora2_species, agora2_genera)
            elif taxa == "order":
                _, _ = pipeline.pipeline(df, order_df, levels, "Order", agora2_orders, agora2_species, agora2_genera)
            elif taxa == "family":
                _, _ = pipeline.pipeline(df, family_df, levels, "Family", agora2_families, agora2_species, agora2_genera)
            elif taxa == "strain":
                print("Strain pipeline still under construction...")
                pass
            else:
                print(f"\"{taxa}\" did not match any of the optional taxonomic levels. Please check spelling")
        except SyntaxError:
            print("A syntax error was found in your arguments. Please check that you inputted a string.")
            pass
    
    # Phylum data for general stats
    species_phylum_list, genus_phylum_list = pipeline.pipeline(df, phylum_df, levels, "Phylum", agora2_phyla, agora2_species, agora2_genera)

    # agora_sepecies_normed - just saved?

    [agora_species_normed_cut, agora_species_cut_renormed] = normalisation.normalise_and_cut(present_species_df, "Species")
    [agora_genus_normed_cut, agora_genus_cut_renormed] = normalisation.normalise_and_cut(present_genus_df, "Genus")

    pan_phylum_species = panname_2_phylum.panname_2_phylum(df, agora_species_cut_renormed, 'Species')
    pan_phylum_genus = panname_2_phylum.panname_2_phylum(df, agora_genus_cut_renormed, 'Genus')

    species_df_list = [present_species_df, species_df, agora_species_normed_cut, pan_phylum_species]
    genus_df_list = [present_genus_df, genus_df, agora_genus_normed_cut,pan_phylum_genus]

    # Get stats
    species_stats = general_overview.general_overview(df, species_phylum_list, species_df_list, "Species")
    genus_stats = general_overview.general_overview(df, genus_phylum_list, genus_df_list, "Genus")

    if path_to_stratification_file is not None:
        stratification.general_statistics_on_groups(species_stats, path_to_stratification_file, 'Species')
        stratification.general_statistics_on_groups(genus_stats, path_to_stratification_file, 'Genus')

    visuals.visualise(species_stats, 'Species')
    # return present_genus_df, present_species_df


# if __name__ == "__main__":

#     main(taxonomy_table=r"C:\Users\THuls\Documents\python_projects\Test\files\taxonomyWoL.tsv",
#                           feature_table=r"C:\Users\THuls\Documents\python_projects\Test\files\feature-table.txt")
