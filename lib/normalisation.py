def normalise_and_cut(present_df, level):

    """ Function to normalise agora2 mapped reads data. Removes relative abundance <1e-5 and re-normalizes. Also
    transforms the genus/species name to its agora2 name
       Parameters
       ----------
       present_df : pandas dataframe, required
           dataframe with genera / species that are present in the agora2 database
       level : string, required
            species or genus, used as file name to distinguish files from each other
       Returns
       -------
       agora_normed_cut : pandas dataframe

       Authors
       -------
       Bram Nap - 04/2022
       """

    # Rename the genus/species name to its agora2 identifier
    present_df = present_df.set_index(present_df.index.str.replace(" ", "_", regex=True))
    present_df = present_df.set_index('pan' + present_df.index.astype(str))
    present_df = present_df.sort_index()

    # Calculate the relative abundances
    total_reads = present_df.sum()
    agora_normed = present_df.loc[:].div(total_reads)

    # Set all abundance lower than 1e-5 to 0
    agora_normed_cut = agora_normed.copy()
    agora_normed_cut[agora_normed_cut < 1e-5] = 0

    agora_normed_cut = agora_normed_cut.loc[(agora_normed_cut.sum(axis=1) != 0), (agora_normed_cut.sum(axis=0) != 0)]

    # Renormalize
    total_rel_abund = agora_normed_cut.sum()
    agora_renormed = agora_normed_cut.loc[:].div(total_rel_abund)

    # Fix some naming issues
    agora_renormed.index = agora_renormed.index.str.replace(".", "", regex=True)
    agora_renormed.index = agora_renormed.index.str.replace("-", "_", regex=True)
    agora_renormed.index = agora_renormed.index.str.replace("(", "_", regex=True)
    agora_renormed.index = agora_renormed.index.str.replace(")", "_", regex=True)
    agora_renormed.index = agora_renormed.index.str.replace("/", "_", regex=True)

    # Save
    agora_normed.to_csv(f'MARS_output/{level}/agora_normed_{level}.csv')
    agora_normed_cut.to_csv(f'MARS_output/{level}/agora_normed_cut_{level}.csv')
    agora_renormed.to_csv(f'MARS_output/{level}/agora_renormed_cut_{level}.csv')

    return agora_normed_cut, agora_renormed
