
def panname_2_phylum(total_df, panname_df, level):

    """ Takes a dataframe with pan names and converts it its corresponding phylum

        Parameters
        ----------
        total_df : pandas dataframe, required
            dataframe with all taxonomic information
        panname_df : pandas dataframe, required
            dataframe which contains the pannames as index
        level : string, required
            string stating Species or Genus

        Returns
        -------
        phyla_from_pan : pandas dataframe
            dataframe with pannames converted to phyla and grouped by phylum name

        Authors
        -------
        Bram Nap - 04/2022
        """

    # Initialise df
    phyla_from_pan = panname_df.copy()

    # Store all the pannames in a list
    pannames = []
    for name in phyla_from_pan.index:
        pannames += [name]

    # Replace the _ and the first pan instance in all the names
    alt_names = pannames.copy()
    alt_names = [sub.replace('_', ' ') for sub in alt_names]
    alt_names = [sub.replace('pan', '', 1) for sub in alt_names]

    # Replace all the _ in the reference df to ensure matches
    total_df[level] = total_df[level].str.replace('_', ' ', regex=True)

    # Create a dictionary with panname:phylum
    pan_phylum_dict = {}
    for i in range(len(alt_names)):
        row = total_df.loc[total_df[level] == alt_names[i]]
        phylum = row['Phylum'].values[0]
        pan_phylum_dict[pannames[i]] = phylum

    # Change all pannames to phyla and group and sum the same phyla
    phyla_from_pan = phyla_from_pan.rename(index=pan_phylum_dict)
    phyla_from_pan = phyla_from_pan.groupby(level).sum()

    return phyla_from_pan
