
def panname_2_phylum(total_df, panname_df, level):

    phyla_from_pan = panname_df.copy()
    pannames = []
    for name in phyla_from_pan.index:
        pannames += [name]

    alt_names = pannames.copy()
    alt_names = [sub.replace('_', ' ') for sub in alt_names]
    alt_names = [sub.replace('pan', '', 1) for sub in alt_names]

    total_df[level] = total_df[level].str.replace('_', ' ', regex=True)

    pan_phylum_dict = {}
    for i in range(len(alt_names)):
        row = total_df.loc[total_df[level] == alt_names[i]]
        phylum = row['Phylum'].values[0]
        pan_phylum_dict[pannames[i]] = phylum

    phyla_from_pan = phyla_from_pan.rename(indephyla_from_pan=pan_phylum_dict)
    phyla_from_pan = phyla_from_pan.groupby(level).sum()

    return phyla_from_pan
