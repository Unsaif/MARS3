import pandas as pd


def panname_2_phylum(total_df, agora2_finished_df, level):

    x = agora2_finished_df.copy()
    pannames = []
    for i in x.index:
        pannames += [i]

    alt_names = pannames.copy()
    alt_names = [sub.replace('_', ' ') for sub in alt_names]
    alt_names = [sub.replace('pan', '', 1) for sub in alt_names]

    # print(df.loc[df['Species'] == 'Clostridiales bacterium 1_7_47FAA'])

    total_df[level] = total_df[level].str.replace('_', ' ', regex=True)

    pan_phylum_dict = {}
    for i in range(len(alt_names)):
        row = total_df.loc[total_df[level] == alt_names[i]]
        phylum = row['Phylum'].values[0]
        pan_phylum_dict[pannames[i]] = phylum

    x = x.rename(index = pan_phylum_dict)
    x = x.groupby(level).sum()

    return x
