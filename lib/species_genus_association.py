import pandas as pd

def association(df, levels, level):
    
    df_associated_species = df.loc[df['Species'] != pd.NA] #present rows in general dataframe
    df_associated_genus = df.loc[df['Genus'] != pd.NA] #present rows in general dataframe
    levels_copy = levels.copy()
    levels_copy.remove(level)
    levels_copy.append("Kingdom") #for dropping
    df_associated_species = df_associated_species.drop(columns=levels_copy).groupby(level).sum()
    df_associated_genus = df_associated_genus.drop(columns=levels_copy).groupby(level).sum()

    return  df_associated_species, df_associated_genus

