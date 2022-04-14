import pandas as pd
import requests

def preprocessing(relative = False, **kwargs):
    
    """ Takes a taxonomy table and feature table or an already combined table to 
    produce an abundance table for every taxonomic level. 

    Parameters 
    ---------- 
    taxonomy_table : str, optional
        Where taxonomy table is located. 
    feature_table : str, optional
        Where feature table is located. 
    combined : str, optional
        Where combined table is located. 
    relative : Bool, optional 
        Signifies if the feature table or combined table has already gone through the relative abundace process.

    Returns 
    ------- 
    kingdom_df, phylum_df, class_df, order_df, family_df, genus_df, species_df, strain_df
        The function returns an abundance dataframe for every taxonomy level from kingdom to strain. 

    Authors
    -------
    Tim Hulshof
    """
    # CAG
    alterations = [
        'Candidatus_', 
        'Candidatus ',
        'Candidatus', 
        '_[0-9]{1}$', 
        '_family', 
        '_order', 
        '_sensu_stricto', 
        ' [0-9]{1}$', 
        ' family', 
        ' order', 
        ' sensu stricto',
        '\[', '\]', 
        '\'',
        ' CAG.*$'
    ]

    specific_alterations = {
        #'_': '__',
        #'\. ': '_',
        #' ': '_',
        '\|': ';',
        '_Family_XI_Incertae_Sedis': ' Incertae Sedis XI',
        '_Family_XIII_Incertae_Sedis': ' Incertae Sedis XIII',
        ' Family XI Incertae Sedis': ' Incertae Sedis XI',
        ' Family XIII Incertae Sedis': ' Incertae Sedis XIII',
        'typhimurium': 'enterica',
        'Ruminococcus_gauvreauii_group': 'Ruminococcus',
        'Ruminococcus_gnavus_group': 'Blautia',
        'Ruminococcus_torques_group': 'Blautia',
        'Eubacterium_coprostanoligenes_group': 'Eubacterium',
        'Eubacterium_hallii_group': 'Anaerobutyricum',
        'Eubacterium_ruminantium_group': 'Eubacterium',
        'Eubacterium_ventriosum_group': 'Eubacterium',
        'Eubacterium_xylanophilum_group': 'Eubacterium',
        'Clostridium_innocuum_group': 'Erysipelatoclostridium',
        'Ruminococcus gauvreauii group': 'Ruminococcus',
        'Ruminococcus gnavus group': 'Blautia',
        'Ruminococcus torques group': 'Blautia',
        'Eubacterium coprostanoligenes group': 'Eubacterium',
        'Eubacterium hallii group': 'Anaerobutyricum',
        'Eubacterium ruminantium group': 'Eubacterium',
        'Eubacterium ventriosum group': 'Eubacterium',
        'Eubacterium xylanophilum group': 'Eubacterium',
        'Clostridium innocuum group': 'Erysipelatoclostridium', 
        ' group': '',
        '_group': ''
    }
    
    levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Strain'] 
    
    #homosynonym retrieval from api
    address = "https://marsagora2api.herokuapp.com/agora2/"
    response_homosynonyms = requests.get(address + "homosynonyms")
    homosynonyms = response_homosynonyms.json()
    
    def pad_or_truncate(some_list, target_len):
        return some_list[:target_len] + [""]*(target_len - len(some_list))
    
    def get_grouped_tax_level(level, df):
        copy = levels.copy() #shallow copy of taxonomic levels
        if combined and relative:
            copy.reverse()
            pos = copy.index(level) #get postional int for finding what rows correspond with taxa level
            copy.remove(level)
            copy = copy + ["occurences"]
            df = df[df["occurences"] == pos].drop(columns=copy)
            df = df.set_index(level)
        else:
            copy.remove(level)
            df = df.drop(columns = copy).groupby(level).sum()
        return df
    
    def separator(fname):
        if fname.endswith(".tsv") or fname.endswith(".txt"):
            return "\t"
        elif fname.endswith(".csv"):
            return ","
        elif fname.endswith(".xlsx"):
            return "excel"
        else:
            raise ValueError("Sorry, only file types .csv, .tsv, .txt or .xlsx are acceptable.")
            
    if len(kwargs) == 2: #taxonomy table and feature table
        combined = False
        try:
            taxonomy_table = kwargs["taxonomy_table"]
            feature_table = kwargs["feature_table"]
        except KeyError as err:
            raise KeyError("Please ensure variables are named \"taxonomy_table\" and \"feature_table\".")
            
    else:
        combined = True
        try:
            combined_table = kwargs["combined"]
        except KeyError as err:
            raise KeyError("Please ensure variable is named \"combined\".")
                   
    if combined:
        taxonomy_table = combined_table
    
    ##Taxonomy Read in | Combined read in
    try:
        tax_sep = separator(taxonomy_table)
        
        if tax_sep == "excel":
            tax = pd.read_excel(taxonomy_table)
        else:
            tax = pd.read_csv(taxonomy_table, sep=tax_sep)
            
        for key in specific_alterations:
            tax = tax.replace(key, specific_alterations[key], regex=True)
        for alteration in alterations:
            tax = tax.replace(alteration, "", regex=True)         
    except FileNotFoundError as err:
        raise FileNotFoundError("Taxonomy table not found. Please check file path.")
    
    if not combined:
        ##Feature Read in 
        try:
            feat_sep = separator(feature_table)
            
            if feat_sep == "excel":
                feat = pd.read_excel(feature_table)
            else:
                feat = pd.read_csv(feature_table, sep=feat_sep, low_memory=False)
                
            if len(feat.columns) == 1: #Deals with potential "# (header)" at top of file
                if feat_sep == "excel":
                    feat = pd.read_excel(feature_table, header = 1)
                else:
                    feat = pd.read_csv(feature_table, sep=feat_sep, header=1) 
                    
            feat = feat.set_index(feat.columns[0])
        except FileNotFoundError as err:
            raise FileNotFoundError("Feature table not found. Please check file path.")
            
    if combined and relative:
        df = tax.iloc[:, 0].replace(".__", "", regex=True).str.split(';', expand=False)
        df = df.apply(lambda x: pad_or_truncate(x, 8)) #ensure every taxonomic level has empty spaces for levels that are not represented
        df = pd.DataFrame(df.to_list(), columns=levels)
        df = df.merge(tax.iloc[: , 1:], left_index=True, right_index=True, how='inner')
        df = df.replace("", pd.NA)
        df = df.replace("_", " ", regex=True)
        df["occurences"] = df.apply(lambda x: x.isna().sum(), axis=1) #Occurences of <NA>
    elif combined and not relative:
        df = tax.iloc[:, 0].replace(".__", "", regex=True).str.split(';', expand=False)
        df = df.apply(lambda x: pad_or_truncate(x, 8)) #ensure every taxonomic level has empty spaces for levels that are not represented
        df = pd.DataFrame(df.to_list(), columns=levels)
        df = df.merge(tax.iloc[: , 1:], left_index=True, right_index=True, how='inner')
        df = df.replace("", pd.NA)
        df = df.replace("_", " ", regex=True)
        df["occurences"] = df.apply(lambda x: x.isna().sum(), axis=1) #Occurences of <NA>

        ##current method for knowing if tax file contains prefix-less species 
        ##possibly should add some methodology e.g. if almost no genus appears in species name then we move on
        ##might need to be stand alone function if it appears in combined files as well
        naming_convention = True
        for i, entry in df.iterrows(): #are Genus and Species combined already?
            if entry["Genus"] == "": 
                pass
            elif entry["Genus"] in entry["Species"]: #if *a* genus is found in a species entry we say yes
                naming_convention = False
                break

        df["Species"] = df["Genus"].str.strip() + "_" + df["Species"].str.strip()
        #df["Species"] = df["Species"].str.replace("(.*?)_(?!\S)", "", regex=True)

    else:
        df = tax.iloc[:, 1].replace(".__", "", regex=True).str.split(';', expand=False) #access second column

        if len(df.iloc[0]) == 7: #if no strain column
            for taxa in df:
                taxa.append("")

        df = pd.DataFrame(df.to_list(), columns=levels)

        for level in levels:
            df[level] = df[level].str.strip() #remove potential whitespace

        df = df.fillna("") #eliminate noneTypes that appear
        
        ##current method for knowing if tax file contains prefix-less species 
        ##possibly should add some methodology e.g. if almost no genus appears in species name then we move on
        ##might need to be stand alone function if it appears in combined files as well
        naming_convention = True
        for i, entry in df.iterrows(): #are Genus and Species combined already?
            if entry["Genus"] == "": 
                pass
            elif entry["Genus"] in entry["Species"]: #if *a* genus is found in a species entry we say yes
                naming_convention = False
                break

        if naming_convention:
            df["Species"] = df["Genus"] + "_" + df["Species"]
            df["Species"] = df["Species"].str.replace("(.*?)_(?!\S)", "", regex=True) #replace everything before and up to "_" if nothing after

        df["Index"] = tax.iloc[:, 0]
        df = df.set_index("Index")
        df = df.merge(feat, left_index=True, right_index=True, how='inner')
        df = df.replace("", pd.NA) #Non applicable instead of empty string
        
    ##currently not adjusting on genus level -- rational: we do not always know if genus part of species name is true genus name    
    df = df.set_index("Species") 
    df = df.rename(index=homosynonyms) #homosynonym rename on species level
    df = df.reset_index()

    kingdom_df = get_grouped_tax_level("Kingdom", df)
    phylum_df = get_grouped_tax_level("Phylum", df)
    class_df = get_grouped_tax_level("Class", df)
    order_df = get_grouped_tax_level("Order", df)
    family_df = get_grouped_tax_level("Family", df)
    genus_df = get_grouped_tax_level("Genus", df)
    species_df = get_grouped_tax_level("Species", df)
    strain_df = get_grouped_tax_level("Strain", df)

    return df, kingdom_df, phylum_df, class_df, order_df, family_df, genus_df, species_df, strain_df