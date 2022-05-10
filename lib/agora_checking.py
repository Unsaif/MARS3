def agora_checking(df_level, agora2_level_set):
    present = []
    absent = []
    for taxa in df_level.index:
        # print(taxa)
        if taxa.strip() in agora2_level_set or taxa in agora2_level_set:
            present.append(taxa)
        else:
            absent.append(taxa)
            
    return df_level.loc[absent], df_level.loc[present]