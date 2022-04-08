
from lib import species_genus_association, agora_checking

def pipeline(df, total_df, levels, level, agora2_level_set, agora2_species, agora2_genera):
    associated_species, associated_genus = species_genus_association.association(df, levels, level)

    absent, present = agora_checking.agora_checking(total_df, agora2_level_set)

    total_df_species_in_agora = df[df["Species"].isin(agora2_species)]
    total_df_genera_in_agora = df[df["Genus"].isin(agora2_genera)]

    associated_species_agora2, _ = species_genus_association.association(total_df_species_in_agora, levels, level)
    _, associated_genus_agora2 = species_genus_association.association(total_df_genera_in_agora, levels, level)

    #save these files
    total_df.to_csv(f'MARS_output/total_{level}.csv')
    associated_species.to_csv(f'MARS_output/associated_species_{level}.csv')
    associated_genus.to_csv(f'MARS_output/associated_genus_{level}.csv')
    absent.to_csv(f'MARS_output/absent_{level}.csv')
    present.to_csv(f'MARS_output/present_{level}.csv')
    associated_species_agora2.to_csv(f'MARS_output/associated_species_agora2_{level}.csv')
    associated_genus_agora2.to_csv(f'MARS_output/associated_genus_agora2_{level}.csv')