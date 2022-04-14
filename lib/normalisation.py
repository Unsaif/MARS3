def normalise_and_cut(present_df, level):

    total_reads = present_df.sum()
    agora_normed = present_df.loc[:].div(total_reads)

    agora_normed_cut = agora_normed.copy()
    agora_normed_cut[agora_normed_cut < 1e-5] = 0

    # Renormalize
    total_rel_abund = agora_normed_cut.sum()
    agora_renormed = agora_normed_cut.loc[:].div(total_rel_abund)

    #TODO: saving files 
    # save agora normed 
    # agora normed cut -not mgpipe ready
    # save agora renormed 
    agora_normed.to_csv(f'MARS_output/agora_normed_{level}.csv')
    agora_normed_cut.to_csv(f'MARS_output/agora_normed_cut_{level}.csv')
    agora_renormed.to_csv(f'MARS_output/agora_renormed_cut_{level}.csv')

    return agora_normed_cut
