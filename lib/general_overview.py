import pandas as pd

def general_overview(initial_df, list_phylum_df, list_species_genus_dfs, level, extra_phyla=None):

    """ Takes the dataframe from general_overview and the user-defined path to the stratification file to run general
       statistics on the different values in the general_overview dataframe. It saves the statistics dataframe as a .csv
       at the end of the script

       Parameters
       ----------
       initial_df : pandas dataframe, required
           dataframe from general_overview function
       list_phylum_df: string, required
           path to the stratification file
       list_species_genus_dfs
       extra_phyla

       Returns
       -------
       None

       Authors
       -------
       Bram Nap - 04/2022
       """

    # TODO: Add loss after cut and renorm for bac, firm, and ratio
    total_phylum, associated_phylum, associated_agora_phylum = list_phylum_df
    reads_afteragora_df, reads_beforeagora_df, normalised_cutoff_df = list_species_genus_dfs

    sum_initial_df = sum_rename_sort(initial_df.iloc[:, 8:], 'Total Reads')

    sum_associated_reads = sum_rename_sort(associated_phylum, 'Associated Reads')
    final_df = pd.merge(sum_initial_df, sum_associated_reads, left_index=True, right_index=True)

    sum_associated_agora_reads = sum_rename_sort(associated_agora_phylum, 'AGORA2 associated reads')
    final_df = pd.merge(final_df, sum_associated_agora_reads, left_index=True, right_index=True)

    associated_reads_coverage = sum_associated_reads/sum_initial_df
    associated_reads_coverage = associated_reads_coverage.rename('% Reads that have species')
    final_df = pd.merge(final_df, associated_reads_coverage, left_index=True, right_index=True)

    agora_total_reads_coverage = sum_associated_agora_reads / sum_initial_df
    agora_total_reads_coverage = agora_total_reads_coverage.rename('% Agora2 total reads coverage')
    final_df = pd.merge(final_df, agora_total_reads_coverage, left_index=True, right_index=True)

    agora_associated_reads_coverage = sum_associated_agora_reads / sum_associated_reads
    agora_associated_reads_coverage = agora_associated_reads_coverage.rename('% Agora2 associated reads coverage')
    final_df = pd.merge(final_df, agora_associated_reads_coverage, left_index=True, right_index=True)

    before_agora_taxa = reads_beforeagora_df[reads_beforeagora_df > 0].count()
    before_agora_taxa = before_agora_taxa.rename('# taxa before agora2 mapping')
    final_df = pd.merge(final_df, before_agora_taxa, left_index=True, right_index=True)

    after_agora_taxa = reads_afteragora_df[reads_afteragora_df > 0].count()
    after_agora_taxa = after_agora_taxa.rename('# taxa after agora2 mapping')
    final_df = pd.merge(final_df, after_agora_taxa, left_index=True, right_index=True)

    after_agora_normcut_taxa = normalised_cutoff_df[normalised_cutoff_df > 0].count()
    after_agora_normcut_taxa = after_agora_normcut_taxa.rename('# taxa after agora2 mapping and normalizing')
    final_df = pd.merge(final_df, after_agora_normcut_taxa, left_index=True, right_index=True)

    loss_to_cutoff = sum_rename_sort(normalised_cutoff_df, 'Sum of rel. abund. 1st norm')
    final_df = pd.merge(final_df, loss_to_cutoff, left_index=True, right_index=True)

    agora_cutoff_total_coverage = agora_total_reads_coverage*loss_to_cutoff
    agora_cutoff_total_coverage = agora_cutoff_total_coverage.rename('Coverage agora2 total reads after norm and cutoff')
    final_df = pd.merge(final_df, agora_cutoff_total_coverage, left_index=True, right_index=True)

    agora_cutoff_associated_coverage = agora_associated_reads_coverage*loss_to_cutoff
    agora_cutoff_associated_coverage = agora_cutoff_associated_coverage.rename('Coverage agora2 associated reads after norm and cutoff')
    final_df = pd.merge(final_df, agora_cutoff_associated_coverage, left_index=True, right_index=True)

    total_phylum_of_interest = find_phylum_reads(total_phylum, associated_phylum, associated_agora_phylum, 'Bacteroidetes')

    total_phylum_of_interest = find_phylum_reads(total_phylum, associated_phylum, associated_agora_phylum, 'Firmicutes',
                                                 total_phylum_of_interest)

    if extra_phyla is not None:
        for phyl in extra_phyla:
            total_phylum_of_interest = find_phylum_reads(total_phylum, associated_phylum, associated_agora_phylum,
                                                total_phylum_of_interest, phyl)
    else:
        pass

    final_df = pd.merge(final_df, total_phylum_of_interest, left_index=True, right_index=True)

    ratio_df = ratio_calc(final_df)

    final_df = pd.merge(final_df, ratio_df, left_index=True, right_index=True)
    final_df.to_csv(f'MARS_output/general_statistics_{level}.csv')

    return final_df


def sum_rename_sort(df, name):
    df_sum = df.sum(axis=0)
    df_sumname = df_sum.rename(name)
    df_sumname = df_sumname.sort_index()

    return df_sumname


def find_phylum_reads(total_reads_df, associated_reads_df, agora_reads_df, phylum, final_df=None):

    phylum = phylum.lower()
    phylum = phylum.capitalize()

    try:
        total_phylum_reads = total_reads_df.loc[[phylum]]
        associated_phylum_reads = associated_reads_df.loc[[phylum]]
        agora_phylum_reads = agora_reads_df.loc[[phylum]]

        reads_loss_due_assocation = (total_phylum_reads-associated_phylum_reads)/total_reads_df
        reads_loss_due_mapping = (associated_phylum_reads-agora_phylum_reads)/total_phylum_reads

        combined_df = pd.concat([total_phylum_reads, associated_phylum_reads, agora_phylum_reads, reads_loss_due_assocation, reads_loss_due_mapping], axis=0)

        combined_df = combined_df.transpose()
        combined_df.columns = ['Total reads ' + phylum], ['Associated reads ' + phylum], \
                          ['Associated reads after agora ' + phylum], ['loss reads due to unassociation ' + phylum], \
                          ['loss reads due to agora2 mapping ' + phylum]

        if final_df is not None:
            final_df = pd.merge(final_df, combined_df, left_index=True, right_index=True)
        else:
            final_df = combined_df

    except KeyError:
        print(f'You have misspelled {phylum} or it is not present, please check the spelling or remove')

    return final_df


def ratio_calc(final_df):
    # F/B ratio

    # TODO: change to explicit names

    bacter_total_reads = final_df.iloc[:, 12]
    firmic_total_reads = final_df.iloc[:, 17]

    bacter_associated_reads = final_df.iloc[:, 13]
    firmic_associated_reads = final_df.iloc[:, 18]

    bacter_agora_reads = final_df.iloc[:, 14]
    firmic_agora_reads = final_df.iloc[:, 19]

    total_ratio = firmic_total_reads/bacter_total_reads
    total_ratio = total_ratio.rename('Total Fir/Bac ratio')

    associated_ratio = firmic_associated_reads/bacter_associated_reads
    associated_ratio = associated_ratio.rename('Associated Fir/Bac ratio')

    agora_ratio = firmic_agora_reads/bacter_agora_reads
    agora_ratio = agora_ratio.rename('Agora Fir/Bac ratio')

    comb_ratio = pd.merge(total_ratio, associated_ratio, left_index=True, right_index=True)
    comb_ratio = pd.merge(comb_ratio, agora_ratio, left_index=True, right_index=True)

    return comb_ratio
