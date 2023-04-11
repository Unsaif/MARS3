import pandas as pd
import warnings
import numpy as np


def general_overview(initial_df, list_phylum_df, list_species_genus_dfs, level, extra_phyla=None):

    """ Takes information from all previously created dataframes and gives an overview of general values of the data
    These include: reads, coverages, loss of reads and specific phylum information

       Parameters
       ----------
       initial_df : pandas dataframe, required
           reads dataframe with no changes
       list_phylum_df : list, required
            list with variables containing various phylum-level information dataframes
       list_species_genus_dfs : list, required
            list with variables containing various genus/species - level information dataframes
       level : string, required
            species or genus, used as file name to distinguish files from each other
       extra_phyla : list, optional | defaults to None
            list with the names as strings of phyla other than bacteroidetes or firmicutes the user wants information on
       Returns
       -------
       None

       Authors
       -------
       Bram Nap - 04/2022
       """

    # Unpack the lists
    total_phylum, associated_phylum, associated_agora_phylum = list_phylum_df
    reads_afteragora_df, reads_beforeagora_df, normalised_cutoff_df, pan_phylum_df = list_species_genus_dfs

    # Calculate the sum of reads for each sample of respective dfs and add join them
    sum_initial_df = sum_rename_sort(initial_df.iloc[:, 8:], 'tot_Reads')

    sum_associated_reads = sum_rename_sort(associated_phylum, level + '_Reads')
    final_df = pd.merge(sum_initial_df, sum_associated_reads, left_index=True, right_index=True)

    sum_associated_agora_reads = sum_rename_sort(associated_agora_phylum, 'AGORA2_' + level + '_Reads')
    final_df = pd.merge(final_df, sum_associated_agora_reads, left_index=True, right_index=True)

    # Calculate the coverages for each sample and add to the final_df variable
    associated_reads_coverage = sum_associated_reads / sum_initial_df
    associated_reads_coverage = associated_reads_coverage.rename(level + '_Reads_FracTot')
    final_df = pd.merge(final_df, associated_reads_coverage, left_index=True, right_index=True)

    agora_total_reads_coverage = sum_associated_agora_reads / sum_initial_df
    agora_total_reads_coverage = agora_total_reads_coverage.rename('AGORA2_'+level+'_Reads_FracTot')
    final_df = pd.merge(final_df, agora_total_reads_coverage, left_index=True, right_index=True)

    agora_associated_reads_coverage = sum_associated_agora_reads / sum_associated_reads
    agora_associated_reads_coverage = agora_associated_reads_coverage.rename('AGORA2_' + level + '_Reads_Frac'+level)
    final_df = pd.merge(final_df, agora_associated_reads_coverage, left_index=True, right_index=True)

    # Obtain the amount of the species / genera that have, at least, 1 read for each sample and add to final_df variable
    before_agora_taxa = reads_beforeagora_df[reads_beforeagora_df > 0].count()
    before_agora_taxa = before_agora_taxa.rename('Tot_'+level)
    final_df = pd.merge(final_df, before_agora_taxa, left_index=True, right_index=True)

    after_agora_taxa = reads_afteragora_df[reads_afteragora_df > 0].count()
    after_agora_taxa = after_agora_taxa.rename('ABORA2_' + level)
    final_df = pd.merge(final_df, after_agora_taxa, left_index=True, right_index=True)

    after_agora_normcut_taxa = normalised_cutoff_df[normalised_cutoff_df > 0].count()
    after_agora_normcut_taxa = after_agora_normcut_taxa.rename('AGORA2_Cutoff_'+level)
    final_df = pd.merge(final_df, after_agora_normcut_taxa, left_index=True, right_index=True)

    # Calculate the total rel. abundance for each sample after cutoff and add to final_df variable
    loss_to_cutoff = sum_rename_sort(normalised_cutoff_df, 'AGORA2_Cutoff_RelAb')
    final_df = pd.merge(final_df, loss_to_cutoff, left_index=True, right_index=True)

    # Calculate the effect of cutoff on coverages and add to final_df variable
    agora_cutoff_total_coverage = agora_total_reads_coverage*loss_to_cutoff
    agora_cutoff_total_coverage = agora_cutoff_total_coverage.rename('AGORA2_Cutoff_'+level+'Reads_FracTot')
    final_df = pd.merge(final_df, agora_cutoff_total_coverage, left_index=True, right_index=True)

    agora_cutoff_associated_coverage = agora_associated_reads_coverage*loss_to_cutoff
    agora_cutoff_associated_coverage = agora_cutoff_associated_coverage.rename('AGORA2_Cutoff_'+level+'Reads_Frac'+level)
    final_df = pd.merge(final_df, agora_cutoff_associated_coverage, left_index=True, right_index=True)

    # Total amount of reads for the phylum of bacteroidetes for each sample in different dataframes
    total_phylum_of_interest = find_phylum_reads(total_phylum, associated_phylum, associated_agora_phylum, 'Bacteroidetes', level)

    # Total amount of reads for the phylum of firmicutes for each sample in different dataframes
    total_phylum_of_interest = find_phylum_reads(total_phylum, associated_phylum, associated_agora_phylum, 'Firmicutes',
                                                 level, total_phylum_of_interest)

    # Total amount of reads for user defined phyla for each sample in different dataframes
    if extra_phyla is not None:
        for phyl in extra_phyla:
            total_phylum_of_interest = find_phylum_reads(total_phylum, associated_phylum, associated_agora_phylum,
                                                         phyl, level, total_phylum_of_interest)
    else:
        pass

    # add phylum information dataframe to the final_df variable
    final_df = pd.merge(final_df, total_phylum_of_interest, left_index=True, right_index=True)

    # Calculate the fir/bac ratio for each sample and add it to the final_df variable
    ratio_df = ratio_calc(final_df, level)
    final_df = pd.merge(final_df, ratio_df, left_index=True, right_index=True)

    # Calculate the fir/bac ratio after mapping, cutoff and renormalisation
    final_ratio = pan_phylum_df.loc['Firmicutes']/pan_phylum_df.loc['Bacteroidetes']
    final_ratio = final_ratio.rename('FB_Ratio_AGORA2_Cutoff')
    final_df = pd.merge(final_df, final_ratio, left_index=True, right_index=True)

    # Save the final_df variable as a .csv file
    if level == 'Sp':
        fol = 'Species'
    elif level == 'Ge':
        fol = 'Genus'
    final_df.to_csv(f'MARS_output/{fol}/general_overview_{fol}.csv')

    return final_df


def sum_rename_sort(df, name):

    """ Takes a dataframe, and sums all the columns. The resulting Series is given a name, is sorted based on the index
    names and returned

          Parameters
          ----------
          df : pandas dataframe, required
              Dataframe with summable entries
          name : string, required
               The name of the Series resulting from the summation

          Returns
          -------
          df_sumname : pandas series
            The summed, renamed and sorted series based on the input dataframe

          Authors
          -------
          Bram Nap - 04/2022
          """

    df_sum = df.sum(axis=0)
    df_sumname = df_sum.rename(name)
    df_sumname = df_sumname.sort_index()

    return df_sumname


def find_phylum_reads(total_reads_df, associated_reads_df, agora_reads_df, phylum, level, final_df=None):

    """ Finds out the total amount of reads for a given phylum for different dataframes and calculates where and how much
    information is lost

          Parameters
          ----------
          total_reads_df : pandas dataframe, required
              Dataframe with total phylum reads
          associated_reads_df : pandas dataframe, required
              Dataframe with total phylum reads that are associated to a genus or species
          agora_reads_df : pandas dataframe, required
              Dataframe with total phylum reads that are associated to a AGORA2 mapped genus or species
          phylum : string, required
              The phylum for which the summary is made
          final_df : pandas dataframe, optional
              A dataframe that can be used to add subsequent dataframes created in this function to. Usually is the
              output of this function when re-called

          Returns
          -------
          final_df : pandas dataframe
            overview of the phylum in question in combination with the final_df if that is given as input

          Authors
          -------
          Bram Nap - 04/2022
          """

    # Make sure that format of the phylum name will correspond to the phyla names in the dataframe
    
    phylum = phylum.lower()
    phylum = phylum.capitalize()
    
    # A try statement to catch misspelled or non-existent phyla names
    try:
        # Obtain relevant data from the respective dataframes
        total_phylum_reads = total_reads_df.loc[phylum]
        associated_phylum_reads = associated_reads_df.loc[phylum]
        agora_phylum_reads = agora_reads_df.loc[phylum]
        
        # Calculate the fractions lost due to no association and agora2 mapping
        reads_loss_due_assocation = (total_phylum_reads-associated_phylum_reads)/total_reads_df
        reads_loss_due_mapping = (associated_phylum_reads-agora_phylum_reads)/total_phylum_reads

        # Combine different Series into one dataframe
        frame = { 
            f'{phylum}_Tot_Reads': total_phylum_reads,
            f'{phylum}_{level}_Reads': associated_phylum_reads,
            f'{phylum}_AGORA2_{level}_Reads': agora_phylum_reads,
            f'{phylum}_Reads_no{level}': reads_loss_due_assocation.loc[phylum].T,
            f'{phylum}_{level}_Reads_noAGORA2': reads_loss_due_mapping
            }

        combined_df = pd.DataFrame(frame)

        # If there is final_df input, combine the two dataframes
        if final_df is not None:
            final_df = pd.merge(final_df, combined_df, left_index=True, right_index=True)

        # If there is no final_df input, final_df is the combined_df
        else:
            final_df = combined_df

    # If there is no exact match in the dataframe raise a warning
    except KeyError:
        warnings.warn(f'You have misspelled {phylum} or it is not present, please check the spelling or remove')

    return final_df


def ratio_calc(final_df, level):
    """ Calculates the firmicutes/bacteroidetes ratios

          Parameters
          ----------
          final_df : pandas dataframe, required
              Takes the complete dataframe from the general_overview function

          Returns
          -------
          comb_ratio : pandas dataframe
              A dataframe from all the ratio calculation Series created in this function

          Authors
          -------
          Bram Nap - 04/2022
          """

    # Select the relevant data
    bacter_total_reads = final_df.loc[:, 'Bacteroidetes_Tot_Reads']
    firmic_total_reads = final_df.loc[:, 'Firmicutes_Tot_reads']

    bacter_associated_reads = final_df.loc[:, f'Bacteroidetes_{level}_Reads']
    firmic_associated_reads = final_df.loc[:, f'Firmicutes_{level}_Reads']

    bacter_agora_reads = final_df.loc[:, f'Bacteroidetes_AGORA2_{level}_Reads']
    firmic_agora_reads = final_df.loc[:, f'Firmicutes_AGORA2_{level}_Reads']

    # Calculate the ratios
    total_ratio = firmic_total_reads/bacter_total_reads
    total_ratio = total_ratio.rename('FB_Ratio_tot')

    associated_ratio = firmic_associated_reads/bacter_associated_reads
    associated_ratio = associated_ratio.rename('FB_Ratio_' + level)

    agora_ratio = firmic_agora_reads/bacter_agora_reads
    agora_ratio = agora_ratio.rename('FB_Ratio_AGORA2')

    # Adds to one single dataframe
    comb_ratio = pd.merge(total_ratio, associated_ratio, left_index=True, right_index=True)
    comb_ratio = pd.merge(comb_ratio, agora_ratio, left_index=True, right_index=True)

    return comb_ratio
