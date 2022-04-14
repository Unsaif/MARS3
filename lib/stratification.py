import pandas as pd
import scipy.stats as stats


# TODO: change general_info_df to maybe the new name of general_stats ? CHANGE TO GENERAL OVERVIEW
def general_statistics_on_groups(general_info_df, path_to_stratification_file):

    """ Takes the dataframe from general_overview and the user-defined path to the stratification file to run general
    statistics on the different values in the general_overview dataframe. It saves the statistics dataframe as a .csv
    at the end of the script

    Parameters
    ----------
    general_info_df : pandas dataframe, required
        dataframe from general_overview function
    path_to_stratification_file: string, required
        path to the stratification file

    Returns
    -------
    None

    Authors
    -------
    Bram Nap - 04/2022
    """

    # Based on the extension of the stratification path, read the stratification file
    # TODO: Maybe use a try statement to throw an error when there is no file extension in the string name?
    if path_to_stratification_file.endswith('.xlsx'):
        stratification_df = pd.read_excel(path_to_stratification_file, index_col=0)
    else:
        stratification_df = pd.read_csv(path_to_stratification_file)
    # Set the index of the df to strings, required for merging with the general_overview df
    stratification_df.index = stratification_df.index.astype(str)

    # Combine the general_overview_df and the stratification_df
    combined_df = pd.merge(general_info_df, stratification_df, left_index=True, right_index=True)

    # Find all the unique groups in the stratification column
    groups = combined_df.iloc[:, -1].unique()

    # Needed to make sure the first df merge is carried out
    first_iteration = 1
    # Required as placeholder
    statistics_df = 0

    # AVERAGES AND STANDARD DEVIATIONS
    for group in groups:
        # Create a sub_df for each unique group
        sub_df = combined_df.loc[combined_df.iloc[:, -1] == group].copy()
        # Drop the group information column, required to avoid warning during statistics
        sub_df.drop(sub_df.columns[[-1]], axis=1, inplace=True)

        # Calculate the averages and standard deviations
        avg = sub_df.mean()
        avg = avg.rename(f'Average_{group}')

        std = sub_df.std()
        std = std.rename(f'Std_dev_{group}')

        # If no statistics_df has been made yet
        if first_iteration:
            # Merge the average and standard devations
            statistics_df = pd.merge(avg, std, left_index=True, right_index=True)
            # As first merge has been carried out, set to 0
            first_iteration = 0
        else:
            # Merge averages then standard deviations of subsequent groups the the statistics df
            statistics_df = pd.merge(statistics_df, avg, left_index=True, right_index=True)
            statistics_df = pd.merge(statistics_df, std, left_index=True, right_index=True)

    # BASIC STATISTICS

    # Initialise a list to store values in
    pvalue_list = []

    # If there are only 2 groups, t-test
    if len(groups) == 2:

        # Obtain the group names
        group1 = groups[0]
        group2 = groups[1]

        # Create dataframes per group, drop the last column as it only contains group info
        group1_df = combined_df.loc[combined_df.iloc[:, -1] == group1].copy()
        group1_df.drop(group1_df.columns[[-1]], axis=1, inplace=True)

        group2_df = combined_df.loc[combined_df.iloc[:, -1] == group2].copy()
        group2_df.drop(group2_df.columns[[-1]], axis=1, inplace=True)

        # For each column (general info) in the dataframes
        for column in group1_df:
            # Obtain group specific values
            group1_col = group1_df[column]
            group2_col = group2_df[column]

            # Translate into array type, required for t-test
            group1_array = group1_col.to_numpy()
            group2_array = group2_col.to_numpy()

            # Run t-test and save the p_value
            t_value, p_value = stats.ttest_ind(group1_array, group2_array)
            pvalue_list += [p_value]

    # If there are more than 2 groups, perform ANOVA
    elif len(groups) > 2:

        # Create a sub_df for the first group that is used to loop over the columns for
        sub_df = combined_df.loc[combined_df.iloc[:, -1] == groups[0]].copy()
        sub_df.drop(sub_df.columns[[-1]], axis=1, inplace=True)

        # For each column (general info) in the dataframe
        for column in sub_df:
            # Initialse list to store data in
            data_list = []

            # For each group
            for group in groups:
                # Create a group specific dataframe, obtain data from the relevant column and store in a list
                specific_df = combined_df.loc[combined_df.iloc[:, -1] == group].copy()
                specific_data = specific_df[column]
                data_list += [specific_data]

            # Run ANOVA and store the p-value, * is used to deconstruct the list into separate variables
            f_anova, p_anova = stats.f_oneway(*data_list)
            pvalue_list += [p_anova]

    else:
        pass

    # Create a series from the list, with the same indexes as the statistics_df
    p_series = pd.Series(pvalue_list, name='p_values', index=statistics_df.index)
    statistics_df = statistics_df.join(p_series)

    # Change it so that the p-values are the first column
    cols = list(statistics_df.columns)
    cols = [cols[-1]] + cols[:-1]
    stat_df = statistics_df[cols]

    # TODO:Save file here
