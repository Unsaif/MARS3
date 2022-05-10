
def visualise(overview, stats):

    # EVEN IF WE HAVE GROUPS, THESE MENTIONED BELOW ARE JUST TAKEN AS 1 GROUP AS GENERAL FIGURE OVERVIEW OF ALL THE DATA

    # Always make boxplots of
    # Total reads

    # Together
    # % Reads that have species/genus
    # % Agora2reads on total reads
    # End togheterh


    # Taxa columns 3 boxplots in 1 fig.

    # File with sum of rel .abun after 1st norm and cutoff that are <0.95

    # All ratio fir/bac in one figure (4 boxplots)

    # This would than in a figure be a boxplot for each group for the sig value
    significant_stats = stats.copy()
    # select all rows with p-value <0.05. (always the first column)
    # Create boxplots with for all sig values based on the groups

    return