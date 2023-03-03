import matplotlib.pyplot as plt

def visualise(overview, level, stats=None):

    # EVEN IF WE HAVE GROUPS, THESE MENTIONED BELOW ARE JUST TAKEN AS 1 GROUP AS GENERAL FIGURE OVERVIEW OF ALL THE DATA

    # Always make boxplots of
    # Total reads
    if stats is not None:
        x = 2
    else:
        plt.boxplot(overview.totalReads)

    plt.ylabel('total reads')
    plt.savefig('MARS_output/Figs/TotalReadsBoxplot.png')

    plt.subplot(1, 2, 1)
    plt.boxplot(overview.loc[:, 'percentage' + level + 'Reads'])
    plt.title(f'Fraction {level}-level reads')
    plt.xticks([1], ['all data'])

    plt.subplot(1, 2, 2)
    plt.boxplot(overview.loc[:, 'percentageAgoraTotalReads'])
    plt.title(f'Fraction {level}-level \n AGORA2 reads')
    plt.xticks([1], ['all data'])
    plt.savefig('MARS_output/Figs/coverageBoxPlots.png')

    # Taxa columns 3 boxplots in 1 fig.

    # All ratio fir/bac in one figure (4 boxplots) DONE
    plt.subplot(2, 2, 1)
    plt.boxplot(overview.loc[:, 'FBRatioTotal'])
    plt.title('Total reads')
    plt.xticks([1], ['all data'])

    plt.subplot(2, 2, 2)
    plt.boxplot(overview.loc[:, f'FBRatio{level}'])
    plt.title(f'{level} Reads')
    plt.xticks([1], ['all data'])

    plt.subplot(2, 2, 3)
    plt.boxplot(overview.loc[:, 'FBRatioAgora'])
    plt.title('Agora reads')
    plt.xticks([1], ['all data'])

    plt.subplot(2, 2, 4)
    plt.boxplot(overview.loc[:, 'FBRatioAfterCutOff'])
    plt.title('After Cut-Off')
    plt.xticks([1], ['all data'])

    plt.suptitle('Firmicutes/Bacteroidetes Ratio')
    plt.savefig('MARS_output/Figs/FirmicutesBacteroidetesRatios.png')

    # Figure of additional taxa if statement

    # This would than in a figure be a boxplot for each group for the sig value
    if stats != None:
        significant_stats = stats.copy()
        # select all rows with p-value <0.05. (always the first column)
        # Create boxplots with for all sig values based on the groups

    return
