import pandas as pd

def sanityChecker(initial_df, species_df):

    totalSpecies = species_df[species_df > 0].count()
    noSpecies = totalSpecies[totalSpecies == 0]
    noSpecies = noSpecies.replace(0, 1)
    noSpecies = noSpecies.rename('noSpecies')

    fewSpecies = totalSpecies[totalSpecies.between(1, 10)]
    fewSpecies = fewSpecies.rename('fewSpecies')

    dfmerged = pd.concat([noSpecies, fewSpecies],axis=1)
    dfmerged = dfmerged.fillna(0)

    # Check the reads that have a lower amount of 3x stdev of the mean.
    summedReads = initial_df[:, 8:].sum()

    with open('Samples2LookAt.txt', 'w') as file:
        for i, row in dfmerged.iterrows():
            ID = i
            if dfmerged.iloc[i, 0] == 1:
                file.write(f'Sample {ID} does not seem to have any species associated reads. It can not be used to '
                           f'create panSpecies microbiome models\n')
            if dfmerged.iloc[i, 1] != 0:
                file.write(f'Sample {ID} does seems to have reads for 10 or less species. This can impact the '
                           f'reliability of the panSpecies microbiome models')







