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
    with open('Samples2LookAt.txt', 'w') as file:
        file.write('hey')







