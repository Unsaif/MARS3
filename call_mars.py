from main import main

# taxonomy_table = r'/mnt/c/Users/THuls/Downloads/Silva16SU19_taxonomy.tsv'
# feature_table = r'/mnt/c/Users/THuls/Downloads/Silva16SU19test.txt'

taxonomy_table = r'C:\Users\MSPG\Documents\GitHub\MARS3\testFiles\taxonomy.tsv'
feature_table = r'C:\Users\MSPG\Documents\GitHub\MARS3\testFiles\feature-table.tsv'
path_to_stratification_file = r'/mnt/c/Users/THuls/Documents/python_projects/Test/files/strat_mult_group.xlsx'

main(taxonomy_table=taxonomy_table, feature_table=feature_table, path_to_stratification_file=None)