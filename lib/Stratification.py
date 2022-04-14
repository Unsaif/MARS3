import pandas as pd
import scipy.stats as stats


def split_df(df, stratpath):

    if stratpath.endswith('.xlsx'):
        strat = pd.read_excel(stratpath, index_col=0)
    else:
        strat = pd.read_csv(stratpath)
    strat.index = strat.index.astype(str)

    comb = pd.merge(df, strat, left_index=True, right_index=True)

    groups = comb.iloc[:,-1].unique()

    init = 1
    stat_df = 0

    for group in groups:
        sub_df = comb.loc[comb.iloc[:, -1] == group].copy()
        sub_df.drop(sub_df.columns[[-1]], axis=1, inplace=True)

        avg = sub_df.mean()
        avg = avg.rename(f'Average_{group}')

        std = sub_df.std()
        std = std.rename(f'Std_dev_{group}')

        if init:
            stat_df = pd.merge(avg, std, left_index=True, right_index=True)
            init = 0
        else:
            stat_df = pd.merge(stat_df, avg, left_index=True, right_index=True)
            stat_df = pd.merge(stat_df, std, left_index=True, right_index=True)

    p_val_list = []
    if len(groups) == 2:

        base_group = groups[0]
        test_group = groups[1]

        df_base = comb.loc[comb.iloc[:, -1] == base_group].copy()
        df_base.drop(df_base.columns[[-1]], axis=1, inplace=True)

        df_test = comb.loc[comb.iloc[:, -1] == test_group].copy()
        df_test.drop(df_test.columns[[-1]], axis=1, inplace=True)


        for column in df_base:
            base = df_base[column]
            test = df_test[column]

            base_a = base.to_numpy()
            test_a = test.to_numpy()

            t_stat, p_value = stats.ttest_ind(base_a, test_a)
            p_val_list += [p_value]

    elif len(groups) > 2:
        #ANOVA

        sub_df = comb.loc[comb.iloc[:, -1] == groups[0]].copy()
        sub_df.drop(sub_df.columns[[-1]], axis=1, inplace=True)
        for column in sub_df:
            data_list = []
            for group in groups:
                spef_df = comb.loc[comb.iloc[:, -1] == group].copy()
                spef_data = spef_df[column]
                data_list += [spef_data]

            f_anova, p_anova = stats.f_oneway(*data_list)
            p_val_list += [p_anova]

    else:
        pass

    # BONJOUR
    p_series = pd.Series(p_val_list, name='P_values', index=stat_df.index)
    stat_df = stat_df.join(p_series)

    cols = list(stat_df.columns)
    cols = [cols[-1]] + cols[:-1]
    stat_df = stat_df[cols]
    print(stat_df)

    # Save file here 
