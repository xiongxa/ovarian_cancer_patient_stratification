import pandas as pd
import os
from lifelines import CoxPHFitter
from config.constants import CLINICAL_CHARA, IMMU_CHARA
from sklearn.model_selection import train_test_split
import copy
# 获取环境变量
PARP_DIR = os.environ.get('PARP_DIR')

def splite_dataframe(df, split_ratio, line_type, immu_cand, bio_cand, cli_cand):
    # df_shuffled = df.sample(frac=1, random_state=18).reset_index(drop=True)
    # df_shuffled = df.sample(frac=1, random_state=3).reset_index(drop=True)  # almost ok
    df_shuffled = df.sample(frac=1, random_state=99).reset_index(drop=True)  # almost ok
    train_size = int(len(df_shuffled) * split_ratio)
    # 分割 DataFrame
    df_train = df_shuffled[:train_size]
    df_test = df_shuffled[train_size:]
    
    default_features = ['PFS(month)', 'PFS_type', 'OS(month)', 'OS_type']
    df_train_gen = remove_missing_data(df_train[['BRCA_HRD status']+default_features])
    df_train_immu = remove_missing_data(df_train[immu_cand+default_features])
    df_train_bio = remove_missing_data(df_train[bio_cand+default_features])
    df_train_cli = remove_missing_data(df_train[cli_cand+default_features])
    df_train_combine = remove_missing_data(df_train[['BRCA_HRD status']+immu_cand+bio_cand+cli_cand+default_features])
    
    df_test_gen = remove_missing_data(df_test[['BRCA_HRD status']+default_features])
    df_test_immu = remove_missing_data(df_test[immu_cand+default_features])
    df_test_bio = remove_missing_data(df_test[bio_cand+default_features])
    df_test_cli = remove_missing_data(df_test[cli_cand+default_features])
    df_test_combine = remove_missing_data(df_test[['BRCA_HRD status']+immu_cand+bio_cand+cli_cand+default_features])
    
    print(F'Dataframe train shape: {df_train.shape} {df_train_gen.shape}, {df_train_immu.shape}, {df_train_bio.shape}, {df_train_cli.shape}, {df_train_combine.shape}')
    print(F'Dataframe test shape: {df_test.shape} {df_test_gen.shape}, {df_test_immu.shape}, {df_test_bio.shape}, {df_test_cli.shape}, {df_test_combine.shape}')
    
    df_train_gen.to_excel(F'{PARP_DIR}/feature_selection/data/datasets/{line_type}_gen_dataframe_train.xlsx', index=False)
    df_train_immu.to_excel(F'{PARP_DIR}/feature_selection/data/datasets/{line_type}_immu_dataframe_train.xlsx', index=False)
    df_train_bio.to_excel(F'{PARP_DIR}/feature_selection/data/datasets/{line_type}_bio_dataframe_train.xlsx', index=False)
    df_train_cli.to_excel(F'{PARP_DIR}/feature_selection/data/datasets/{line_type}_clinical_dataframe_train.xlsx', index=False)
    df_train_combine.to_excel(F'{PARP_DIR}/feature_selection/data/datasets/{line_type}_combine_dataframe_train.xlsx', index=False)
    
    df_test_gen.to_excel(F'{PARP_DIR}/feature_selection/data/datasets/{line_type}_gen_dataframe_test.xlsx', index=False)
    df_test_immu.to_excel(F'{PARP_DIR}/feature_selection/data/datasets/{line_type}_immu_dataframe_test.xlsx', index=False)
    df_test_bio.to_excel(F'{PARP_DIR}/feature_selection/data/datasets/{line_type}_bio_dataframe_test.xlsx', index=False)
    df_test_cli.to_excel(F'{PARP_DIR}/feature_selection/data/datasets/{line_type}_clinical_dataframe_test.xlsx', index=False)
    df_test_combine.to_excel(F'{PARP_DIR}/feature_selection/data/datasets/{line_type}_combine_dataframe_test.xlsx', index=False)

def predict_score_for_combine_data(model, omic, cand_features, line_type):
    df_train_combine = pd.read_excel(io=F'{PARP_DIR}/feature_selection/data/datasets/{line_type}_combine_dataframe_train.xlsx', engine='openpyxl')
    df_test_combine = pd.read_excel(io=F'{PARP_DIR}/feature_selection/data/datasets/{line_type}_combine_dataframe_test.xlsx', engine='openpyxl')
    
    partial_hazard = model.predict_partial_hazard(df_train_combine[cand_features+['PFS(month)', 'PFS_type']])
    
    partial_hazard = partial_hazard.to_frame(name='partial_hazard')
    partial_hazard.to_excel(F'{PARP_DIR}/feature_selection/data/partial_hazard/{line_type}_{omic}_cox_partial_hazard_train_combine.xlsx', index=False)
    
    partial_hazard = model.predict_partial_hazard(df_test_combine[cand_features+['PFS(month)', 'PFS_type']])
    partial_hazard = partial_hazard.to_frame(name='partial_hazard')
    partial_hazard.to_excel(F'{PARP_DIR}/feature_selection/data/partial_hazard/{line_type}_{omic}_cox_partial_hazard_test_combine.xlsx', index=False)
    
def remove_missing_data(df):
    # for feature in df.columns:
    #     df = df[df[feature]!= -1]
    df_copy = copy.deepcopy(df)
    # if 'BRCA_HRD status' in df_copy.keys():
    #     df_copy = df_copy[df_copy['BRCA_HRD status']!= -1]
        
    # # 统计每一列中 -1 出现的比例
    # percentage_minus_one = df_copy.apply(lambda col: (col == -1).mean())
    # # 挑选出比例大于 0.5 的列名并转为 list
    # selected_columns = percentage_minus_one[percentage_minus_one > 0.8].index.tolist()
    

    # # 删掉空缺比例大于 0.5 的列
    # print('Origin Shape:', df_copy.shape)
    # print('Drop coloumns:', len(selected_columns))
    # df_copy = df_copy.drop(selected_columns, axis=1)
    # print('Final Shape:', df_copy.shape)
    # print("Process finished!")
    return df_copy

def stepwise_selection(data, df_test_data, duration_col, event_col, included, threshold_in=0.05, threshold_out=0.05):
    while True:
        changed = False
        model = CoxPHFitter().fit(data[included + [duration_col, event_col]], duration_col, event_col)
        pvalues = model.summary['p']
        print('Pvalues', pvalues)
        worst_pval = pvalues.max()
        if worst_pval > threshold_out:
            worst_feature = pvalues.idxmax()
            print('Remove feature:', worst_feature, worst_pval)
            included.remove(worst_feature)
            changed = True
        excluded = list(set(data.columns) - set(included) - {duration_col, event_col})
        new_pval = pd.Series(index=excluded)
        for new_column in excluded:
            print("new_column:", new_column)
            model = CoxPHFitter().fit(data[included + [new_column] + [duration_col, event_col]], duration_col, event_col, step_size=0.1)
            new_pval[new_column] = model.summary.loc[new_column, 'p']
        best_pval = new_pval.min()
        print('New pval in each step:', new_pval)
        if best_pval < threshold_in:
            best_feature = new_pval.idxmin()
            included.append(best_feature)
            print('Add in feature:', best_feature, best_pval)
            changed = True
        if not changed:
            break
    # Predict risk score.
    final_model = CoxPHFitter().fit(data[included + [duration_col, event_col]], duration_col, event_col)
    partial_hazard = final_model.predict_partial_hazard(data[included + [duration_col, event_col]])
    partial_hazard_test = final_model.predict_partial_hazard(df_test_data[included + [duration_col, event_col]])
    df_pvalues = pvalues.reset_index()
    df_pvalues.columns = ['feature', 'p']
    
    df_ex_pvalues = new_pval.reset_index()
    df_ex_pvalues.columns = ['feature', 'p']
    
    return final_model, df_pvalues, df_ex_pvalues, partial_hazard.to_frame(name='partial_hazard'), partial_hazard_test.to_frame(name='partial_hazard'), included

def coxFit(data, df_test_data, duration_col, event_col, included):
    # Predict risk score.
    final_model = CoxPHFitter().fit(data[included + [duration_col, event_col]], duration_col, event_col)
    partial_hazard = final_model.predict_partial_hazard(data[included + [duration_col, event_col]])
    partial_hazard_test = final_model.predict_partial_hazard(df_test_data[included + [duration_col, event_col]])
    pvalues = final_model.summary['p']
    df_pvalues = pvalues.reset_index()
    df_pvalues.columns = ['feature', 'p']
    return final_model, df_pvalues, partial_hazard.to_frame(name='partial_hazard'), partial_hazard_test.to_frame(name='partial_hazard')

def coxFit_for_each_modality(omic, line_type, cand_features):
    print(F'omic {omic}, line_type {line_type}, cand_features {cand_features}.')
    df_train = pd.read_excel(io=F'{PARP_DIR}/feature_selection/data/datasets/{line_type}_{omic}_dataframe_train.xlsx', engine='openpyxl')
    df_test = pd.read_excel(io=F'{PARP_DIR}/feature_selection/data/datasets/{line_type}_{omic}_dataframe_test.xlsx', engine='openpyxl')
    final_model, df_pvalues, partial_hazard_train, partial_hazard_test = coxFit(df_train, df_test, 'PFS(month)', 'PFS_type', cand_features)
    
    # df_pvalues.to_excel(F'{PARP_DIR}/feature_selection/data/p_value/{line_type}_{omic}_cox_pvalues_include.xlsx', index=False)
    # new_pval.to_excel(F'{PARP_DIR}/feature_selection/data/{line_type}_clinical_cox_pvalues_exclude.xlsx', index=False)
    partial_hazard_train.to_excel(F'{PARP_DIR}/feature_selection/data/partial_hazard/{line_type}_{omic}_cox_partial_hazard_train.xlsx', index=False)
    partial_hazard_test.to_excel(F'{PARP_DIR}/feature_selection/data/partial_hazard/{line_type}_{omic}_cox_partial_hazard_test.xlsx', index=False)
    final_model.summary.to_excel(F'{PARP_DIR}/feature_selection/data/p_value/{line_type}_{omic}_cox_summary_train.xlsx')
    
    predict_score_for_combine_data(final_model, omic, cand_features, line_type)
    
def select_clinical_features(line_type, p_thresh):
    candidate_features = pd.read_excel(io=F'{PARP_DIR}/feature_selection/data/feature_stability_{line_type}.xlsx', engine='openpyxl')
    cand_features = [f for f in candidate_features['features'] if f in CLINICAL_CHARA]
    
    # if 'Platinum sensitive' in cand_features:
    #     cand_features.remove('Platinum sensitive')
    df_basic = pd.read_excel(io=F'{PARP_DIR}/data/parp_stats_eng_{line_type}_800.xlsx', engine='openpyxl')
    
    df_basic_filter = remove_missing_data(df_basic[cand_features+['PFS(month)', 'PFS_type']])
    df_shuffled = df_basic_filter.sample(frac=1, random_state=18).reset_index(drop=True)
    train_size = int(len(df_shuffled) * 0.7)
    # 分割 DataFrame
    df_train = df_shuffled[:train_size]
    df_test = df_shuffled[train_size:]
    
    final_model, df_pvalues, new_pval, partial_hazard_train, partial_hazard_test, included_features = stepwise_selection(df_train, df_test, 'PFS(month)', 'PFS_type', cand_features[:1], threshold_in=p_thresh, threshold_out=p_thresh)
    df_pvalues.to_excel(F'{PARP_DIR}/feature_selection/data/p_value/{line_type}_clinical_cox_pvalues_selected.xlsx', index=False)
    new_pval.to_excel(F'{PARP_DIR}/feature_selection/data/p_value/{line_type}_clinical_cox_pvalues_exclude.xlsx', index=False)
    
    return included_features

def select_bio_features(line_type, p_thresh):
    candidate_features = pd.read_excel(io=F'{PARP_DIR}/feature_selection/data/feature_stability_{line_type}.xlsx', engine='openpyxl')
    cand_features = [f for f in candidate_features['features'] if f not in CLINICAL_CHARA+IMMU_CHARA+['OS(month)', 'OS_type']]
    
    df_basic = pd.read_excel(io=F'{PARP_DIR}/data/parp_stats_eng_{line_type}_800.xlsx', engine='openpyxl')
    
    df_basic_filter = remove_missing_data(df_basic[cand_features+['PFS(month)', 'PFS_type']])
    df_shuffled = df_basic_filter.sample(frac=1, random_state=18).reset_index(drop=True)
    train_size = int(len(df_shuffled) * 0.7)
    # 分割 DataFrame
    df_train = df_shuffled[:train_size]
    df_test = df_shuffled[train_size:]
    
    final_model, df_pvalues, new_pval, partial_hazard_train, partial_hazard_test, included_features = stepwise_selection(df_train, df_test, 'PFS(month)', 'PFS_type', cand_features[:1], threshold_in=p_thresh, threshold_out=p_thresh)
    df_pvalues.to_excel(F'{PARP_DIR}/feature_selection/data/p_value/{line_type}_bio_cox_pvalues_selected.xlsx', index=False)
    new_pval.to_excel(F'{PARP_DIR}/feature_selection/data/p_value/{line_type}_bio_cox_pvalues_exclude.xlsx', index=False)
    
    return included_features

def combine_features_data(line_type):
    df_train_combine = pd.read_excel(io=F'{PARP_DIR}/feature_selection/data/datasets/{line_type}_combine_dataframe_train.xlsx', engine='openpyxl')
    df_cli_ph_train = pd.read_excel(io=F'{PARP_DIR}/feature_selection/data/partial_hazard/{line_type}_clinical_cox_partial_hazard_train_combine.xlsx', engine='openpyxl')
    df_immu_ph_train = pd.read_excel(io=F'{PARP_DIR}/feature_selection/data/partial_hazard/{line_type}_immu_cox_partial_hazard_train_combine.xlsx', engine='openpyxl')
    df_bio_ph_train = pd.read_excel(io=F'{PARP_DIR}/feature_selection/data/partial_hazard/{line_type}_bio_cox_partial_hazard_train_combine.xlsx', engine='openpyxl')
    df_base_data_train = df_train_combine[['PFS(month)', 'PFS_type', 'BRCA_HRD status']]
    df_base_data_train['clinical_ph'] = list(df_cli_ph_train['partial_hazard'])
    df_base_data_train['immu_ph'] = list(df_immu_ph_train['partial_hazard'])
    df_base_data_train['bio_ph'] = list(df_bio_ph_train['partial_hazard'])
    
    df_test_combine = pd.read_excel(io=F'{PARP_DIR}/feature_selection/data/datasets/{line_type}_combine_dataframe_test.xlsx', engine='openpyxl')
    df_cli_ph_test = pd.read_excel(io=F'{PARP_DIR}/feature_selection/data/partial_hazard/{line_type}_clinical_cox_partial_hazard_test_combine.xlsx', engine='openpyxl')
    df_immu_ph_test = pd.read_excel(io=F'{PARP_DIR}/feature_selection/data/partial_hazard/{line_type}_immu_cox_partial_hazard_test_combine.xlsx', engine='openpyxl')
    df_bio_ph_test = pd.read_excel(io=F'{PARP_DIR}/feature_selection/data/partial_hazard/{line_type}_bio_cox_partial_hazard_test_combine.xlsx', engine='openpyxl')
    df_base_data_test = df_test_combine[['PFS(month)', 'PFS_type', 'BRCA_HRD status']]
    df_base_data_test['clinical_ph'] = list(df_cli_ph_test['partial_hazard'])
    df_base_data_test['immu_ph'] = list(df_immu_ph_test['partial_hazard'])
    df_base_data_test['bio_ph'] = list(df_bio_ph_test['partial_hazard'])
    return df_base_data_train, df_base_data_test
    
def combine_feature_model(df_data_train, df_data_test, line_type):
    combine_dict = {
        "combine_GC": ['BRCA_HRD status', 'clinical_ph'],
        "combine_GI": ['BRCA_HRD status', 'immu_ph'],
        "combine_GB": ['BRCA_HRD status', 'bio_ph'],
        "combine_CI": ['clinical_ph', 'immu_ph'],
        "combine_CB": ['clinical_ph', 'bio_ph'],
        "combine_IB": ['immu_ph', 'bio_ph'],
        "combine_GCI": ['BRCA_HRD status', 'clinical_ph', 'immu_ph'],
        "combine_GCB": ['BRCA_HRD status', 'clinical_ph', 'bio_ph'],
        "combine_GIB": ['BRCA_HRD status', 'immu_ph', 'bio_ph'],
        "combine_CIB": ['immu_ph', 'clinical_ph', 'bio_ph'],
        "combine_GCIB": ['BRCA_HRD status', 'clinical_ph', 'immu_ph', 'bio_ph'],
    }
    for combine, combine_fea in combine_dict.items():
        final_model = CoxPHFitter().fit(df_data_train[combine_fea + ['PFS(month)', 'PFS_type']], 'PFS(month)', 'PFS_type')
        partial_hazard = final_model.predict_partial_hazard(df_data_train[combine_fea + ['PFS(month)', 'PFS_type']])
        partial_hazard_test = final_model.predict_partial_hazard(df_data_test[combine_fea + ['PFS(month)', 'PFS_type']])
        partial_hazard = partial_hazard.to_frame(name='partial_hazard')
        partial_hazard_test = partial_hazard_test.to_frame(name='partial_hazard')
        partial_hazard.to_excel(F'{PARP_DIR}/feature_selection/data/partial_hazard/{line_type}_{combine}_cox_partial_hazard_train_combine.xlsx', index=False)
        partial_hazard_test.to_excel(F'{PARP_DIR}/feature_selection/data/partial_hazard/{line_type}_{combine}_cox_partial_hazard_test_combine.xlsx', index=False)
        final_model.summary.to_excel(F'{PARP_DIR}/feature_selection/data/p_value/{line_type}_{combine}_cox_summary_train_combine.xlsx')
    
if __name__ == '__main__':
    cli_features = [
        'Pathological type',
        'Primary surgery hospital',
        'Treatment regimen',
        'PARPi type',]
    immu_features = ['P53','ER','Ki67']
    bio_features = ['Lipid routine test - High-Density Lipoprotein Cholesterol', 'Lipid routine test - Triglycerides', 'Routine liver function tests-Total Bile Acids', 'Routine renal function program-Creatinine', 'Routine serological examination of hepatitis B virus-Hepatitis B Virus e Antigen', 'Tumor Markers (Protein microarray C-12)-Alpha-Fetoprotein', 'Blood routine + reticulocyte count - Absolute Basophil Count']
    # bio_features = [
    #     'Routine renal function program-Creatinine',
    #     'Human immunodeficiency virus antibody antigen combined assay -HIV antigen antibody',]
    for line_type in [
        'first_line',
        'post_line'
        ]:
        # 读取 Excel 数据
        df_basic = pd.read_excel(io=F'{PARP_DIR}/data/parp_stats_eng_{line_type}_800.xlsx', engine='openpyxl')
        p_thresh = 0.15 if line_type == 'first_line' else 0.15
        cli_features = select_clinical_features(line_type, p_thresh)
        # bio_features = select_bio_features(line_type, 0.3) # first version
        bio_features = select_bio_features(line_type, 0.05) # first version
        # bio_features = select_bio_features(line_type, p_thresh)
        print('bio_features:', bio_features)
        if 'Platinum sensitive' in cli_features:
            cli_features.remove('Platinum sensitive')
        # Splite data into training and testing set
        splite_dataframe(df_basic, 0.7, line_type, immu_cand=immu_features, bio_cand=bio_features, cli_cand=cli_features)
        for omics, cand_features in {
            'clinical': cli_features,
            'immu': immu_features,
            'bio': bio_features
        }.items():
            coxFit_for_each_modality(omics, line_type, cand_features)
        
        df_fl_data_train, df_fl_data_test = combine_features_data(line_type)
        combine_feature_model(df_fl_data_train, df_fl_data_test, line_type)
    
