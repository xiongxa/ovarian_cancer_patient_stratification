from lifelines import CoxPHFitter
import lifelines
import pandas as pd
import os

# 获取环境变量
PARP_DIR = os.environ.get('PARP_DIR')

def evaluate_feature_cox(df, feat_col, target_type):
    """
    :param df: pandas DataFrame with each row being a patient and each column being a feature or outcome, containing columns "duration" (float) of survival time and "observed" (float) to delineate observed [1.0] or censored [0.0] outcomes 
    :param feat_col: column name (str) of feature values (float)
    :return: single-entry dict with {feat_col (str): [p (float), log(partial_hazard_ratio) (float)]}
    """

    col_list = [F'{target_type}(month)', F'{target_type}_type']
    feature_name = ''
    feature_result = {}
    if type(feat_col) is list:
        col_list += feat_col
        feature_name = '+'.join(feat_col)
    else:
        raise RuntimeError(F"Unknown feature type {feat_col}")

    model = CoxPHFitter(penalizer=0.0)
    try:
        model.fit(df[col_list], duration_col=F'{target_type}(month)', event_col=F'{target_type}_type')
    except (lifelines.exceptions.ConvergenceError, lifelines.exceptions.ConvergenceWarning) as e:
        try:
            model = CoxPHFitter(penalizer=0.2)
            model.fit(df[col_list], duration_col=F'{target_type}(month)', event_col=F'{target_type}_type')
        except (lifelines.exceptions.ConvergenceError, lifelines.exceptions.ConvergenceWarning) as e:
            return None

    for feature in feat_col:
        coef = model.summary.coef[feature]
        hr = model.summary['exp(coef)'][feature]
        p = model.summary.p[feature]
        feature_result[feature] = model.summary
        # print(feature_result)
    return feature_result

def feature_ana_first_line(df_basic_info_first_line, output_file):
    chra_list = df_basic_info_first_line.keys()
    print(df_basic_info_first_line.shape)
    for target_type in ['PFS', 'OS']:
        df_cox_result = None
        for feature in chra_list:
            if feature not in ['PFS(month)', 'PFS_type', 'OS(month)', 'OS_type','Height', 'Weight']:
                cox_result = pd.DataFrame()
                df_for_cox = df_basic_info_first_line[df_basic_info_first_line[feature]!= -1]
                print(F'Feature {feature}, After delete -1: {df_for_cox.shape}')
                cox_result = evaluate_feature_cox(df_for_cox, [feature], target_type)
                # cox_result = evaluate_feature_cox(df_basic_info_first_line, [feature])
                if cox_result is not None:
                    cox_result = cox_result[feature]
                    cox_result['feature'] = [feature]
                    if df_cox_result is None:
                        df_cox_result = cox_result
                    else:
                        df_cox_result = pd.concat([df_cox_result, cox_result], axis=0, ignore_index=True)
        df_cox_result.to_excel(output_file + target_type + '.xlsx', index=False)

def feature_ana_post_line(df_basic_info_post_line, output_file):
    chra_list = df_basic_info_post_line.keys()
    print(df_basic_info_post_line.shape)
    for target_type in ['PFS', 'OS']:
        df_cox_result = None
        for feature in chra_list:
            if feature not in ['PFS(month)', 'PFS_type', 'OS(month)', 'OS_type', 'Height', 'Weight']:
                cox_result = pd.DataFrame()
                df_for_cox = df_basic_info_post_line[df_basic_info_post_line[feature]!= -1]
                print(F'Feature {feature}, After delete -1: {df_for_cox.shape}')
                cox_result = evaluate_feature_cox(df_for_cox, [feature], target_type)
                # cox_result = evaluate_feature_cox(df_basic_info_post_line, [feature])
                if cox_result is not None:
                    cox_result = cox_result[feature]
                    cox_result['feature'] = [feature]
                    if df_cox_result is None:
                        df_cox_result = cox_result
                    else:
                        df_cox_result = pd.concat([df_cox_result, cox_result], axis=0, ignore_index=True)
        df_cox_result.to_excel(output_file + target_type + '.xlsx', index=False)

if __name__ == '__main__':
    # 读取 Excel 数据
    df_basic_first_line = pd.read_excel(io=F'{PARP_DIR}/data/parp_stats_eng_first_line_800.xlsx', engine='openpyxl', index_col=None)
    df_basic_post_line = pd.read_excel(io=F'{PARP_DIR}/data/parp_stats_eng_post_line_800.xlsx', engine='openpyxl', index_col=None)
    feature_ana_first_line(df_basic_first_line, F'{PARP_DIR}/feature_selection/data/feature_ana_first_line')
    feature_ana_post_line(df_basic_post_line, F'{PARP_DIR}/feature_selection/data/feature_ana_post_line')
    