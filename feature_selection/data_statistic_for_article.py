import pandas as pd
from scipy.stats import mannwhitneyu
from data_preprocess import combine_brca_hrd, get_label_type
from config.constants import CLINICAL_CHARA_FOR_STAT, IMMU_CHARA
import statistics
import os

# 获取环境变量
PARP_DIR = os.environ.get('PARP_DIR')

def cal_fraction_for_one_feature(result_df, df, category_column):
    counted_cate_num = 0
    # 获取所有唯一的分类型变量值并去掉 -1
    unique_values = sorted([value for value in df[category_column].unique() if value != -1])
    
    feature_category_match = {
        'PS': {'0': '0 score', '1': '1 score', '2': '2 score', 'Missing': 'Missing'},
        'Pathological type': {'0': 'Serous', '1': 'Mucinous', '2': 'Clear cell', '3': 'Endometrioid', 'Missing': 'Missing'},
        'Stage': {'1': 'I', '2': 'II', '3': 'III', '4': 'IV', 'Missing': 'Missing'},
        'Tumor metastasis type': {'0': 'Without metastasis', '1': 'Organ metastasis', '2': 'Abdominal cavity, uterus and intestine metastasis', 'Missing': 'Missing'},
        'Primary surgery hospital': {'0': 'Cancer Hospital', '1': 'Non Cancer Hospital', 'Missing': 'Missing'},
        'Side reaction': {'0': 'No', '1': 'Yes', 'Missing': 'Missing'},
        'Chronic endocrine history': {'0': 'No', '1': 'Yes', 'Missing': 'Missing'},
        'History of cardiovascular disease': {'0': 'No', '1': 'Yes', 'Missing': 'Missing'},
        'History of infectious diseases': {'0': 'No', '1': 'Yes', 'Missing': 'Missing'},
        'History of other tumors': {'0': 'No', '1': 'Yes', 'Missing': 'Missing'},
        'Platinum sensitive': {'0': 'Yes', '1': 'No', 'Missing': 'Missing'},
        'Treatment before PARP': {'0': 'Cisplatin + Paclitaxel', '1': 'Cisplatin + Paclitaxel + Bevacizumab', '2': 'Others', 'Missing': 'Missing'},
        'PARP type': {'0': 'Olaparib', '1': 'Niraparib', '2': 'Olaparib + Niraparib', '3': 'Fluzoparib', 'Missing': 'Missing'},
        'First/second line and posterior line': {'0': 'First-line', '1': 'Second-line', '2': 'Post-line', 'Missing': 'Missing'},
        'Second Surgery Hospital': {'0': 'Others', '1': 'Hunan Cancer Hospital', 'Missing': 'Missing'},
        'BRCA_HRD status': {'0': 'No', '1': 'Yes', 'Missing': 'Missing'},
        'P53': {'0': 'No', '1': 'Yes', 'Missing': 'Missing'},
        'ER': {'0': 'Weak', '1': 'Medium', '2': 'Strong', 'Missing': 'Missing'},
        'P16': {'0': 'No', '1': 'Yes', 'Missing': 'Missing'},
    }
    continue_features = ['Age', 'Marriage age', 'Height', 'Weight', 'BMI', 'Ki67']
    half_cate_feature = {
        'Number of pregnancies': {'thresh': 3},
        'Number of births': {'thresh': 3},   
    }
    
    if category_column in continue_features:
        data_group1 = list(df[(df[category_column] != -1)][category_column])
        result_df = result_df.append({
        'Characteristics': category_column,
        'Count': str(round(statistics.mean(data_group1), 2)) + '±' + str(round(statistics.stdev(data_group1),2))
    }, ignore_index=True)
        # if category_column == 'Age':
        #     best_threshold = 60
        #     group1 = '<60'
        #     group2 = '>=60'
        # else:
        #     best_threshold = 0.5
        #     group1 = '<0.5'
        #     group2 = '>=0.5'
        # # 将当前分类型变量值的数据作为 Group 1
        # data_group1 = df[(df[category_column] < best_threshold) & (df[category_column] != -1)]
        # data_group2 = df[(df[category_column] >= best_threshold) & (df[category_column] != -1)]
        # # 将结果添加到 result_df 中
        # result_df = result_df.append({
        #     'Characteristics': str(group1),
        #     'Count': str(len(data_group1)) + '(' + str(round(len(data_group1)/float(len(df))*100,2)) + ')',
        # }, ignore_index=True)
        # result_df = result_df.append({
        #     'Characteristics': str(group2),
        #     'Count': str(len(data_group2)) + '(' + str(round(len(data_group2)/float(len(df))*100,2)) + ')',
        # }, ignore_index=True)
    elif category_column in half_cate_feature.keys():
        result_df = result_df.append({
        'Characteristics': category_column,
        'Count': ''
    }, ignore_index=True)
        thresh = half_cate_feature[category_column]['thresh']
        # 遍历每个唯一的分类型变量值
        for group1 in unique_values:
            if group1 < thresh:
                # 将当前分类型变量值的数据作为 Category
                data_group1 = df[(df[category_column] == group1) & (df[category_column] != -1)]
                N_group1 = len(data_group1)
                counted_cate_num += N_group1
                # 将结果添加到 result_df 中
                result_df = result_df.append({
                    'Characteristics': str('\quad ') + str(group1),
                    'Count': str(N_group1) + '(' + str(round(N_group1/float(len(df))*100,2)) + ')',
                }, ignore_index=True)
            else:
                data_group1 = df[(df[category_column] >= thresh) & (df[category_column] != -1)]
                N_group1 = len(data_group1)
                counted_cate_num += N_group1
                # 将结果添加到 result_df 中
                result_df = result_df.append({
                    'Characteristics':str('\quad ') + str('>=') + str(thresh),
                    'Count':str(N_group1) + '(' + str(round(N_group1/float(len(df))*100,2)) + ')',
                }, ignore_index=True)
                break
        if counted_cate_num != len(df):
            result_df = result_df.append({
                'Characteristics':str('\quad ') + 'Missing',
                'Count':str(len(df)-counted_cate_num) + '(' + str(round((len(df)-counted_cate_num)/float(len(df))*100,2)) + ')',
            }, ignore_index=True)
    else:
        result_df = result_df.append({
        'Characteristics': category_column,
        'Count': ''
    }, ignore_index=True)
        print(category_column)
        # 遍历每个唯一的分类型变量值
        for group1 in unique_values:
            # 将当前分类型变量值的数据作为 Category
            data_group1 = df[(df[category_column] == group1) & (df[category_column] != -1)]
            N_group1 = len(data_group1)
            counted_cate_num += N_group1
            # 将结果添加到 result_df 中
            result_df = result_df.append({
                'Characteristics':str('\quad ') + feature_category_match[category_column][str(group1)],
                'Count': str(N_group1) + '(' + str(round(N_group1/float(len(df))*100,2)) + ')',
            }, ignore_index=True)
        if counted_cate_num != len(df):
            result_df = result_df.append({
                'Characteristics':str('\quad ') + 'Missing',
                'Count':str(len(df)-counted_cate_num) + '(' + str(round((len(df)-counted_cate_num)/float(len(df))*100,2)) + ')',
            }, ignore_index=True)
    return result_df

def cal_fraction_for_features(df, features, output_file):
    # 创建一个空的 DataFrame 用于存储结果
    result_df = pd.DataFrame(columns=['Characteristics', 'Count'])
    for feature in features:
        df_feature_PFS = df[[feature]]
        result_df = cal_fraction_for_one_feature(result_df, df_feature_PFS, feature)
    result_df.to_excel(output_file, index=False)

if __name__=='__main__':
    input_file = F'{PARP_DIR}/data/parp_stats_eng_800.xlsx'
    df_patients_all = pd.read_excel(input_file, engine='openpyxl', index_col=None)
    print('Original:', df_patients_all.shape)
    df_patients_af_select = df_patients_all[
     ((df_patients_all['PFS(month)']>6) & ((df_patients_all['PARP type']==0)|(df_patients_all['PARP type']==1)))
        ]
    print('After delete:', df_patients_af_select.shape)
    combine_brca_hrd(df_patients_af_select)
    get_label_type(df_patients_af_select, 'PFS(month)', 'PFS_status')
    
    
    cal_fraction_for_features(df_patients_af_select, CLINICAL_CHARA_FOR_STAT, F'{PARP_DIR}/feature_selection/data/clinical_statistic.xlsx')
    cal_fraction_for_features(df_patients_af_select, IMMU_CHARA, F'{PARP_DIR}/feature_selection/data/immu_statistic.xlsx')