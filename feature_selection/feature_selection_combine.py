import pandas as pd
from statsmodels.stats.outliers_influence import variance_inflation_factor
from statsmodels.tools.tools import add_constant
import config.constants as CONSTANTS
import os

# 获取环境变量
PARP_DIR = os.environ.get('PARP_DIR')

CLI_THRESH_F, BIO_THRESH_F = 0.1, 0.1  # First line
CLI_THRESH_P, BIO_THRESH_P = 0.1, 0.1  # Post line

def calculate_VIF(df, output_file):
    # 添加常数列，因为VIF计算需要截距
    df_with_const = add_constant(df)

    # 计算VIF
    vif_data = pd.DataFrame()
    vif_data["Variable"] = df_with_const.columns
    vif_data["VIF"] = [variance_inflation_factor(df_with_const.values, i) for i in range(df_with_const.shape[1])]

    # 打印VIF数据
    print(vif_data)
    vif_data.to_excel(output_file, index=False)
    return vif_data


def feature_selection_first_line(df_first_line_feature):
    df_cox_result = None
    multi_feature = []
    exclude_features = ['BRCA_HRD status', 'BRCA status', 'HRD', 'First/second line and posterior line', 'Second Surgery Hospital', 'PFS(month)', 'PFS_status', 'OS(month)', 'OS_status']
    for _, row in df_first_line_feature.iterrows():
        if row['feature'] in exclude_features:
            continue
        if row['feature'] in CONSTANTS.CB_CHRAR and row['p'] < CLI_THRESH_F:
            multi_feature.append(row['feature'])
        elif row['feature'] not in CONSTANTS.CB_CHRAR and row['p'] < BIO_THRESH_F:
            multi_feature.append(row['feature'])
    return multi_feature

def feature_selection_post_line(df_post_line_feature):
    df_cox_result = None
    multi_feature = []
    exclude_features = ['BRCA_HRD status', 'BRCA status', 'HRD', 'First/second line and posterior line', 'Second Surgery Hospital', 'PFS(month)', 'PFS_status', 'OS(month)', 'OS_status']
    for _, row in df_post_line_feature.iterrows():
        if row['feature'] in exclude_features:
            continue
        if row['feature'] in CONSTANTS.CB_CHRAR and row['p'] < CLI_THRESH_P:
            multi_feature.append(row['feature'])
        elif row['feature'] not in CONSTANTS.CB_CHRAR and row['p'] < BIO_THRESH_P:
            multi_feature.append(row['feature'])
    return multi_feature

def feature_selection_via_VIF(df_VIF):
    df_after_VIF = df_VIF[df_VIF['VIF']<10]
    return list(df_after_VIF['Variable'])


def elastic_net_stability_ana(df, features, output_file):
    import pandas as pd
    from sklearn.linear_model import ElasticNetCV
    from sklearn.model_selection import train_test_split
    from sklearn.preprocessing import StandardScaler
    import numpy as np

    # 划分特征和目标变量
    if 'const' in features:
        features.remove('const')
    X = df[features]
    y = df['PFS(month)']
    
    print('features:', features)
    print(df)
    print(X)
    print(y)
    
    # 划分训练集和测试集
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # 数据标准化
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    # 创建Elastic Net模型，使用交叉验证选择alpha和l1_ratio
    elastic_net = ElasticNetCV(cv=5, random_state=42)
    elastic_net.fit(X_train_scaled, y_train)

    # 打印最佳的alpha和l1_ratio
    print("Best alpha:", elastic_net.alpha_)
    print("Best l1_ratio:", elastic_net.l1_ratio_)

    # 获取选择的特征
    selected_features = X.columns[elastic_net.coef_ != 0]

    # 打印选中的特征
    print("Selected features:", selected_features)

    # 进行稳定性分析
    # X_stab_test = df[selected_features]
    X_stab_test = df[features]
    num_iterations = 10
    selected_features_matrix = np.zeros((num_iterations, X_stab_test.shape[1]), dtype=bool)

    for i in range(num_iterations):
        X_train, X_test, y_train, y_test = train_test_split(X_stab_test, y, test_size=0.2, random_state=i)
        X_train_scaled = scaler.fit_transform(X_train)
        elastic_net.fit(X_train_scaled, y_train)
        selected_features_matrix[i, :] = elastic_net.coef_ != 0

    # 统计特征的稳定性
    feature_stability = np.mean(selected_features_matrix, axis=0)

    # 打印特征稳定性
    print("Feature stability:", feature_stability)
    # print('len:', len())
    df_selected_features = pd.DataFrame({
        'features': list(features),
        'feature stability': list(feature_stability),})
    print(df_selected_features)
    df_selected_features = df_selected_features[df_selected_features['feature stability']>=0.7]
    df_selected_features.to_excel(output_file, index=False)

    return selected_features, feature_stability

def combine_VIF_stability(df_VIF, feature_sel_VIF, feature_stability, output_file):
    print(df_VIF)
    df_VIF = df_VIF[df_VIF['Variable']!='const']
    df_sel_VIF_stab = pd.DataFrame({
        'Variable': feature_sel_VIF,
        'feature stability': feature_stability,
    })
    print(df_VIF)
    print(df_sel_VIF_stab)
    df_VIF_stab = pd.merge(df_VIF, df_sel_VIF_stab, on='Variable', how='left')
    df_VIF_stab.to_excel(output_file, index=False)
    
if __name__ == '__main__':
    # First line patients
    df_first_line_feature = pd.read_excel(io=F'{PARP_DIR}/feature_selection/data/feature_ana_first_linePFS.xlsx', engine='openpyxl')
    features_first_line = feature_selection_first_line(df_first_line_feature)
    print(F"Features in first line analysis {features_first_line}")
    df_patients_first_line = pd.read_excel(F'{PARP_DIR}/data/parp_stats_eng_first_line_800.xlsx', engine='openpyxl')
     # Save Patients baisc statistic information.
    df_patients_first_line[features_first_line].describe().to_excel(F'{PARP_DIR}/feature_selection/data/feature_stats_first_line.xlsx', index=False)
    # Calculate VIF for each variables.
    df_VIF_first_line = calculate_VIF(df_patients_first_line[features_first_line], F'{PARP_DIR}/feature_selection/data/feature_VIF_first_line.xlsx')
    # Stability analysis for each variables.
    features_VIF_first_line = feature_selection_via_VIF(df_VIF_first_line)
    print(features_VIF_first_line)
    _, feature_stability_fl = elastic_net_stability_ana(df_patients_first_line, features_VIF_first_line, F'{PARP_DIR}/feature_selection/data/feature_stability_first_line.xlsx')
    # Combine VIF and stability value
    combine_VIF_stability(df_VIF_first_line, features_VIF_first_line, feature_stability_fl, F'{PARP_DIR}/feature_selection/data/feature_VIF_stability_first_line.xlsx')
    
    
    # Post line patients
    df_post_line_feature = pd.read_excel(io=F'{PARP_DIR}/feature_selection/data/feature_ana_post_linePFS.xlsx', engine='openpyxl')
    features_post_line = feature_selection_post_line(df_post_line_feature)
    print(F"Features in post line analysis {features_post_line}")
    df_patients_post_line = pd.read_excel(F'{PARP_DIR}/data/parp_stats_eng_post_line_800.xlsx', engine='openpyxl')
    # Save Patients baisc statistic information.
    df_patients_post_line[features_post_line].describe().to_excel(F'{PARP_DIR}/feature_selection/data/feature_stats_post_line.xlsx', index=False)
    # Calculate VIF for each variables.
    df_VIF_post_line = calculate_VIF(df_patients_post_line[features_post_line], F'{PARP_DIR}/feature_selection/data/feature_VIF_post_line.xlsx')
    # Stability analysis for each variables.
    features_VIF_post_line = feature_selection_via_VIF(df_VIF_post_line)
    print(features_VIF_post_line)
    _, feature_stability_pl = elastic_net_stability_ana(df_patients_post_line, features_VIF_post_line, F'{PARP_DIR}/feature_selection/data/feature_stability_post_line.xlsx')
    # Combine VIF and stability value
    combine_VIF_stability(df_VIF_post_line, features_VIF_post_line, feature_stability_pl, F'{PARP_DIR}/feature_selection/data/feature_VIF_stability_post_line.xlsx')
