from feature_selection.calculate_km_curve_c_index import combine_partial_hazard_data
import pandas as pd
import os
import numpy as np
from scipy.stats import kendalltau
import matplotlib.pyplot as plt
import seaborn as sns

# 获取环境变量
PARP_DIR = os.environ.get('PARP_DIR')

def get_partial_hazard_for_all_model():
    scores = {'first_line': {}, 'post_line': {}}
    for line_type in ['first_line', 'post_line']:
        df_patients_test = pd.read_excel(io=F'{PARP_DIR}/feature_selection/data/datasets/{line_type}_combine_dataframe_test.xlsx', engine='openpyxl')
        scores[line_type]['G'] = list(df_patients_test['BRCA_HRD status'])
        for omics in ['clinical', 'immu', 'bio']:
            df_omics_patients, df_omics_patients_test = combine_partial_hazard_data(line_type, omics)
            scores[line_type][omics] = list(df_omics_patients_test['partial_hazard'])
    return scores

def draw_kendall_heat_map(scores, line_type):
    # 计算两两之间的肯德尔等级相关系数
    models = list(scores.keys())
    models_rename = {'G': 'G', 'clinical': 'C', 'immu': 'P', 'bio': 'B'}
    models_label = [models_rename[name] for name in list(scores.keys())]
    n_models = len(models)
    tau_matrix = np.zeros((n_models, n_models))

    for i in range(n_models):
        for j in range(i, n_models):
            if i == j:
                tau_matrix[i, j] = 1.0  # 自相关为1
            else:
                tau, _ = kendalltau(scores[models[i]], scores[models[j]])
                tau_matrix[i, j] = abs(tau)
                tau_matrix[j, i] = abs(tau)

    # 创建 DataFrame 以便使用 Seaborn 绘图
    tau_df = pd.DataFrame(tau_matrix, index=models_label, columns=models_label)

    # 创建掩码来屏蔽上半部分和对角线
    mask = np.triu(np.ones_like(tau_matrix, dtype=bool), k=1)

    # 画出热力图
    # plt.clf()
    fig, ax = plt.subplots(figsize=(10, 10))
    sns.heatmap(tau_df, mask=mask, annot=True, annot_kws={"size": 20},cmap='coolwarm', vmin=0, vmax=1)
    # 修改坐标轴刻度字体大小
    plt.xticks(fontsize=20)  # 修改 x 轴标签字体大小
    plt.yticks(fontsize=20)  # 修改 y 轴标签字体大小
    
    # 图例字体大小
    cbar = plt.gca().collections[0].colorbar
    cbar.ax.tick_params(labelsize=20)
    fig.savefig(F'{PARP_DIR}/feature_selection/figures/{line_type}_kendall_heatmap.png',dpi=400,format='png')
    
scores = get_partial_hazard_for_all_model()
draw_kendall_heat_map(scores['first_line'], 'first_line')
draw_kendall_heat_map(scores['post_line'], 'post_line')