import pandas as pd
from scipy.stats import mannwhitneyu
from data_preprocess import combine_brca_hrd, get_label_type
from config.constants import CLINICAL_CHARA_FOR_STAT, IMMU_CHARA
import statistics
import os
import matplotlib.pyplot as plt

# 获取环境变量
PARP_DIR = os.environ.get('PARP_DIR')

def draw_pie(line_type, df, target_col, count_col):
    # 筛选出 col1 为 0 的数据，并统计 col2 中 0 和 1 的数量
    group0 = df[df[target_col] == 0][count_col].value_counts(normalize=True)

    # 筛选出 col1 为 1 的数据，并统计 col2 中 0 和 1 的数量
    group1 = df[df[target_col] == 1][count_col].value_counts(normalize=True)

    fig, ax = plt.subplots(figsize=(10, 10))
    # 绘制第一个饼图（col1=0）
    plt.figure(figsize=(10, 5))

    print(group0)
    plt.pie(group0, labels=[F'{count_col}=-1', F'{count_col}=0', F'{count_col}=1'], autopct='%1.1f%%', startangle=90)
    plt.title(F'{target_col}=0')
    fig.savefig(F'{PARP_DIR}/feature_selection/figures/km_curves/{line_type}_{target_col}_{count_col}_0pie_curve.png',dpi=400,format='png')

    print(group1)
    fig, ax = plt.subplots(figsize=(10, 10))
    # 绘制第二个饼图（col1=1）
    plt.figure(figsize=(10, 5))
    plt.pie(group1, labels=[F'{count_col}=-1', F'{count_col}=0', F'{count_col}=1'], autopct='%1.1f%%', startangle=90)
    plt.title(F'{target_col}=1')

    fig.savefig(F'{PARP_DIR}/feature_selection/figures/km_curves/{line_type}_{target_col}_{count_col}_1pie_curve.png',dpi=400,format='png')

if __name__=='__main__':
    df_patients_first_line = pd.read_excel(F'{PARP_DIR}/data/parp_stats_eng_first_line_800.xlsx', engine='openpyxl')
    df_patients_post_line = pd.read_excel(F'{PARP_DIR}/data/parp_stats_eng_post_line_800.xlsx', engine='openpyxl')
    
    
    draw_pie('first_line', df_patients_first_line, 'PARPi type', 'BRCA_HRD status')
