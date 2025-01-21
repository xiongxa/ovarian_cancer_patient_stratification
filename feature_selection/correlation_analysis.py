import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import config.constants as CONSTANTS
import os
from matplotlib.patches import Rectangle

# 获取环境变量
PARP_DIR = os.environ.get('PARP_DIR')


def plot_heatmap(df, line_type, modal):
    print(df)
    print(df.keys())
    # df = df.rename(columns={"Treatment before PARP":"Chemotherapy regimen",})
    if not list(df.keys()):
        return
    df_sorted = df[sorted(df.columns)]
    data_corr = df_sorted.corr(method='spearman')
    # if modal == 'bio':
    #     # 设置阈值
    #     threshold = 0.5
    #     # 根据阈值筛选列
    #     selected_columns = data_corr.columns[(data_corr.abs() >= threshold).any(axis=0)]
    #     # 重新构建相关性矩阵，只包括筛选后的列
    #     data_corr = data_corr[selected_columns].loc[selected_columns]
    plt.cla()
    # fig, _ = plt.subplots(figsize=(10, 10))
    # 使用Seaborn绘制热力图
    fig = plt.figure(figsize=(28, 28))
    heatmap = sns.heatmap(data_corr, annot=False, cmap='coolwarm', fmt='.2f', linewidths=.01)
    
    # 定义特征组
    # 每组包含 10 个特征，设置对应的颜色和图例标签
    groups_fl = {
        "Blood type analysis": {"range": (0, 7), "color": "red", "label_color": "red"},
        "Blood cell analysis": {"range": (8, 59), "color": "blue", "label_color": "blue"},
        "Blood clotting routine": {"range": (60, 65), "color": "green", "label_color": "green"},
        "Blood routine": {"range": (66, 93), "color": "purple", "label_color": "purple"},
        "Liver function routine": {"range": (113, 121), "color": "orange", "label_color": "orange"},
        "Renal function routine": {"range": (122, 125), "color": "brown", "label_color": "brown"},
        "Tumor Markers": {"range": (139, 148), "color": "magenta", "label_color": "magenta"},
        "Urine analysis": {"range": (149, 159), "color": "cyan", "label_color": "cyan"},
        # "Group 9": {"range": (80, 89), "color": "magenta", "label_color": "magenta"},
        # "Group 10": {"range": (90, 99), "color": "lime", "label_color": "lime"},
    }

    groups_pl = {
        "Blood type analysis": {"range": (0, 3), "color": "red", "label_color": "red"},
        "Blood clotting routine": {"range": (4, 8), "color": "green", "label_color": "green"},
        "Blood routine": {"range": (9, 58), "color": "purple", "label_color": "purple"},
        "Lipid routine": {"range": (68, 71), "color": "blue", "label_color": "blue"},
        "Liver function routine": {"range": (73, 81), "color": "orange", "label_color": "orange"},
        "Renal function routine": {"range": (82, 85), "color": "brown", "label_color": "brown"},
        "Stool routine": {"range": (86, 89), "color": "lime", "label_color": "lime"},
        "Tumor Markers": {"range": (90, 100), "color": "magenta", "label_color": "magenta"},
        "Urine analysis": {"range": (101, 111), "color": "cyan", "label_color": "cyan"},
    }

    groups = groups_fl if line_type == 'first_line' else groups_pl
    # 为每个组添加边框
    for group, properties in groups.items():
        start, end = properties["range"]
        color = properties["color"]
        # 绘制矩形框
        rect = Rectangle(
            (start, start),  # 左下角坐标
            end - start + 1,  # 宽度
            end - start + 1,  # 高度
            linewidth=4,
            edgecolor=color,  # 边框颜色
            facecolor="none"  # 无填充
        )
        heatmap.add_patch(rect)

    # 添加图例
    legend_handles = [
        Rectangle((0, 0), 1, 1, edgecolor=properties["label_color"], facecolor="none", label=group)
        for group, properties in groups.items()
    ]
    plt.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(1.2, 1), title="Feature Groups", fontsize=20)


    cbar = heatmap.collections[0].colorbar
    
    
    # 倾斜横坐标标签
    heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation=45, ha='right')
    if modal != 'bio':
        plt.xticks(fontsize=30)
        plt.yticks(fontsize=30)
        cbar.ax.tick_params(labelsize=30)
    else:
        cbar.ax.tick_params(labelsize=20)
    plt.tight_layout()
    plt.show()
    fig.savefig(F'{PARP_DIR}/feature_selection/figures/correlations/correlation_{line_type}_{modal}.png',dpi=400,format='png')

def delete_empty_column(df):
    # 统计每一列中 -1 出现的比例
    percentage_minus_one = df.apply(lambda col: (col == -1).mean())

    # 挑选出比例大于 0.5 的列名并转为 list
    selected_columns = percentage_minus_one[percentage_minus_one > 0.5].index.tolist()

    # 删掉空缺比例大于 0.5 的列
    print('origin:', df.shape)
    print('drop:', len(selected_columns))
    data_for_immu_after_del = df.drop(selected_columns, axis=1)
    print('final:', data_for_immu_after_del.shape)
    return data_for_immu_after_del
    
def corre_ana_by_modal(df, line_type):
    # 读取 Excel 文件
    label_file = F"{PARP_DIR}/config/label_match.xlsx"  # 替换为你的 Excel 文件路径
    df_label_file = pd.read_excel(label_file, engine='openpyxl')
    bio_feature = list(df_label_file['english'])
    clinical_feature = CONSTANTS.CLINICAL_CHARA
    immu_feature = CONSTANTS.IMMU_CHARA
    # bio_feature = CONSTANTS.BIO_CHARA
    all_features = list(df.keys())
    df_clinical = df[list(set(clinical_feature).intersection(set(all_features)))]
    df_immu = df[list(set(immu_feature).intersection(set(all_features)))]
    df_bio = df[list(set(bio_feature).intersection(set(all_features)))]
    df_bio = delete_empty_column(df_bio)
    # plot_heatmap(df_clinical, line_type, 'clinical')
    # plot_heatmap(df_immu, line_type, 'immu')
    plot_heatmap(df_bio, line_type, 'bio')
    
    
if __name__=='__main__':
    df_patients_first_line = pd.read_excel(F'{PARP_DIR}/data/parp_stats_eng_first_line_800.xlsx', engine='openpyxl')
    df_patients_post_line = pd.read_excel(F'{PARP_DIR}/data/parp_stats_eng_post_line_800.xlsx', engine='openpyxl')
    corre_ana_by_modal(df_patients_first_line, 'first_line')
    corre_ana_by_modal(df_patients_post_line, 'post_line')