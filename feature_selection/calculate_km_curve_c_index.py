import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, statistics
import pandas as pd
from scipy.stats import mannwhitneyu
import os
import numpy as np
from lifelines.utils import concordance_index
import random
import matplotlib.colors as mcolors
import copy

# 获取环境变量
PARP_DIR = os.environ.get('PARP_DIR')

def draw_survival_curve(df_patients_first_line, df_patients_post_line, fig_name, group1_name, group2_name, figure_type='PFS'):
    print(df_patients_first_line.shape, df_patients_post_line.shape)
    # 创建一个示例数据集
    data_first_line = {
        'time': list(df_patients_first_line[F'{figure_type}(month)']),
        'event': list(df_patients_first_line[F'{figure_type}_type'])
    }
    data_post_line = {
        'time': list(df_patients_post_line[F'{figure_type}(month)']),
        'event': list(df_patients_post_line[F'{figure_type}_type'])
    }
    
    # 进行Log-rank检验
    results = statistics.logrank_test(
        data_first_line['time'],
        data_post_line['time'],
        event_observed_A=data_first_line['event'],
        event_observed_B=data_post_line['event']
    )
    # 输出结果
    p_value = results.p_value
    
    # 创建 Kaplan-Meier 估计器
    kmf_first_line = KaplanMeierFitter()
    kmf_post_line = KaplanMeierFitter()

    # 适配数据
    kmf_first_line.fit(durations=data_first_line['time'], event_observed=data_first_line['event'], label='G1')
    kmf_post_line.fit(durations=data_post_line['time'], event_observed=data_post_line['event'], label='G2')

    # 输出中位生存时间
    median_survival_time_1 = round(kmf_first_line.median_survival_time_, 2) if str(kmf_first_line.median_survival_time_) != 'inf' else 'NR'
    median_survival_time_2 = round(kmf_post_line.median_survival_time_, 2) if str(kmf_post_line.median_survival_time_) != 'inf' else 'NR'
    print(f"Median survival time for Group 1: {median_survival_time_1}")
    # print(f"Group 1: Confidence interval: {kmf_first_line.confidence_interval_}")
    
    print(f"Median survival time for Group 2: {median_survival_time_2}") 
    # print(f"Group 2: Confidence interval: {kmf_post_line.confidence_interval_}")   
    fig, ax = plt.subplots(figsize=(10, 10))
    # 画生存曲线
    ax = kmf_first_line.plot_survival_function(ci_show=False, show_censors=True, censor_styles={"marker": "+", "ms": 10, "mew": 2}, linewidth=3, color='blue')
    ax = kmf_post_line.plot_survival_function(ci_show=False, show_censors=True, censor_styles={"marker": "+", "ms": 10, "mew": 2}, linewidth=3, color='red')
    ax.set_ylim(-0.1, 1.1)
    
    from matplotlib.lines import Line2D
    # 自定义图例条目
    custom_legend = [
        Line2D([0], [0], color='blue', lw=3, marker='+', markersize=10, label=F'{group1_name}(m{figure_type}: {median_survival_time_1},n={len(df_patients_first_line)})'),
        Line2D([0], [0], color='red', lw=3, marker='+', markersize=10, label=F'{group2_name}(m{figure_type}: {median_survival_time_2},n={len(df_patients_post_line)})'),
    ]

    # 添加图例
    plt.legend(handles=custom_legend, fontsize=20)

    if round(p_value, 3) >= 0.001:
        plt.text(2, 0.1, F'P={round(p_value, 3)}', fontsize=20)
    else:
        plt.text(2, 0.1, F'P<0.001', fontsize=20)
    plt.xlabel('Time(Month)', fontdict={'fontsize': 20})
    plt.ylabel(F'Proportion({figure_type})', fontdict={'fontsize': 20})
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.show()
    fig.savefig(F'{PARP_DIR}/feature_selection/figures/km_curves/{fig_name}_{figure_type}_survive_curve.png',dpi=400,format='png')

def cal_percentile(data):
    import numpy as np
    data = list(data)
    # 计算第一四分位数（Q1）
    q1 = np.percentile(sorted(data), 25)

    # 计算第三四分位数（Q3）
    q3 = np.percentile(sorted(data), 75)

    print(F"mean:{np.mean(np.array(data))}; std:{np.std(np.array(data))}; Q1:{q1}; Q3:{q3}.")
    return [np.mean(np.array(data)), np.std(np.array(data)), q1, q3]

def combine_partial_hazard_data(line_type, omics):
    if 'combine' in omics:
        df_patients = pd.read_excel(io=F'{PARP_DIR}/feature_selection/data/datasets/{line_type}_combine_dataframe_train.xlsx', engine='openpyxl')
        df_omics_patients = pd.read_excel(io=F'{PARP_DIR}/feature_selection/data/partial_hazard/{line_type}_{omics}_cox_partial_hazard_train_combine.xlsx', engine='openpyxl')
    else:
        df_patients = pd.read_excel(io=F'{PARP_DIR}/feature_selection/data/datasets/{line_type}_{omics}_dataframe_train.xlsx', engine='openpyxl')
        df_omics_patients = pd.read_excel(F'{PARP_DIR}/feature_selection/data/partial_hazard/{line_type}_{omics}_cox_partial_hazard_train.xlsx', engine='openpyxl')
        
        
    df_omics_patients['PFS(month)'] = df_patients['PFS(month)']
    df_omics_patients['PFS_type'] = df_patients['PFS_type']
    df_omics_patients['OS(month)'] = df_patients['OS(month)']
    df_omics_patients['OS_type'] = df_patients['OS_type']
    print(line_type, omics, max(list(df_omics_patients['partial_hazard'])))
    
    if 'combine' in omics:
        df_patients_test = pd.read_excel(io=F'{PARP_DIR}/feature_selection/data/datasets/{line_type}_combine_dataframe_test.xlsx', engine='openpyxl')
        df_omics_patients_test = pd.read_excel(io=F'{PARP_DIR}/feature_selection/data/partial_hazard/{line_type}_{omics}_cox_partial_hazard_test_combine.xlsx', engine='openpyxl')
    else:
        df_patients_test = pd.read_excel(io=F'{PARP_DIR}/feature_selection/data/datasets/{line_type}_{omics}_dataframe_test.xlsx', engine='openpyxl')
        df_omics_patients_test = pd.read_excel(F'{PARP_DIR}/feature_selection/data/partial_hazard/{line_type}_{omics}_cox_partial_hazard_test.xlsx', engine='openpyxl')

    df_omics_patients_test['PFS(month)'] = df_patients_test['PFS(month)']
    df_omics_patients_test['PFS_type'] = df_patients_test['PFS_type']
    df_omics_patients_test['OS(month)'] = df_patients_test['OS(month)']
    df_omics_patients_test['OS_type'] = df_patients_test['OS_type']
    
    return df_omics_patients, df_omics_patients_test

def calculate_best_thresh(line_type, omics):
    df_omics_patients, df_omics_patients_test = combine_partial_hazard_data(line_type, omics)
    partial_hazard = list(df_omics_patients['partial_hazard'])
    min_p = 100
    best_thresh = None
    best_data_group1, best_data_group2 = None, None
    print('Max partial hazard:', min(partial_hazard), max(partial_hazard))
    for thresh in np.arange(round(min(partial_hazard),2)+0.1, round(max(partial_hazard),2)-0.1 if max(partial_hazard) < 50 else 50, 0.1):
        df_group1 = df_omics_patients[df_omics_patients['partial_hazard']<thresh]
        df_group2 = df_omics_patients[df_omics_patients['partial_hazard']>=thresh]
        
        # 进行Log-rank检验
        results = statistics.logrank_test(
            df_group1['PFS(month)'],
            df_group2['PFS(month)'],
            event_observed_A=df_group1['PFS_type'],
            event_observed_B=df_group2['PFS_type']
        )
        # 输出结果
        p_value = results.p_value
    
        print(F'p_value:{p_value}, thresh:{thresh}')
        if p_value < min_p and p_value != 0:
            df_group1_test = df_omics_patients_test[df_omics_patients_test['partial_hazard']<thresh]
            df_group2_test = df_omics_patients_test[df_omics_patients_test['partial_hazard']>=thresh]
            # print(df_omics_patients_test.shape, df_group1_test.shape, df_group2_test.shape)
            if len(df_group1_test) > len(df_omics_patients_test)/10 and len(df_group2_test) > len(df_omics_patients_test)/10:
                best_thresh = thresh
                min_p = p_value
                best_data_group1 = df_group1
                best_data_group2 = df_group2
    print('Best pars:', round(min(partial_hazard),2), round(max(partial_hazard),2), best_thresh, min_p)
    df_group1_test = df_omics_patients_test[df_omics_patients_test['partial_hazard']<best_thresh]
    df_group2_test = df_omics_patients_test[df_omics_patients_test['partial_hazard']>=best_thresh]
    return best_data_group1, best_data_group2, df_group1_test, df_group2_test

def seperate_2catogery_group(df_data, column):
    df_data = df_data[(df_data['PARPi type']==0)|(df_data['PARPi type']==1)]
    # df_data = df_data[(df_data['PARPi type']==1)]
    # df_data = df_data[(df_data['BRCA_HRD status']!=1)]
    df_group1 = df_data[df_data[column]==1]
    df_group2 = df_data[df_data[column]==0]
    return df_group1, df_group2

def calculate_c_index_for_BRCA_group(df_all_patients, line_type, omics):
    df_patients = df_all_patients[df_all_patients['BRCA_HRD status']!=-1]
    from scipy import stats
    mean_c_indexs = []
    c_index_err_l, c_index_err_u  = [], []
    for data in [df_patients]:
        # 进行100次leave-one-out bootstrapping
        n_iterations = 100
        c_indices = []

        for _ in range(n_iterations):
            # Leave-one-out bootstrap sampling
            sampled_data = data.sample(n=len(data), replace=True, random_state=random.randint(0, 10000))
            
            # 计算c-Index
            c_index = concordance_index(sampled_data['PFS(month)'], sampled_data['BRCA_HRD status'], sampled_data['PFS_type'])
            c_indices.append(c_index)

        # 计算95%置信区间
        c_indices = np.array(c_indices)
        mean_c_index = np.mean(c_indices)
        
        # 计算标准误差
        std_dev = np.std(c_indices, ddof=1)
        std_err = std_dev / np.sqrt(len(c_indices))

        # 计算置信区间
        ci = stats.t.interval(0.95, df=len(c_indices)-1, loc=mean_c_index, scale=std_err)

        ci_lower = ci[0]
        ci_upper = ci[1]

        # print(f"Mean c-Index: {mean_c_index}")
        # print(f"95% CI: ({ci_lower}, {ci_upper})")
        
        mean_c_indexs.append(mean_c_index)
        c_index_err_l.append(mean_c_index- ci_lower)
        c_index_err_u.append(ci_upper- mean_c_index)
    
    # plt.clf()
    fig, ax = plt.subplots(figsize=(8, 12))
    # 画图
    print(mean_c_indexs, c_index_err_l, c_index_err_u)
    x_labels = ["All"]
    x = np.arange(len(x_labels))  # 转换为数值型的 X 轴标签
    yerr = [c_index_err_l, c_index_err_u]
    # plt.errorbar(x_labels, mean_c_indexs, yerr=yerr, fmt='o', ecolor='g', capthick=2)
    plt.bar(x_labels, mean_c_indexs, yerr=yerr, color='b', alpha=0.3, align='center', capsize=10)
    plt.ylabel('c-Index', fontdict={'fontsize': 20})
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(min(mean_c_indexs)-max(c_index_err_l)-0.35, max(mean_c_indexs)+max(c_index_err_u)+0.01)
    plt.gca().yaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))
    plt.show()
    fig.savefig(F'{PARP_DIR}/feature_selection/figures/c_index/{line_type}_{omics}_c_index.png',dpi=400,format='png')
    return mean_c_indexs, c_index_err_l, c_index_err_u

def calculate_c_index(line_type, omics, df_patients_BRCA=None, df_patients_BRCA_test=None):
    from scipy import stats
    df_patients, df_patients_test = combine_partial_hazard_data(line_type, omics)
    mean_c_indexs = []
    c_index_err_l, c_index_err_u  = [], []
    for data in [df_patients, df_patients_test]:
        # 进行100次leave-one-out bootstrapping
        n_iterations = 100
        c_indices = []

        for _ in range(n_iterations):
            # Leave-one-out bootstrap sampling
            sampled_data = data.sample(n=len(data), replace=True, random_state=random.randint(0, 10000))
            
            # 计算c-Index
            c_index = concordance_index(sampled_data['PFS(month)'], -sampled_data['partial_hazard'], sampled_data['PFS_type'])
            c_indices.append(c_index)

        # 计算95%置信区间
        c_indices = np.array(c_indices)
        mean_c_index = np.mean(c_indices)
        
        # 计算标准误差
        std_dev = np.std(c_indices, ddof=1)
        std_err = std_dev / np.sqrt(len(c_indices))

        # 计算置信区间
        ci = stats.t.interval(0.95, df=len(c_indices)-1, loc=mean_c_index, scale=std_err)

        ci_lower = ci[0]
        ci_upper = ci[1]

        # print(f"Mean c-Index: {mean_c_index}")
        # print(f"95% CI: ({ci_lower}, {ci_upper})")
        
        mean_c_indexs.append(mean_c_index)
        c_index_err_l.append(mean_c_index- ci_lower)
        c_index_err_u.append(ci_upper- mean_c_index)
    
    # plt.clf()
    fig, ax = plt.subplots(figsize=(8, 12))
    # 画图
    print(mean_c_indexs, c_index_err_l, c_index_err_u)
    x_labels = ["Train", "Test"]
    x = np.arange(len(x_labels))  # 转换为数值型的 X 轴标签
    yerr = [c_index_err_l, c_index_err_u]
    # plt.errorbar(x_labels, mean_c_indexs, yerr=yerr, fmt='o', ecolor='g', capthick=2)
    plt.bar(x_labels, mean_c_indexs, yerr=yerr, color='b', alpha=0.3, align='center', capsize=10)
    plt.ylabel('c-Index', fontdict={'fontsize': 20})
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(min(mean_c_indexs)-max(c_index_err_l)-0.005, max(mean_c_indexs)+max(c_index_err_u)+0.01)
    plt.gca().yaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))
    plt.show()
    fig.savefig(F'{PARP_DIR}/feature_selection/figures/c_index/{line_type}_{omics}_c_index.png',dpi=400,format='png')
    return mean_c_indexs, c_index_err_l, c_index_err_u

# 将标签分割为每个字母一行
def split_label(label):
    return '\n'.join(label)

def draw_c_index(x_labels, mean_c_indexs, c_index_err_l, c_index_err_u, figure_name):
    import matplotlib.cm as cm  
    # 将数据按照 values 排序
    sorted_data = sorted(zip(mean_c_indexs, x_labels, c_index_err_l, c_index_err_u))
    mean_c_indexs, x_labels, c_index_err_l, c_index_err_u = zip(*sorted_data)
    print(x_labels)
    # plt.clf()
    fig, ax = plt.subplots(figsize=(7, 10))
    
    color_map = cm.get_cmap('tab20', len(x_labels))  # 'tab20' 是一个有 20 种颜色的 colormap
    colors = [color_map(i) for i in range(len(x_labels))]

    # 画图
    yerr = [c_index_err_l, c_index_err_u]
    # plt.errorbar(x_labels, mean_c_indexs, yerr=yerr, fmt='o', ecolor='g', capthick=2)
    plt.bar(x_labels, mean_c_indexs, yerr=yerr, color=colors, alpha=1, align='center', capsize=10)
    plt.ylabel('c-Index', fontdict={'fontsize': 20})
    # plt.xticks(rotation=90)
    # 分割后的标签
    split_labels = [split_label(label) for label in x_labels]
    plt.xticks(range(len(x_labels)), split_labels, fontsize=12)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.ylim(min(mean_c_indexs)-max(c_index_err_l)-0.005, max(mean_c_indexs)+max(c_index_err_u)+0.005)
    plt.gca().yaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))
    plt.show()
    fig.savefig(F'{PARP_DIR}/feature_selection/figures/c_index/{figure_name}_c_index.png',dpi=400,format='png')

if __name__ == '__main__':
    # Draw first line and post line survival curve.
    df_patients_first_line = pd.read_excel(F'{PARP_DIR}/data/parp_stats_eng_first_line_800.xlsx', engine='openpyxl')
    df_patients_post_line = pd.read_excel(F'{PARP_DIR}/data/parp_stats_eng_post_line_800.xlsx', engine='openpyxl')
    draw_survival_curve(df_patients_first_line, df_patients_post_line, 'line_type', 'Newly diagnosed', 'Recurrent', figure_type='PFS')
    draw_survival_curve(df_patients_first_line, df_patients_post_line, 'line_type', 'Newly diagnosed', 'Recurrent', figure_type='OS')
    
    # Draw survival curve for BRCA patients.
    df_group1, df_group2 = seperate_2catogery_group(df_patients_first_line, 'BRCA_HRD status')
    draw_survival_curve(df_group1, df_group2, 'first_line_BRCA', 'HRD', 'HRP', figure_type='PFS')
    draw_survival_curve(df_group1, df_group2, 'first_line_BRCA', 'HRD', 'HRP', figure_type='OS')
    
    df_group1, df_group2 = seperate_2catogery_group(df_patients_post_line, 'BRCA_HRD status')
    draw_survival_curve(df_group1, df_group2, 'post_line_BRCA', 'HRD', 'HRP', figure_type='PFS')
    draw_survival_curve(df_group1, df_group2, 'post_line_BRCA', 'HRD', 'HRP', figure_type='OS')
    
    # Draw survival curve for PARPi type patients.
    df_group1, df_group2 = seperate_2catogery_group(df_patients_first_line, 'PARPi type')
    draw_survival_curve(df_group1, df_group2, 'first_line_PARPi', 'Niraparib', 'Olaparib', figure_type='PFS')
    draw_survival_curve(df_group1, df_group2, 'first_line_PARPi', 'Niraparib', 'Olaparib', figure_type='OS')
    
    df_group1, df_group2 = seperate_2catogery_group(df_patients_post_line, 'PARPi type')
    draw_survival_curve(df_group1, df_group2, 'post_line_PARPi', 'Niraparib', 'Olaparib', figure_type='PFS')
    draw_survival_curve(df_group1, df_group2, 'post_line_PARPi', 'Niraparib', 'Olaparib', figure_type='OS')

    exit()
    fl_mean_c_indexs, fl_c_index_err_l, fl_c_index_err_u = calculate_c_index_for_BRCA_group(df_patients_first_line, 'first_line', 'BRCA')
    pl_mean_c_indexs, pl_c_index_err_l, pl_c_index_err_u = calculate_c_index_for_BRCA_group(df_patients_post_line, 'post_line', 'BRCA')
    
    # exit()
    c_index_label = ['G', 'C', 'P', 'B', 'GC', 'GP', 'GB', 'CP', 'CB', 'PB', 'GCP', 'GCB', 'GPB', 'CPB', 'GCPB']
    mean_c_indexs_fl, c_index_err_l_fl, c_index_err_u_fl = copy.deepcopy(fl_mean_c_indexs), copy.deepcopy(fl_c_index_err_l), copy.deepcopy(fl_c_index_err_u)
    train_mean_c_indexs_fl, train_c_index_err_l_fl, train_c_index_err_u_fl = copy.deepcopy(fl_mean_c_indexs), copy.deepcopy(fl_c_index_err_l), copy.deepcopy(fl_c_index_err_u)
    mean_c_indexs_pl, c_index_err_l_pl, c_index_err_u_pl = copy.deepcopy(pl_mean_c_indexs), copy.deepcopy(pl_c_index_err_l), copy.deepcopy(pl_c_index_err_u)
    train_mean_c_indexs_pl, train_c_index_err_l_pl, train_c_index_err_u_pl = copy.deepcopy(pl_mean_c_indexs), copy.deepcopy(pl_c_index_err_l), copy.deepcopy(pl_c_index_err_u)

    # Draw survival curve for diffrent line type and omics.
    for line_type in [
        'first_line', 
        'post_line'
                      ]:
        for omics in [
                    'clinical',
                      'immu',
                      'bio',
                      'combine_GC', 'combine_GI', 'combine_GB',
                      'combine_CI', 'combine_CB', 'combine_IB',
                      'combine_GCI', 'combine_GCB', 'combine_GIB', 'combine_CIB',
                      'combine_GCIB'
                      ]:
            print(line_type, omics)
            df_group1_train, df_group2_train, df_group1_test, df_group2_test = calculate_best_thresh(line_type, omics)
            mean_c_indexs, c_index_err_l, c_index_err_u = calculate_c_index(line_type, omics)
            draw_survival_curve(df_group1_train, df_group2_train, F'{line_type}_{omics}_train', 'Lower risk', 'Higher risk')
            draw_survival_curve(df_group1_train, df_group2_train, F'{line_type}_{omics}_train', 'Lower risk', 'Higher risk', figure_type='OS')
            draw_survival_curve(df_group1_test, df_group2_test, F'{line_type}_{omics}_test', 'Lower risk', 'Higher risk')
            draw_survival_curve(df_group1_test, df_group2_test, F'{line_type}_{omics}_test', 'Lower risk', 'Higher risk', figure_type='OS')
            if line_type == 'first_line':
                train_mean_c_indexs_fl.append(mean_c_indexs[0])
                train_c_index_err_l_fl.append(c_index_err_l[0])
                train_c_index_err_u_fl.append(c_index_err_u[0])
                
                mean_c_indexs_fl.append(mean_c_indexs[1])
                c_index_err_l_fl.append(c_index_err_l[1])
                c_index_err_u_fl.append(c_index_err_u[1])
            else:
                mean_c_indexs_pl.append(mean_c_indexs[1])
                c_index_err_l_pl.append(c_index_err_l[1])
                c_index_err_u_pl.append(c_index_err_u[1])
                
                train_mean_c_indexs_pl.append(mean_c_indexs[0])
                train_c_index_err_l_pl.append(c_index_err_l[0])
                train_c_index_err_u_pl.append(c_index_err_u[0])

    fl_data_df = pd.DataFrame(
        {
            'c_index_label': c_index_label,
            'mean_c_indexs_fl': mean_c_indexs_fl,
            'c_index_err_l_fl': c_index_err_l_fl,
            'c_index_err_u_fl': c_index_err_u_fl,
            'mean_c_indexs_pl': mean_c_indexs_pl,
            'c_index_err_l_pl': c_index_err_l_pl,
            'c_index_err_u_pl': c_index_err_u_pl,
            
            'train_mean_c_indexs_fl': train_mean_c_indexs_fl,
            'train_c_index_err_l_fl': train_c_index_err_l_fl,
            'train_c_index_err_u_fl': train_c_index_err_u_fl,
            'train_mean_c_indexs_pl': train_mean_c_indexs_pl,
            'train_c_index_err_l_pl': train_c_index_err_l_pl,
            'train_c_index_err_u_pl': train_c_index_err_u_pl,
         }
    )
    fl_data_df.to_excel(F'{PARP_DIR}/feature_selection/data/all_model_c_index.xlsx', index=False)
    draw_c_index(c_index_label, mean_c_indexs_fl, c_index_err_l_fl, c_index_err_u_fl, 'first_line_all')
    draw_c_index(c_index_label, mean_c_indexs_pl, c_index_err_l_pl, c_index_err_u_pl, 'post_line_all')