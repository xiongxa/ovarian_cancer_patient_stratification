import pandas as pd
import os

# Getting environment variables
PARP_DIR = os.environ.get('PARP_DIR')

def drop_unuse_coloums(input_file):
    """Map chinese column name to english.

    Args:
        input_file (Excel): Patients file in xlsx format.
    """
    # Reading Excel Data
    df_basic_info_sheet = pd.read_excel(io=input_file, engine='openpyxl', index_col=False)
    # Drop columns which are unuseless for analysis.
    df_basic_info_sheet = df_basic_info_sheet.drop(
        ['序号', '患者', '最近一次就诊时间', '病案号', '出生年月日',
        '第一次确诊时间', '第一次治疗结束时间', '第一次确诊时间',
        '第一次治疗结束时间', '最近一次治疗评价（0=PR；1=NA；2=CR）',
        '是否二次手术（仅限于二线维持治疗；0=否；是=1）', '维持治疗开始时间',
        '维持用药结束时间', '再次复发时间（未复发以随访时间为终点）', '死亡时间',
        '随访终点时间', 'PFI类别（未复发=0；复发=1）'], axis=1)
    print(F"Patients with {df_basic_info_sheet.shape[0]} rows and {df_basic_info_sheet.shape[1]} columns.")
    return df_basic_info_sheet

def map_chinese_to_eng(df, output_file):
    """Map chinese column name to english.

    Args:
        df (Excel): Patients informations in dataframe format.
        output_file(Dataframe): Output as excel file.
    """
    #Data preprocessing, transforming Chinese into English
    name_match_dict = {
        '年龄': 'Age',
        'BRCA突变（0=﹣，1=﹢）': 'BRCA status',
        'HRD':'HRD',
        '总胎数': 'Number of pregnancies',
        '流产数': 'Number of abortions',
        '生产数': 'Number of births',
        'p53评分（p53＞30%为﹢，0=﹣；﹢=1）': 'P53',
        'ER（0=弱；1=中；2=强）': 'ER',
        'p16（评分：0=-；1=+）': 'P16',
        'Ki67（所占百分比,阴性=0）': 'Ki67',
        'PS评分': 'PS',
        'PFS类别（未复发=0；复发=1）':'PFS_type',
        '病理分型（Serous=0;Mucinous=1;Clear cell=2;Endometrioid=3）': 'Pathological type', 
        '分期': 'Stage',
        '转移分类（未转移=0；脏器转移=1；腹腔、子宫、肠转移=2）': 'Tumor metastasis type', 
        '身高(m）': 'Height',
        '体重（kg）':'Weight',
        'BMI': 'BMI',
        '初次手术医院（0=肿瘤医院；1=非肿瘤医院）': 'Primary surgery hospital',
        '副反应情况（0=无；1=有）': 'Side reaction',
        '内分泌慢性病史（0=无；1=有）': 'Chronic endocrine history',
        '心血管疾病史（0=无；1=有）': 'History of cardiovascular disease',
        '传染性疾病史（0=无；1=有）': 'History of infectious diseases',
        '其他肿瘤史（0=无；1=有）': 'History of other tumors',
        '结婚年龄': 'Marriage age',
        '一线/二线及后线(0=一线；1=二线；2=后线）': 'First/second line and posterior line',
        '第二次手术医院（1=本院）': 'Second Surgery Hospital',
        '是否铂敏感(铂敏感=0，铂耐药=1)': 'Platinum sensitive',
        'parp前治疗方案（TP=0;TP+联合=1,Others2）': 'Treatment regimen',
        '维持药物名称（奥拉帕利=0；尼拉帕利=1；奥拉帕利+尼拉帕利=2；氟唑帕利=3）': 'PARPi type',
        'PFS（month）': 'PFS(month)',
        'OS（month）': 'OS(month)',
        'OS类别': 'OS_type',
    }
    # Reading Excel Files
    label_file = F"{PARP_DIR}/config/label_match.xlsx" 
    df_label_file = pd.read_excel(label_file, engine='openpyxl')
    # Converts the data frame to a dictionary
    translation_dict = dict(zip(df_label_file.iloc[:, 0], df_label_file.iloc[:, 1]))
    name_match_dict.update(translation_dict)
    df_basic_info_sheet_after_select = df[name_match_dict.keys()]
    df_basic_info_sheet_after_select.rename(columns=name_match_dict, inplace=True)
    df_basic_info_sheet_after_select.to_excel(output_file, index=False)
    print(F"Patients information with keys {list(df_basic_info_sheet_after_select.keys())}")
    return df_basic_info_sheet_after_select

def abnormal_value_process(df):
    # Handling outliers
    df['年龄'].replace('-1', -1, inplace=True)
    df.fillna(-1, inplace=True)

def only_keep_ola_nla(df_data):
    df_data = df_data[(df_data['PARPi type']==0)|(df_data['PARPi type']==1)]
    print(F'After delete non-ola-nila patients:{df_data.shape}')
    return df_data

# The status of BRCA and HRD are merged
def combine_brca_hrd(df_info):
    brca_status = []
    # Iterate over each row and get the value of any column
    for index, row in df_info.iterrows():
        if row['BRCA status'] ==1 or row['HRD'] == 1:
            brca_status.append(1)
        elif row['BRCA status'] ==0 or row['HRD'] == 0:
            brca_status.append(0)
        else:
            brca_status.append(-1)
    df_info['BRCA_HRD status'] = brca_status


def delete_abnormal_column(df):
    percentage_minus_one = df.apply(lambda col: (col == -1).mean())
    selected_columns = percentage_minus_one[percentage_minus_one > 0.5].index.tolist()
    
    if 'BRCA status' in selected_columns:
        selected_columns.remove('BRCA status')
    if 'HRD' in selected_columns:
        selected_columns.remove('HRD')

    print('Origin Shape:', df.shape)
    print('Drop coloumns:', len(selected_columns))
    df_after_del = df.drop(selected_columns, axis=1)
    print('Final Shape:', df_after_del.shape)
    print("Process finished!")
    return df_after_del
        
def output_first_line_data(df, output_file):
    df_basic_info_first_line = df[((df['First/second line and posterior line']==0) & (df['PFS(month)']>6))]
    df_basic_info_first_line = delete_abnormal_column(df_basic_info_first_line)
    combine_brca_hrd(df_basic_info_first_line)
    df_basic_info_first_line.to_excel(output_file, index=False)

def output_post_line_data(df, output_file):
    df_basic_info_post_line = df[((df['First/second line and posterior line']!=0) & (df['PFS(month)']>6))]
    df_basic_info_post_line = delete_abnormal_column(df_basic_info_post_line)
    combine_brca_hrd(df_basic_info_post_line)
    df_basic_info_post_line.to_excel(output_file, index=False)

if __name__ == '__main__':
    input_file = F'{PARP_DIR}/data/parp_stats_800.xlsx'
    output_file = F'{PARP_DIR}/data/parp_stats_eng_800.xlsx'
    df_patients = drop_unuse_coloums(input_file)
    abnormal_value_process(df_patients)
    df_patients_eng = map_chinese_to_eng(df_patients, output_file)
    print('df_patients_eng:', df_patients_eng.shape)
    df_patients_eng = only_keep_ola_nla(df_patients_eng)
    # get_label_type(df_patients_eng, 'PFS(month)', 'PFS_status')
    output_first_line_data(df_patients_eng, F'{PARP_DIR}/data/parp_stats_eng_first_line_800.xlsx')
    output_post_line_data(df_patients_eng, F'{PARP_DIR}/data/parp_stats_eng_post_line_800.xlsx')
    