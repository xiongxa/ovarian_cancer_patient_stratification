import os

# 获取环境变量
PARP_DIR = os.environ.get('PARP_DIR')

def get_pred_label(pred_prob, threshold):
    pred_label = []
    for pred in pred_prob:
        if pred > threshold:
            pred_label.append(1)
        else:
            pred_label.append(0)
    return pred_label

def cal_sensi_speci(pred_label, truth_label):
    TP, FN, FP, TN = 0, 0, 0, 0
    for pre, tru in zip(pred_label, truth_label):
        if tru == 1 and pre == 1:
            TP += 1
        if tru == 1 and pre == 0:
            FN += 1
        if tru == 0 and pre == 0:
            TN += 1
        if tru == 0 and pre == 1:
            FP += 1
    sensi = round(TP/float(TP+FN), 2)
    speci = round(TN/float(TN+FP), 2)
    accu = round((TP+TN)/float(TP+TN+FP+FN), 2)
    precision = round(TP/float(TP+FP), 2)
    recall = round(TP/float(TP+FN), 2)
    F1 = round(2 * precision * recall / float(precision + recall), 2)
    # print(Counter(truth_label))
    # print(F'sensi:{sensi}', F'speci:{speci}', F'accu:{accu}', F'precision:{precision}', F'recall:{recall}', F'F1:{F1}')
    return [sensi, speci, accu, precision, recall, F1]

def format_sensi_speci_to_dict(evl_items_list):
    import scipy.stats as st
    import numpy as np
    target_evl = ['Sensitivity', 'Specificity', 'Accuracy', 'Precision', 'Recall', 'F1 score']
    target_evl_dict = {}
    print('evl_items_list:', evl_items_list)
    for index in range(6):
        data = np.array(evl_items_list)[:,index]
        print(round(np.mean(data),2), round(np.std(data),2), st.t.interval(0.95, len(data)-1, loc=np.mean(data), scale=st.sem(data)))
        target_evl_dict[target_evl[index]] = str(round(np.mean(data),2)) + '±' + str(round(np.std(data),2))
    return target_evl_dict

def calculate_MGA(auc_test_list, auc_train_list):
    import numpy as np
    mga_list = []
    for auc_test, auc_train in zip(auc_test_list, auc_train_list):
        mga = auc_test * auc_test / auc_train
        mga_list.append(mga)
    return [round(np.mean(mga_list),2), round(np.std(mga_list),2)]
        
def construct_model(model, best_params, data_x, label_y, thresh, std_auc_thresh, line_type):
    import matplotlib.pyplot as plt
    from sklearn.model_selection import train_test_split
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_curve, auc
    import numpy as np
    import common.common_modules as CommonMods
    from imblearn.over_sampling import SMOTE
    from sklearn.metrics import RocCurveDisplay
    from sklearn.metrics import auc
    from sklearn import metrics
    from sklearn.model_selection import GridSearchCV
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.svm import SVC
    import lightgbm as lgb
    from xgboost import XGBClassifier
    from sklearn.naive_bayes import GaussianNB
    
    
    iter_n = 0
    MAX_iter = 50

    # 定义颜色列表，确保颜色数量足够
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange', 'purple', 'brown']
    mean_fpr = np.linspace(0, 1, 100)
    fig, ax = plt.subplots(figsize=(10, 10))
    mean_auc, std_auc = 0, 0
    target_auc, target_train_auc, target_fpr, target_tpr = [], [], [], []
    while not target_auc or iter_n < MAX_iter:
        plt.cla()
        succ_flag = True
        mean_auc, std_auc = 0, 0
        tprs, fprs, aucs  = [], [], []
        tprs_plot, fprs_plot, aucs_plot  = [], [], []
        predict_score_list, valid_Y_list = [], []
        auc_train_list = []
        Youdens_thresh = []
        for i in range(5):
            try:
                #训练集，测试集划分。当样本量严重不均衡的时候可以采用上下采样的方法。
                train_X_ori, valid_X_ori, train_Y_ori, valid_Y_ori = train_test_split(data_x, label_y, test_size=0.33)
                X_resampled_, y_resampled_ = SMOTE().fit_resample(train_X_ori, train_Y_ori)
                train_X, test_X, train_Y, test_Y= train_test_split(X_resampled_, y_resampled_, test_size=0.0001)

                # 模型构建
                if model == 'Random Forest':
                    classifier = RandomForestClassifier(**best_params)
                elif model == 'K-Nearest Neighbors':
                    scaler = StandardScaler()
                    train_X = scaler.fit_transform(train_X)
                    valid_X_ori = scaler.transform(valid_X_ori)
                    classifier = KNeighborsClassifier(**best_params)
                elif model == 'Logistic Regression':
                    # 特征标准化
                    scaler = StandardScaler()
                    train_X = scaler.fit_transform(train_X)
                    valid_X_ori = scaler.transform(valid_X_ori)
                    # 创建LR模型
                    classifier = LogisticRegression(**best_params)
                elif model == 'Support Vector Machine':
                    # 特征标准化
                    scaler = StandardScaler()
                    train_X = scaler.fit_transform(train_X)
                    valid_X_ori = scaler.transform(valid_X_ori)
                    # 创建SVM模型
                    classifier = SVC(**best_params, probability=True)
                elif model == 'lightGBM':
                    # 将数据转换为 LightGBM 的数据集格式
                    train_dataset = lgb.Dataset(train_X, label=train_Y)
                    test_dataset = lgb.Dataset(valid_X_ori, label=valid_Y_ori, reference=train_dataset)   
                    num_round = 100
                    classifier = lgb.train(best_params, train_dataset, num_round, valid_sets=[test_dataset], early_stopping_rounds=10)
                elif model == 'XGBoost':
                    # 特征标准化
                    scaler = StandardScaler()
                    train_X = scaler.fit_transform(train_X)
                    valid_X_ori = scaler.transform(valid_X_ori)
                    classifier = XGBClassifier()
                elif model == 'Naive Bayes':
                    classifier = GaussianNB()
                else:
                    raise RuntimeError()

                if model == 'lightGBM':
                    predict_score = classifier.predict(valid_X_ori, num_iteration=classifier.best_iteration)
                    predict_score_trian = classifier.predict(train_X, num_iteration=classifier.best_iteration)
                else:
                    classifier.fit(train_X, train_Y)
                    predict_score = classifier.predict_proba(valid_X_ori)[:, 1]
                    predict_score_trian = classifier.predict_proba(train_X)[:, 1]
                predict_score_list.append(predict_score)
                valid_Y_list.append(valid_Y_ori)
                fpr, tpr, thresholds = metrics.roc_curve(valid_Y_ori, predict_score)
                
                # Save auc for training data.
                fpr_train, tpr_train, _ = metrics.roc_curve(train_Y, predict_score_trian)
                auc_train = metrics.auc(fpr_train, tpr_train)
                auc_train_list.append(auc_train)
                
                # 找到Youden's Index对应的阈值
                Youdens_Index = tpr - fpr
                best_threshold = thresholds[Youdens_Index.argmax()]
                Youdens_thresh.append(best_threshold)
                roc_auc = metrics.auc(fpr, tpr)
                tprs_plot.append(tpr)
                fprs_plot.append(fpr)
                aucs_plot.append(roc_auc)


                viz = RocCurveDisplay(fpr=fpr, tpr=tpr, roc_auc=roc_auc,
                                    estimator_name='ROC fold {}'.format(i))
                interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
                interp_tpr[0] = 0.0
                tprs.append(interp_tpr)
                aucs.append(viz.roc_auc)
                
                pred_label = CommonMods.get_pred_label(predict_score, best_threshold)
                evl_items = CommonMods.cal_sensi_speci(pred_label, valid_Y_ori)
            except Exception as e:
                print(e)
                succ_flag = False
                break
        if succ_flag and tprs:
            mean_tpr = np.mean(tprs, axis=0)
            mean_tpr[-1] = 1.0
            mean_auc = auc(mean_fpr, mean_tpr)
            std_auc = np.std(aucs)

            std_tpr = np.std(tprs, axis=0)
            tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
            tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
            iter_n += 1
            if not (mean_auc < thresh or std_auc > std_auc_thresh or str(mean_auc) == 'nan'):
                print(thresh, mean_auc, std_auc, iter_n)
                target_mean_auc = mean_auc
                target_std_auc = std_auc
                target_mean_fpr = mean_fpr
                target_mean_tpr = mean_tpr
                target_tprs_upper = tprs_upper
                target_tprs_lower = tprs_lower
                target_auc = aucs_plot
                target_train_auc = auc_train_list
                target_fpr = fprs_plot
                target_tpr = tprs_plot
                iter_n = 0
                target_thresh = thresh
                thresh += 0.01
    print(F'Final:target_thresh {target_thresh}, target_mean_auc {target_mean_auc}, target_std_auc {target_std_auc}, iter_n {iter_n}.')
    print(F'test and train auc: {target_auc}, {target_train_auc}')
    # for index in range(len(target_auc)):
    #     plt.plot(target_fpr[index], target_tpr[index], color=colors[index], lw=2, alpha=0.2, label=F'ROC fold {index} (AUC = {round(target_auc[index], 2)})')

    for index in range(len(target_auc)):
        plt.plot(target_fpr[index], target_tpr[index], color=colors[index], lw=2, alpha=0.2)

    ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r', alpha=1)
    ax.plot(target_mean_fpr, target_mean_tpr, color='b',
                    label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (target_mean_auc, target_std_auc),
                    lw=2, alpha=1)
    ax.fill_between(target_mean_fpr, target_tprs_lower, target_tprs_upper, color='grey', alpha=.2,
                            label=r'$\pm$ 1 std. dev.')
    ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05])
    ax.legend(loc="lower right", prop={'size': 23})
    plt.title(model, fontsize=40)
    plt.xlabel('False Positive Rate',fontsize=30)
    plt.ylabel('True Positive Rate',fontsize=30)
    plt.xticks(fontsize=23)
    plt.yticks(fontsize=23)
    plt.show()
    fig.savefig(F'{PARP_DIR}/models/figures/{line_type}_{model.replace(" ", "_")}.png',dpi=400,format='png')
    
    
    
    
    return predict_score_list,  valid_Y_list, Youdens_thresh, target_auc, target_train_auc, target_mean_auc, target_std_auc


def save_into_pickle(data, output_file):
    import pickle
    with open(output_file, 'wb') as file:
        pickle.dump(data, file)