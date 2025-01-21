# Project Description

## Background

This project is based on a real-world cohort of ovarian cancer patients in China. It aims to develop an innovative ensemble model by integrating multi-omics data, including genetic, clinical, biochemical, and pathological features, to improve risk stratification and prognosis prediction. The project focuses on:

1. **Feature Selection and Ranking**: Identifying key features that significantly impact the efficacy of PARP inhibitor maintenance therapy using multiple analytical methods.
2. **Model Innovation**: Enhancing risk stratification accuracy through ensemble modeling to support personalized treatment strategies.

This study explores the combination of Cox models with multimodal data stacking, presenting an efficient integration method that improves the accuracy of ovarian cancer prognosis predictions and facilitates better clinical decision-making.

## Main Features

1. **Multi-omics Data Integration**:

   - Integration of genetic, clinical, biochemical, and pathological data.
   - Support for both Cox univariate and multivariate analyses.

2. **Feature Selection**:

   - Preliminary feature selection using correlation analysis.
   - Identification of significant features through Cox univariate analysis.
   - Evaluation of multicollinearity issues using Variance Inflation Factor (VIF).
   - Stability analysis to ensure the reliability of selected features.

3. **Prediction Models**:

   - Development of a multimodal data stacking model based on Cox regression.
   - Introduction of bidirectional methods during single-modal model construction to enhance predictive capabilities.
   - Combination and stacking of outputs from different modalities to select the optimal predictive model.

## Requirements

### Prerequisites

- **Programming Language**: Python 3.6+
- **Key Dependencies**:
  - pandas
  - numpy
  - scikit-learn
  - lifelines
  - seaborn
  - matplotlib
  - statsmodels
  - scipy

### Installation

```bash
# Clone the repository
git clone https://github.com/xiongxa/ovarian_cancer_patient_stratification.git

# Navigate to the project directory
cd ovarian_cancer_patient_stratification

# Install dependencies
pip install -r requirements.txt
```

## Usage

### Data Preparation

Place the multi-omics data files in the `data/` directory. Supported formats include `.xlsx`. An example named 'parp_stats_800.xlsx' is located in data folder.

### Steps to Run

1. **Environment Setup**:

   ```bash
   python3 setup.sh
   ```

2. **Data Preprocessing**:

   ```bash
   python feature_selection/data_preprocess.py
   ```

3. **Correlation Analysis**:

   ```bash
   python feature_selection/correlation_analysis.py
   ```

4. **Univariate Analysis**:

   ```bash
   python feature_selection/uni_cox_analysis.py
   ```
5. **Feature selection**:

   ```bash
   python feature_selection/feature_selection_combine.py
   ```

6. **Single-modal Model Construction**:

   ```bash
   python model_construction/multi_cox_stepwise_analysis.py
   ```

7. **Stacking Model Construction**:

   ```bash
   python model_construction/calculate_km_curve_c_index.py
   ```

## File Structure

```
├── common
│   ├── common_modules.py         # Common modules for this project
├── config
│   ├── constant.py # constant variables for this project
│   ├── label_match.xlsx # label match file from Chinese to English
├── data
│   ├── raw_data.csv          # Raw data
├── feature_selection
│   ├──data
│   ├── figures
│   ├── correlation_analysis.py # Correlation analysis script
│   ├── data_preprocess.py  # Data pre-process
│   ├── data_statistic_for_article.py 
│   ├── feature_selection_combine.py  # Feature selection
│   ├── uni_cox_analysis.py  # Univariate analysis script
├── model_construction
│   ├── calculate_kendalltau_coeff.py # Calculate kendalltau for each modality 
│   ├── calculate_km_curve_c_index.py     # Stacking model training script
│   ├── multi_cox_stepwise_analysis.py     # Single-modal model construction script
├── README.md                 # Project documentation
├── requirements.txt          # Project dependencies
├── setup.sh                  # Environment setup script
```

## Notes

1. **Data Privacy**: Ensure that all data is anonymized and compliant with ethical guidelines.
2. **Experimental Validation**: All model predictions should be validated experimentally.
3. **Environment Compatibility**: The project has been tested on Linux and macOS. Windows users may need additional adjustments.

## Contribution Guidelines

We welcome community contributions! Please follow these steps to submit a Pull Request:

1. Fork this repository.
2. Create a new branch: `git checkout -b feature-branch`.
3. Commit your changes: `git commit -m 'Add new feature'`.
4. Push the branch: `git push origin feature-branch`.
5. Create a Pull Request.

## Contact

For any questions, please contact us at [xiongxian@hnca.org.cn](xiongxian@hnca.org.cn).

