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

- **Programming Language**: Python 3.8+
- **Key Dependencies**:
  - pandas
  - numpy
  - scikit-learn
  - lifelines
  - SHAP
  - seaborn
  - matplotlib
  - statsmodels
  - scipy

### Installation

```bash
# Clone the repository
git clone <repository_url>

# Navigate to the project directory
cd <project_directory>

# Create a virtual environment
python -m venv env
source env/bin/activate  # For Windows, use env\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

## Usage

### Data Preparation

Place the multi-omics data files in the `data/` directory. Supported formats include `.csv` and `.tsv`.

### Steps to Run

1. **Environment Setup**:

   ```bash
   python3 setup.sh
   ```

2. **Data Preprocessing**:

   ```bash
   python common/preprocess.py --input data/raw_data.csv --output data/processed_data.csv
   ```

3. **Correlation Analysis**:

   ```bash
   python feature_selection/correlation_analysis.py --data data/processed_data.csv --output data/correlation_results.csv
   ```

4. **Univariate Analysis**:

   ```bash
   python feature_selection/univariate_analysis.py --data data/processed_data.csv --output data/univariate_results.csv
   ```

5. **Single-modal Model Construction**:

   ```bash
   python model_construction/build_single_modal_model.py --config config/single_modal_config.yaml
   ```

6. **Stacking Model Construction**:

   ```bash
   python model_construction/train_stacking_model.py --config config/stacking_model_config.yaml
   ```

### Example Configuration Files

- `config/single_modal_config.yaml`

```yaml
model: cox_single_modal
parameters:
  input_modal: clinical
  output_dir: results/single_modal
```

- `config/stacking_model_config.yaml`

```yaml
model: cox_stacking
parameters:
  base_models:
    - random_forest
    - svm
    - bayesian
  stacking_method: weighted
output_dir: results/stacking_model
```

## File Structure

```
├── common
│   ├── preprocess.py         # Data preprocessing script
├── config
│   ├── single_modal_config.yaml # Single-modal model configuration file
│   ├── stacking_model_config.yaml # Stacking model configuration file
├── data
│   ├── raw_data.csv          # Raw data
│   ├── processed_data.csv    # Processed data
│   ├── correlation_results.csv # Correlation analysis results
│   ├── univariate_results.csv # Univariate analysis results
├── feature_selection
│   ├── correlation_analysis.py # Correlation analysis script
│   ├── univariate_analysis.py  # Univariate analysis script
├── model_construction
│   ├── build_single_modal_model.py # Single-modal model construction script
│   ├── train_stacking_model.py     # Stacking model training script
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

For any questions, please contact us at [email@example.com](mailto:email@example.com).

