'''Calculate and save the final evaluation metrics.'''

# NOTE: Every evaluator will do this slightly differently depending on how the data is presented

import os
import sys
import json
import pandas as pd
import numpy as np
import itertools
from datetime import datetime, timezone
from scipy.stats import pearsonr
from scipy.spatial import distance
from sklearn.metrics import accuracy_score, classification_report

from config import EVALUATOR_NAME, EVALUATOR_INPUT_PATH

# Mapping from ground truth labels to DeepBICCN2 predictor output classes
# Based on: https://github.com/aertslab/CREsted-paper/blob/main/Figure_5/validated_enhancers_scoring.ipynb
DF_TO_CLASS_MAPPING = {
    'Astro': 'Astro',
    'Endo': 'Endo',
    'L2_3IT': 'L2_3IT',
    'L5IT': 'L5IT',
    'L5ET': 'L5ET',
    'L5_6NP': 'L5_6NP',
    'L6CT': 'L6CT',
    'L6IT': 'L6IT',
    'L6b': 'L6b',
    'Lamp5': 'Lamp5',
    'Micro_PVM': 'Micro_PVM',
    'OPC': 'OPC',
    'Oligo': 'Oligo',
    'Pvalb': 'Pvalb',
    'Sncg': 'Sncg',
    'Sst': 'Sst',
    'SstChodl': 'SstChodl',
    'VLMC': 'VLMC',
    'Vip': 'Vip',
    'L4IT': 'L5IT',
    'L6_IT_Car3': 'L6IT',
    'Lamp5_Lhx6': 'Lamp5',
    'Pvalb_Chc': 'Pvalb',
}

def calculate_and_save_metrics(saved_predictions_path, output_dir):
    try:
        if os.path.exists(saved_predictions_path):
            print("----- Starting Evaluation Calculation and Saving as CSV -----")
            MEASURED_DATA_PATH = EVALUATOR_INPUT_PATH # NOTE: This may not be the same for other evaluators
            print(f"Using measured data from: {MEASURED_DATA_PATH}")
            print(f"Using predictions from: {saved_predictions_path}")
            print(f"Evaluation results will be saved here: {output_dir}")
            
            seq_column = "enhancerID" # This can change depending on data

            evaluation_summary_filename = f"evaluation_summary_{EVALUATOR_NAME}.csv"
            evaluation_summary_filepath = os.path.join(output_dir, evaluation_summary_filename)
            
            # Initialize an empty list to get summary for all tasks
            all_task_evaluation_results = []
            
            try:
                # Load measured data file and predictions file ONCE (not with every function call).
                # NOTE: Evaluator builders: If measured_file_path is not a tab-separated file,
                # this line (pd.read_csv) will need to be adjusted or replaced with the
                # appropriate pandas read function (e.g., pd.read_excel, pd.read_csv with different sep)
                # or custom loading logic (e.g., for .npy files).
                enhancer_labels = pd.read_csv(MEASURED_DATA_PATH, sep='\t', header=0)
                print(enhancer_labels)

                # Now load predictions
                with open(saved_predictions_path, 'r') as f:
                    predictions_file_content = json.load(f)
                
                # Extract Predictor Name
                predictor_name_base = predictions_file_content.get("predictor_name", None) # Resort to None if predictor name is not available
                predictor_name = predictor_name_base.replace(" ", "_").replace("/", "_")
                
                if (
                    "prediction_tasks" not in predictions_file_content or
                    # Also flag cases in case prediction_tasks key is returned empty
                    not predictions_file_content["prediction_tasks"] or
                    # And flag if any 'predictions' keys are empty
                    any(not key.get("predictions") for key in predictions_file_content["prediction_tasks"])
                ):
                    print("WARNING: 'prediction_tasks' key missing, empty, or one of the tasks has empty predictions.")
                else:
                    dfs = []
                    for task in predictions_file_content["prediction_tasks"]:
                        cell_type = task["cell_type_requested"]
                        preds = task["predictions"]

                        # Convert predictions dict → DataFrame (id as index, one column for this cell type)
                        df = pd.DataFrame.from_dict(preds, orient="index", columns=[cell_type])
                        dfs.append(df)

                        # Merge all cell-type columns side by side
                        preds_df = pd.concat(dfs, axis=1)
                    print(preds_df)
                    prediction_task_data_nopredictions = [{k: v for k, v in task.items() if k != "predictions"} for task in predictions_file_content["prediction_tasks"]]
                    print(prediction_task_data_nopredictions)


                    task_evaluation_dict = calculate_enhancer_classAUROC(
                            measured_df=enhancer_labels,
                            prediction_data=preds_df,
                            label_mapping=DF_TO_CLASS_MAPPING,
                            # chromosome_column and chromosomes_to_filter_list can be added as arguments
                    )
                    if task_evaluation_dict:
                        # Extract all metrics (matches CREsted paper notebook metrics)
                        accuracy_value = task_evaluation_dict.get('accuracy')
                        precision_value = task_evaluation_dict.get('precision')
                        recall_value = task_evaluation_dict.get('recall')
                        f1_value = task_evaluation_dict.get('f1_score')

                        # Get UTC timestamp for predictor_name
                        timestamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S.%f")
                        description = "Enhancer Classification"

                        # Create base metric dict
                        base_dict = {
                            "Evaluator": EVALUATOR_NAME,
                            "Description": description,
                            "Predictor_name": predictor_name,
                            "Time_stamp": timestamp,
                        }

                        # Add all metrics to results (matches notebook output)
                        all_task_evaluation_results.append({
                            **base_dict,
                            'Metric': 'accuracy',
                            'Value': str(accuracy_value),
                        })
                        all_task_evaluation_results.append({
                            **base_dict,
                            'Metric': 'precision',
                            'Value': str(precision_value),
                        })
                        all_task_evaluation_results.append({
                            **base_dict,
                            'Metric': 'recall',
                            'Value': str(recall_value),
                        })
                        all_task_evaluation_results.append({
                            **base_dict,
                            'Metric': 'f1_score',
                            'Value': str(f1_value),
                        })

            except Exception as e:
                print(f"An error occurred during evaluation: {e}")
                
        # Once all the data is received, save them all into a summary CSV
            # print(all_task_correlation_results)
            if all_task_evaluation_results:
                summary_df = pd.DataFrame(all_task_evaluation_results)
                csv_file_exists: bool = os.path.isfile(evaluation_summary_filepath)
                try:
                    summary_df.to_csv(evaluation_summary_filepath, mode='a',
                                    sep='\t', header=(not csv_file_exists), index=False)
                    if csv_file_exists:
                        print("Appended to existing summary CSV file")
                    else:
                        print("Created a new summary CSV file")
                    print(f"Saved evaluation summary to {evaluation_summary_filepath}!")
                except IOError as e:
                    print("\nNo evaluation resuls were saved!")

        else:
            print("Evaluator run did not complete successfully.")
            print(f"Predictions file not found in '{saved_predictions_path}'.")
            print("Skipping evaluation calculation!")

    except Exception as e:
        print(f"An unexpected error occurred during evaluation calculations: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()

def calculate_enhancer_classAUROC(
    measured_df: pd.DataFrame,
    prediction_data: pd.DataFrame,
    label_mapping: dict = None,
    ):

    """
    Calculates Accuracy, Precision, Recall, and F1 scores for enhancer classification.
    Matches the metrics from the CREsted paper notebook.

    Args:
        measured_df (pd.DataFrame): The dataframe of Enhancer labels with 'target_ct' column.
        prediction_data (pd.DataFrame): Dataframe of enhancer sequences and their predictions for each of the cell types.
        label_mapping (dict, optional): Mapping from ground truth labels to predictor output classes.

    Returns:
        evaluation_details (dict): Dictionary containing accuracy, precision, recall, and f1_score (all weighted averages).
    """
    # Check for NAs in predictions
    print("Original size of the predictions is:")
    print(prediction_data.shape)
    na_rows = prediction_data[prediction_data.isna().any(axis=1)]
    if not na_rows.empty:
        print("Rows with NaN values:")
        print(na_rows)

    print(measured_df)

    # Apply label mapping to ground truth if provided
    if label_mapping:
        print("\n--- Applying label mapping to ground truth ---")
        print(f"Original unique labels: {sorted(measured_df['target_ct'].unique())}")
        measured_df = measured_df.copy()
        measured_df['target_ct'] = measured_df['target_ct'].map(label_mapping)
        print(f"Mapped unique labels: {sorted(measured_df['target_ct'].unique())}")

        # Check for any unmapped labels (will be NaN after mapping)
        unmapped = measured_df['target_ct'].isna().sum()
        if unmapped > 0:
            print(f"WARNING: {unmapped} labels could not be mapped!")

    # Get predicted cell type (argmax)
    prediction_data['max_cell_type'] = prediction_data.idxmax(axis=1)

    # Calculate accuracy
    acc = accuracy_score(measured_df['target_ct'], prediction_data['max_cell_type'])
    print("Accuracy:", acc)

    # Calculate classification report to get weighted precision, recall, and f1-score
    # This matches the CREsted paper notebook approach
    report = classification_report(
        measured_df['target_ct'],
        prediction_data['max_cell_type'],
        output_dict=True,
        zero_division=0
    )

    # Extract weighted averages (matches the notebook metrics)
    precision_weighted = report['weighted avg']['precision']
    recall_weighted = report['weighted avg']['recall']
    f1_weighted = report['weighted avg']['f1-score']

    print("Precision (weighted):", precision_weighted)
    print("Recall (weighted):", recall_weighted)
    print("F1-score (weighted):", f1_weighted)

    evaluation_details = {
        'accuracy': acc,
        'precision': precision_weighted,
        'recall': recall_weighted,
        'f1_score': f1_weighted
    }

    return evaluation_details