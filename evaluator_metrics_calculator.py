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

from config import EVALUATOR_NAME, EVALUATOR_INPUT_PATH, CELL_TYPES, TASK_NAME_PREFIX

# Mapping from ground truth labels to DeepBICCN2 predictor output classes.
# Based on: https://github.com/aertslab/CREsted-paper/blob/main/Figure_5/validated_enhancers_scoring.ipynb
# NOTE: values are the canonical (full) cell-type names the Predictor now returns
# in cell_type_requested/actual, so the argmax comparison stays aligned.
DF_TO_CLASS_MAPPING = {
    'Astro': 'Astrocyte',
    'Endo': 'Endothelial',
    'L2_3IT': 'L2/3 IT',
    'L5IT': 'L5 IT',
    'L5ET': 'L5 ET',
    'L5_6NP': 'L5/6 NP',
    'L6CT': 'L6 CT',
    'L6IT': 'L6 IT',
    'L6b': 'L6b',
    'Lamp5': 'Lamp5',
    'Micro_PVM': 'Micro-PVM',
    'OPC': 'OPC',
    'Oligo': 'Oligodendrocyte',
    'Pvalb': 'Pvalb',
    'Sncg': 'Sncg',
    'Sst': 'Sst',
    'SstChodl': 'Sst Chodl',
    'VLMC': 'VLMC',
    'Vip': 'Vip',
    'L4IT': 'L5 IT',
    'L6_IT_Car3': 'L6 IT',
    'Lamp5_Lhx6': 'Lamp5',
    'Pvalb_Chc': 'Pvalb',
}

METRIC_NAMES = ['accuracy', 'precision', 'recall', 'f1_score']

def calculate_and_save_metrics(saved_predictions_path, output_dir, is_success_response=True):
    """
    Scores the saved Predictor response and appends one row per metric (accuracy, precision,
    recall, f1_score) to evaluation_summary_<EVALUATOR_NAME>.csv.

    The 4 rows are always written. Every metric is "NaN" when the run is disqualified or
    cannot be scored (GAME spec: HTTP error, per-task error, missing/extra IDs,
    non-finite values, ...), so a failed run never vanishes from the summary.
    """
    print("----- Starting Evaluation Calculation and Saving as CSV -----")
    print(f"Using measured data from: {EVALUATOR_INPUT_PATH}")
    print(f"Using predictions from: {saved_predictions_path}")
    print(f"Evaluation results will be saved here: {output_dir}")

    # Defaults written when anything below fails
    predictor_name = "UnknownPredictor"
    prediction_task_data_nopredictions = []
    metric_values = {metric: "NaN" for metric in METRIC_NAMES}

    try:
        with open(saved_predictions_path, 'r') as f:
            predictions_file_content = json.load(f)

        # Guard the predictor name so a missing key can't crash .replace()
        predictor_name_base = predictions_file_content.get("predictor_name") or "UnknownPredictor"
        predictor_name = predictor_name_base.replace(" ", "_").replace("/", "_")

        prediction_tasks = predictions_file_content.get("prediction_tasks") or []
        prediction_task_data_nopredictions = [
            {k: v for k, v in task.items() if k != "predictions"}
            for task in prediction_tasks if isinstance(task, dict)
        ]

        if not is_success_response:
            print("Predictor did not return a 200 OK status -> not scoring (metrics = NaN).")
        else:
            enhancer_labels = pd.read_csv(EVALUATOR_INPUT_PATH, sep='\t', header=0)
            scores = score_predictions(prediction_tasks, enhancer_labels)
            if scores is not None:
                metric_values = {metric: str(scores[metric]) for metric in METRIC_NAMES}

    except Exception as e:
        print(f"An error occurred during evaluation -> metrics = NaN: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()

    # Always write the 4 metric rows, even if scoring failed above
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S.%f")
    summary_df = pd.DataFrame([{
        "evaluator_name": EVALUATOR_NAME,
        "description": "Enhancer Classification",
        "predictor_name": predictor_name,
        "time_stamp": timestamp,
        "metric": metric,
        "value": metric_values[metric],
        "prediction_task(s)_data": prediction_task_data_nopredictions,
    } for metric in METRIC_NAMES])

    evaluation_summary_filepath = os.path.join(output_dir, f"evaluation_summary_{EVALUATOR_NAME}.csv")
    csv_file_exists = os.path.isfile(evaluation_summary_filepath)
    try:
        summary_df.to_csv(evaluation_summary_filepath, mode='a',
                          sep='\t', header=(not csv_file_exists), index=False)
        print("Appended to existing summary CSV file" if csv_file_exists else "Created a new summary CSV file")
        print(f"Saved evaluation summary to {evaluation_summary_filepath}!")
    except IOError as e:
        print(f"\nNo evaluation results were saved! {e}", file=sys.stderr)

def score_predictions(prediction_tasks, enhancer_labels):
    """
    Checks the Predictor's tasks and, if they are all valid, computes the 4 metrics.

    Returns:
        dict of metric -> value, or None if the run is disqualified (the reason is printed)
    """
    if not prediction_tasks or any(not isinstance(task, dict) or not isinstance(task.get("predictions"), dict)
                                   or not task["predictions"] for task in prediction_tasks):
        print("'prediction_tasks' is missing or empty, or a task's predictions are empty or not an object -> not scoring.")
        return None

    # Identify each task's cell type by its name, which the Predictor must echo back unchanged,
    # not by cell_type_requested (a Predictor could echo that string slightly differently).
    # The response must hold exactly one task per cell type that was requested.
    cell_type_by_task_name = {f"{TASK_NAME_PREFIX}{cell_type}": cell_type for cell_type in CELL_TYPES}
    returned_names = [task.get("name") for task in prediction_tasks]
    unknown_names = sorted({str(n) for n in returned_names if n not in cell_type_by_task_name})
    duplicate_names = sorted({n for n in returned_names if returned_names.count(n) > 1 and n in cell_type_by_task_name})
    missing_names = sorted(set(cell_type_by_task_name) - set(returned_names))
    if unknown_names or duplicate_names or missing_names:
        print("Returned tasks do not match the tasks that were sent -> not scoring:")
        print(f"  unknown task names: {unknown_names}")
        print(f"  duplicate task names: {duplicate_names}")
        print(f"  missing task names: {missing_names}")
        return None

    # Every enhancer is labelled by the argmax across all 19 cell types, so if any
    # task returned an error the comparison is incomplete: disqualify the whole run.
    error_tasks = [task for task in prediction_tasks if "error" in task["predictions"]]
    if error_tasks:
        print(f"{len(error_tasks)} of {len(prediction_tasks)} "
              f"task(s) returned an error -> not scoring (the argmax needs every cell type):")
        for task in error_tasks:
            print(f"  task '{task.get('name')}': {task['predictions']['error']}")
        return None

    # Every task must return exactly the enhancer IDs that were sent (GAME spec,
    # Evaluator requests note 6). The request is built from this same file, so
    # its enhancerIDs are the IDs that were sent.
    sent_ids = set(enhancer_labels["enhancerID"].astype(str))
    id_problems = find_prediction_id_problems(prediction_tasks, sent_ids)
    if id_problems:
        print("Predictions do not hold exactly the sequence IDs that were sent -> not scoring:")
        for problem in id_problems:
            print(f"  {problem}")
        return None

    # Every prediction must be a finite number (GAME spec, Evaluator requests note 6)
    value_problems = find_nonfinite_prediction_problems(prediction_tasks)
    if value_problems:
        print("Predictions contain values that are not finite numbers -> not scoring:")
        for problem in value_problems:
            print(f"  {problem}")
        return None

    # One column per cell type (looked up from the task name), indexed by enhancerID
    preds_df = pd.concat(
        [pd.DataFrame.from_dict(task["predictions"], orient="index",
                                columns=[cell_type_by_task_name[task["name"]]])
         for task in prediction_tasks],
        axis=1,
    )
    print(preds_df)

    return calculate_enhancer_classAUROC(
        measured_df=enhancer_labels,
        prediction_data=preds_df,
        label_mapping=DF_TO_CLASS_MAPPING,
    )

def find_prediction_id_problems(prediction_tasks, sent_ids):
    """
    Checks that each task's predictions hold exactly the sequence IDs that were sent:
    none missing, none extra. Tasks that returned an error object are skipped here.

    Returns:
        list of str: one message per task with mismatched IDs (empty if all match)
    """
    problems = []
    for task in prediction_tasks:
        preds = task.get("predictions", {})
        if "error" in preds:
            continue
        returned_ids = set(preds.keys())
        missing = sent_ids - returned_ids
        extra = returned_ids - sent_ids
        if missing or extra:
            problems.append(
                f"task '{task.get('name')}': {len(missing)} missing, {len(extra)} extra "
                f"(e.g. missing {sorted(missing)[:3]}, extra {sorted(extra)[:3]})"
            )
    return problems

def find_nonfinite_prediction_problems(prediction_tasks):
    """
    Checks that every point prediction is a finite number (not NaN, inf, null, a string
    or a list). Tasks that returned an error object are skipped here.

    Returns:
        list of str: one message per task with bad values (empty if all are finite numbers)
    """
    problems = []
    for task in prediction_tasks:
        preds = task.get("predictions", {})
        if "error" in preds:
            continue
        bad_ids = [
            seq_id for seq_id, value in preds.items()
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not np.isfinite(value)
        ]
        if bad_ids:
            examples = {seq_id: preds[seq_id] for seq_id in bad_ids[:3]}
            problems.append(f"task '{task.get('name')}': {len(bad_ids)} non-finite value(s), e.g. {examples}")
    return problems

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

    # Get predicted cell type (argmax across the cell-type columns)
    predicted_ct = prediction_data.idxmax(axis=1)

    # Pair each enhancer's true label with its prediction by enhancerID, not by row
    # position: the Predictor may return the sequences in a different order than the CSV.
    true_ct = measured_df.set_index('enhancerID')['target_ct']
    paired = pd.concat([true_ct.rename('true_ct'), predicted_ct.rename('predicted_ct')],
                       axis=1, join='inner')
    print(f"Paired {len(paired)} enhancers by enhancerID "
          f"({len(true_ct)} labelled, {len(predicted_ct)} predicted)")

    # Calculate accuracy
    acc = accuracy_score(paired['true_ct'], paired['predicted_ct'])
    print("Accuracy:", acc)

    # Calculate classification report to get weighted precision, recall, and f1-score
    # This matches the CREsted paper notebook approach
    report = classification_report(
        paired['true_ct'],
        paired['predicted_ct'],
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