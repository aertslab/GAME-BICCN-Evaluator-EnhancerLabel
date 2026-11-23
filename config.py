'''Configuration Script for Evaluator Name, Input File, and Preferred Data Format'''

import os
import json
from datetime import datetime

# --- Core Evaluator Settings ---
# Base evaluator name (the container build date is appended automatically below)
evaluator_base = "DeepBICCN2_EnhancerLabel_Evaluator"
# Name of the file being used for predictions
input_file = 'biccn_enhancers_withSequence.csv'

# --- Prediction tasks ---
# Canonical (full) cell-type names, matching what the DeepBICCN2 Predictor advertises in /help.
# One prediction task is sent per cell type, named f"{TASK_NAME_PREFIX}{cell_type}"; the metrics
# step uses that name to know which cell type each returned task belongs to.
CELL_TYPES = ["Astrocyte", "Endothelial", "L2/3 IT", "L5 ET", "L5 IT", "L5/6 NP", "L6 CT", "L6 IT", "L6b", "Lamp5", "Micro-PVM", "OPC", "Oligodendrocyte", "Pvalb", "Sncg", "Sst", "Sst Chodl", "VLMC", "Vip"]
TASK_NAME_PREFIX = "accessibility_"

# --- Directory Settings ---
# Get the absolute path of the script's directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Determine if running inside a container or not
if os.path.exists("/.singularity.d"):
    # Running inside the container
    EVALUATOR_DATA_DIR = "/evaluator_data"

    # Read the build timestamp from Apptainer's auto-generated labels file
    # to version the Evaluator name (Apptainer writes this label at build time).
    try:
        with open("/.singularity.d/labels.json", "r") as f:
            labels = json.load(f)
        raw_build_date = labels.get("org.label-schema.build-date", "")

        # Example format: "Friday_28_November_2025_18:6:29_PST"
        # Strip day-of-week and timezone, keep the core date+time
        parts = raw_build_date.split("_")
        # parts: ['Friday', '28', 'November', '2025', '18:6:29', 'PST']
        date_str = f"{parts[1]}_{parts[2]}_{parts[3]}_{parts[4]}"

        dt = datetime.strptime(date_str, "%d_%B_%Y_%H:%M:%S")
        build_timestamp = dt.strftime("%Y%m%d-%H%M%S")
        timezone_label = parts[5] if len(parts) > 5 else "UNK"
        EVALUATOR_NAME = f"{evaluator_base}_{build_timestamp}_{timezone_label}"

    except Exception as e:
        print(f"Warning: Could not parse build timestamp from labels.json: {e}")
        EVALUATOR_NAME = f"{evaluator_base}_unknown"
else:
    # Running outside the container (development mode)
    EVALUATOR_DATA_DIR = os.path.join(SCRIPT_DIR, "evaluator_data")
    EVALUATOR_NAME = f"{evaluator_base}_dev"

EVALUATOR_INPUT_PATH = os.path.join(EVALUATOR_DATA_DIR, input_file)
output_filename_base = f'{EVALUATOR_NAME}_predictions' # predictor_name will be added to this upon receiving the response

# Debug logs for validation
print(f"Using input file: {EVALUATOR_INPUT_PATH}")

# --- API Communication Settings ---
REQUEST_FORMAT = "application/json"
REQUEST_FORMAT = REQUEST_FORMAT.lower()

RESPONSE_FORMAT = "application/json"
RESPONSE_FORMAT = RESPONSE_FORMAT.lower()

# HTTP request retry
MAX_RETRIES = 50
RETRY_INTERVAL = 30 # Seconds to wait between each retry attempt