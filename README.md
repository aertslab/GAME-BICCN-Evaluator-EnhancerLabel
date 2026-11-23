# BICCN Validated Enhancer Dataset

Evaluator container for the Genomic API for Model Evaluation (GAME). The system assesses computational
predictions of cell type-specific chromatin accessibility against experimentally validated enhancers from the mouse motor cortex.

https://zenodo.org/records/17802620

## Dataset Composition

The evaluation dataset (`evaluator_data/biccn_enhancers_withSequence.csv`) comprises 171 experimentally validated enhancers from the
BRAIN Initiative Cell Census Network (BICCN). Enhancers are classified across 19 mouse motor cortex cell types:

- Astrocytes (Astro)
- Endothelial cells (Endo)
- Excitatory neurons: L2/3 IT, L5 ET, L5 IT, L5/6 NP, L6 CT, L6 IT, L6b
- Inhibitory neurons: Lamp5, Pvalb, Sncg, Sst, Sst Chodl, Vip
- Glial cells: Microglia/PVM, OPC, Oligodendrocytes
- Vascular cells: VLMC

Each enhancer sequence is annotated with its target cell type, specificity classification (on-target only),
and activity strength (strong or weak). Sequences range from approximately 200-500 base pairs and were functionally validated using
massively parallel reporter assays (MPRA).

## Evaluation Metrics

Cell Type Classification: The evaluator computes multiclass classification metrics to assess the model's ability to correctly predict
which cell type each enhancer regulates:

- Accuracy: Overall fraction of correctly classified enhancers
- Precision (weighted): Weighted average precision across cell types, accounting for class imbalance
- Recall (weighted): Weighted average recall across cell types
- F1-score (weighted): Weighted harmonic mean of precision and recall

## Prediction Request

The Evaluator sends one request with `"readout": "point"` and the 171 sequences keyed by `enhancerID`, with no flanks or
`prediction_ranges`. It asks for one task per cell type (the 19 names in `CELL_TYPES` in `config.py`), each with
`type: "accessibility"`, `species: "mus_musculus"`, `scale: "linear"` and the name `accessibility_<cell type>`.
An enhancer's predicted cell type is the task with the highest prediction for it, which is compared with its `target_ct`.

Because every cell type is needed for that comparison, the run is not scored (all metrics are `NaN`) if the Predictor does not
return HTTP 200, if any task is missing, duplicated or failed, if a task's enhancer IDs differ from those sent, or if any
prediction is not a finite number.

## Build and run

Building is optional: you only need to build the container if you don't already have `biccn_evaluator.sif`, or if you have
changed the code or `evaluator_data/`. Build from the repository root (the build pulls `python:3.13-slim` and pip packages,
so it needs internet access):

```bash
apptainer build biccn_evaluator.sif evaluator.def
```

Run the Evaluator against a running Predictor:

```bash
apptainer run biccn_evaluator.sif <PREDICTOR_IP> <PREDICTOR_PORT> <OUTPUT_DIR>
```

The Evaluator calls the Predictor's `/formats` and then `/predict`, and writes to `<OUTPUT_DIR>` (created if missing):

- `<EVALUATOR_NAME>_predictions_from_<predictor_name>.json`: the Predictor's response
- `evaluation_summary_<EVALUATOR_NAME>.csv`: tab-separated, one row per metric, appended on each run (values are `NaN`
  if the run is not scored)

The Evaluator names itself `DeepBICCN2_EnhancerLabel_Evaluator_<build timestamp>` inside the container (from the Apptainer
build date) and `DeepBICCN2_EnhancerLabel_Evaluator_dev` outside it. If the Predictor cannot be reached, the run exits with an
error and writes nothing.

# Citations

This evaluator is based on validated enhancer data from:

Johansen, N.J., Kempynck, N., Zemke, N.R., Somasundaram, S., De Winter, S., Hooper, M., Dwivedi, D., Lohia, R., Wehbe, F., Li, B.,
Abaffyová, D., Armand, E.J., De Man, J., Ekşi, E.C., Hecker, N., Hulselmans, G., Konstantakos, V., Mauduit, D., Mich, J.K., Partel,
G., Daigle, T.L., Levi, B.P., Zhang, K., Tanaka, Y., Gillis, J., Ting, J.T., Ben-Simon, Y., Miller, J., Ecker, J.R., Ren, B., Aerts,
S., Lein, E.S., Tasic, B., and Bakken, T.E. (2025). Evaluating methods for the prediction of cell-type-specific enhancers in the
mammalian cortex. Cell Genomics 5(6), 100879. https://doi.org/10.1016/j.xgen.2025.100756

Ben-Simon, Y., Hooper, M., Narayan, S., Daigle, T. L., Dwivedi, D., Way, S. W., Oster, A., Stafford, D. A., Mich, J. K., Taormina, M. J., Martinez, R. A., Opitz-Araya, X., Roth, J. R., Alexander, J. R., Allen, S., Amster, A., Arbuckle, J., Ayala, A., Baker, P. M., . . . Ransford, S. (2025). A suite of enhancer AAVs and transgenic mouse lines for genetic access to cortical cell types. Cell, 188(11), 3045-3064.e23. https://doi.org/10.1016/j.cell.2025.05.002
