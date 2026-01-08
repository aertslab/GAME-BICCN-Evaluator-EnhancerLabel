# BICCN Validated Enhancer Dataset

Evaluator container for the Genomic API for Model Evaluation (GAME). The system assesses computational
predictions of cell type-specific chromatin accessibility against experimentally validated enhancers from the mouse motor cortex.

https://zenodo.org/records/17802620

## Dataset Composition

The evaluation dataset comprises 171 experimentally validated cis-regulatory elements from the BRAIN Initiative Cell Census Network
(BICCN) across 19 mouse motor cortex cell types:

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

# Citations

This evaluator is based on validated enhancer data from:

Johansen, N.J., Kempynck, N., Zemke, N.R., Somasundaram, S., De Winter, S., Hooper, M., Dwivedi, D., Lohia, R., Wehbe, F., Li, B.,
Abaffyová, D., Armand, E.J., De Man, J., Ekşi, E.C., Hecker, N., Hulselmans, G., Konstantakos, V., Mauduit, D., Mich, J.K., Partel,
G., Daigle, T.L., Levi, B.P., Zhang, K., Tanaka, Y., Gillis, J., Ting, J.T., Ben-Simon, Y., Miller, J., Ecker, J.R., Ren, B., Aerts,
S., Lein, E.S., Tasic, B., and Bakken, T.E. (2025). Evaluating methods for the prediction of cell-type-specific enhancers in the
mammalian cortex. Cell Genomics 5(6), 100879. https://doi.org/10.1016/j.xgen.2025.100756

Ben-Simon, Y., Hooper, M., Narayan, S., Daigle, T. L., Dwivedi, D., Way, S. W., Oster, A., Stafford, D. A., Mich, J. K., Taormina, M. J., Martinez, R. A., Opitz-Araya, X., Roth, J. R., Alexander, J. R., Allen, S., Amster, A., Arbuckle, J., Ayala, A., Baker, P. M., . . . Ransford, S. (2025). A suite of enhancer AAVs and transgenic mouse lines for genetic access to cortical cell types. Cell, 188(11), 3045-3064.e23. https://doi.org/10.1016/j.cell.2025.05.002
