# kong-lab-ml-project
Machine learning research project developed as a research assistant to study immune gene signatures in breast cancer.

## Sample ReadME

### Package Installation

The following libaries are necessary to run these files:
- tidyverse
- GEOquery
- affy
- limma
- survival

etc

### Experimentation Order

The following files were created & run in this order to produce an output:

1. load_data.ipynb
2. preprocess.ipynb
3. deg_analysis.ipynb
4. cox_survival_analysis.ipynb

### Notes

- Files labeled {name}_learning.ipynb are not used to produce output but simply to independently learn & pratice a particular process/concept.

- To load data, simply change the curr variable to the name of the dataset you want to load and rerun the notebook.