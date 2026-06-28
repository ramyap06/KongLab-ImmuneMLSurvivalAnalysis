# Codebase Review Questions

Walk me through the project. Start from the biological question you are trying to answer, then take me through the pipeline end to end. What are you predicting, what is your input data, and how does each stage feed into the next?

## DEG Filtering

In `deg.R`, you apply a triple threshold (p < 0.05, adj.p < 0.05, |logFC| >= 1). Why those three criteria together, and does the ordering matter?

## Univariate to Penalized Cox

In `cox_ai.R`, you go from univariate Cox to penalized Cox (LASSO/Elastic Net). What problem does each one solve that the previous could not?

## Risk Scoring Model Split

You fit a separate gene-only model for risk scoring rather than reading coefficients directly off the final multivariate model. Why?

## Gene Set Construction

In `load_data.r`, how was the immune gene catalog built and what does "immune genes" mean in this context?

## ml.R

`ml.R` is currently a stub. What is the plan for it?
