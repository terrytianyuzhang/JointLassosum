# Synthetic-data based hyperparameter tuning

Here is the source code for synthetic-data based cross-validation (CV). This method is general and we take JLS as an example. We are doing some approximate CV when the original individual-level data is not accessible.

## Step 0: fit a JLS model

One of the inputs for this pipeline is a preliminary estimated lasso regression coefficient with a pre-specified hyperparameter combination. I put it in the [input folder](/input/).






