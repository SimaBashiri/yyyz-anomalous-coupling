# Fitting Instructions

This guide provides instructions on how to perform fits using the Combine tool with the `MultiDimFit` method. You will perform fits at the expected level using the Asimov dataset and fit using an observed smeared toy dataset.

## Steps to Run the Fit

### 1. Move the `Mymodel.py` file

Move the `Mymodel.py` file to the `python` directory using the following command:

```bash
mv Mymodel.py python/  #copy the model to combine physics models
mkdir o  #workspace output
```

### 2. Fit at Expected Level Using the Asimov Dataset
```bash
text2workspace.py -P HiggsAnalysis.CombinedLimit.MyModel:yyyz datacard_Vertical_ZVertexCut.txt -o o/yyyz.root 
combine -M MultiDimFit o/yyyz.root --algo none --setParameters zta=0.1,zta~=0.01 -t -1
```
### 3. Fit Using an Observed Smeared Toy Dataset
```
combine -M MultiDimFit o/yyyz.root --algo none --setParameters zta=0.1,zta~=0.01,r=1 -t 1

```
