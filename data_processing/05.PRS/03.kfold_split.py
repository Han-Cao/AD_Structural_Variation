import pandas as pd
from sklearn.model_selection import RepeatedKFold

fold=10
repeat=10

# k fold for AD
output_prefix = 'input/fold_AD/SRS.analysis_sample.fold'

df_sample = pd.read_csv("input/SRS.pheno_covar.txt", sep='\t')
df_sample = df_sample.dropna(subset="AD").reset_index(drop=True)
df_sample = df_sample.sample(frac=1).reset_index(drop=True)

# split within each diagnosis to ensure balanced case-control
df_disease = df_sample[df_sample['AD'] == 2].reset_index(drop=True)
df_control = df_sample[df_sample['AD'] == 1].reset_index(drop=True)

kf = RepeatedKFold(n_splits=fold, n_repeats=repeat, random_state=42)
for i, (train_index, test_index) in enumerate(kf.split(df_disease)):
    df_fold = df_disease.loc[test_index, ]
    df_fold[['FID','IID']].to_csv(output_prefix + str(i) + '.fam', sep='\t', index=False, header=True)

for i, (train_index, test_index) in enumerate(kf.split(df_control)):
    df_fold = df_control.loc[test_index, ]
    df_fold[['FID','IID']].to_csv(output_prefix + str(i) + '.fam', sep='\t', index=False, header=False, mode="a")
