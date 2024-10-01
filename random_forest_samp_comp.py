#!/usr/bin/python 3

#details about his script here - https://docs.google.com/document/d/1uldwa4ep4Pr1NrD6hmcmLY3AL-vCiniPYVSswhX1kF8/edit

import sklearn as sk
from sklearn import preprocessing
import pandas as pd
import matplotlib.pyplot as plt


#These paths need changing
raw_data = pd.read_csv('../taxonomy_totalread.csv')
seq_info = pd.read_csv('../sequencing_details_samp_comp.csv')

# Create a column of the read counts divided by the read counts for that barcode.
#Don't think we need this for analysis
raw_data['read_norm_by_total'] = raw_data.read_count/raw_data.total_read

# Want to merge these two datasets based on barcode
data = raw_data.merge(seq_info, on = 'barcode', how ='left')

# We want to aggregate information by barcode 
# Number of unique taxaID, also keep the other info

data = data.groupby(['barcode'], as_index=False).aggregate({
    'taxaID': pd.Series.nunique,
    'location': 'first',
    'total_read': 'first',
    'Sampler' : 'first',
    'Time' : 'first',
    'bp' : 'first',
    'MeanLength' : 'first',
    'Shortest' : 'first',
    'Longest' : 'first',
    'N50Length' : 'first'
}).rename(columns={'taxaID': 'num_unique_taxa'})

data = data.replace(',','', regex=True)

data.barcode = data.barcode.astype('category')
data.location = data.location.astype('category')

#normalising the data (making it from 0-1), want each col to be normalised on it's own
data[['num_unique_taxa', 'total_read', 'MeanLength', 'Shortest', 'Longest', 'N50Length']] = preprocessing.normalize(data[['num_unique_taxa', 'total_read', 'MeanLength', 'Shortest', 'Longest', 'N50Length']], norm='max', axis=0)

from sklearn.decomposition import PCA

#get dummies is one hot encoding the data, so it's 1's and 0's, then joining this to our normalised data
data = data[['num_unique_taxa', 'total_read', 'MeanLength', 'Shortest', 'Longest', 'N50Length']].join(pd.get_dummies(data.location)).join(pd.get_dummies(data.Sampler).join(pd.get_dummies(data.Time)))
data = data.rename(columns={'50': '50sec', '25' : '25sec'})
data.columns = data.columns.astype(str)
y = data.num_unique_taxa #target 
X = data.loc[:, data.columns != 'num_unique_taxa'] #features (everything else)
#Would split the data if we had lots more, to test how good it is 
#train_test_split(X, y, test_size=0.2, random_state=0)

#pca = PCA()
#pca_out = pca.fit_transform(X)
#explained_variance = pca.explained_variance_ratio_


# evaluate random forest ensemble for regression
import pydot
from numpy import mean
from numpy import std
from sklearn.datasets import make_regression
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedKFold
from sklearn.ensemble import RandomForestRegressor
from sklearn.tree import export_graphviz


# define the model
model = RandomForestRegressor()
model.fit(X,y)

# evaluate the model
cv = RepeatedKFold(n_splits=10, n_repeats=3, random_state=1)
n_scores = cross_val_score(model, X, y, scoring='neg_mean_absolute_error', cv=cv, n_jobs=-1, error_score='raise')
# report performance
print('MAE: %.3f (%.3f)' % (mean(n_scores), std(n_scores)))

estimator = model.estimators_[0]

# Export as dot file
export_graphviz(estimator, out_file='tree.dot', 
                feature_names = data.columns[1:],
                rounded = True, proportion = False, 
                precision = 2, filled = True)

from subprocess import check_call
check_call(['dot','-Tpng','tree.dot','-o','not_barcode_all_variables.png'])