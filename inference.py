# # With min max scaling
#
#  Coefficients:
#  [[-0.84150197 -0.96470767 -1.05068263]]
#  Intercept:
#  [1.18389693]
#
# # Dividing by max possible val 
#
#  Coefficients:
#  [[-0.74944314 -0.88844054 -1.11174202]]
#  Intercept:
#  [1.12997506]
#
# # # # # # #

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import itertools
from sklearn.model_selection import KFold
from sklearn.svm import SVR
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from scipy.stats import pearsonr

def combinations(array_a, array_b):
  return list(map(lambda tpl: '_'.join(tpl), itertools.product(array_a, array_b)))

input_file = 'full_chunk-8_2023-01-20 21:47:35.174.csv'
#input_file = 'full_chunk-16_2023-01-20 22:43:23.476.csv'
#input_file = 'full_chunk-32_2023-01-20 22:52:03.383.csv'
#input_file = 'full_chunk-64_2023-01-20 22:56:13.456.csv'

data = pd.read_csv("eeg/erp/%s" % input_file)

events_threshold = 0
selected_components = ['P100', 'N100', 'P200', 'N200', 'P300']
selected_attributes = ['mean']

def create_model():
  if model_type == 'linear':
    return LinearRegression()
  elif model_type == 'svm':
    return SVR(kernel='rbf')
  else:
    raise Exception("Unknown model type %s. Use 'linear' or 'svm'" % model_type)

def analyze_neupsilin(aphasia):
  neupsilin_cols = ['digit_ordering', 'listening_span', 'verbal_evocation']

  subset = data[data['has_aphasia'] == aphasia].drop_duplicates(subset=['person'])

  print("Means for aphasia = ", aphasia)
  print(subset[neupsilin_cols].mean())

def check_acceptance_criterion(target_column, model_rmse):
  subset = data[data['has_aphasia'] == False].drop_duplicates(subset=['person'])
  threshold = np.std(subset[target_column].values)

  delta = threshold - model_rmse

  if delta > 0:
    print("Great, the model has less error than there is variance in W!")
  else:
    print("Great, the criterion failed. The error of the model is greater than the variance in W")
  
  print("Threshold: %.5f" % threshold)
  print("Model RMSE: %.5f" % model_rmse)
  print("Difference: %.5f" % delta)

def create_target_column(col_name):
  coefficients = {
    'digit_ordering': 0.74944314,
    'listening_span': 0.88844054,
    'verbal_evocation': 1.11174202
  }

  data[col_name] = 0

  for key, value in coefficients.items():
    data[col_name] += data[key] * value

def analyze_by_channels():
  channels = data['channel'].unique()

  model_name = ({
    'svm': 'SVM',
    'linear': 'Linear'
  })[model_type]

  out_data = []

  for channel in channels:
    rows_of_channel = data[data['channel'] == channel]
    selected_rows = rows_of_channel[rows_of_channel['n_events'] > events_threshold]

    X_cols = combinations(selected_components, selected_attributes)

    X = selected_rows[X_cols].values
    y = selected_rows['ground_truth'].values

    kf = KFold(n_splits=5)

    all_rmse = []

    ch_y_true = []
    ch_y_pred = []

    for i, (train_index, test_index) in enumerate(kf.split(X)):
      print("=== Fold %d (training on %d samples, testing on %d) ===" % (i, len(train_index), len(test_index)))
      X_train, y_train = X[train_index], y[train_index]
      X_test, y_test = X[test_index], y[test_index]

      model = create_model()
      model.fit(X_train, y_train)

      y_pred = model.predict(X_test)

      rmse = mean_squared_error(y_test, y_pred, squared=False)

      print("RMSE for channel %d fold %d: %.6f" % (channel, i, rmse))

      all_rmse.append(rmse)

      ch_y_true.extend(y_test)
      ch_y_pred.extend(y_pred)

    average_rmse = np.mean(all_rmse)
    print("=> Average RMSE for channel %d: %.6f" % (channel, average_rmse))

    check_acceptance_criterion('ground_truth', average_rmse)

    plt.suptitle("[%s] W: real vs estimativa - canal %s" % (model_name, channel))
    plt.title("RMSE = %.5f" % average_rmse)
    plt.plot(ch_y_true, 'bo', label="Real")
    plt.plot(ch_y_pred, 'rx', label="Estimativa")
    plt.grid()
    plt.legend()
    plt.savefig('eeg/erp/img/inference/%s_real_vs_estimate_channel%s.png' % (model_name, channel))
    plt.close()

    out_data.append({
      'channel': channel,
      'rmse': average_rmse,
    })

  out_df = pd.DataFrame.from_records(out_data)
  out_df.to_csv('eeg/erp/inference/result_%s_%s.csv' % (model_type, input_file))


def analyze_combined():
  model_name = ({
    'svm': 'SVM',
    'linear': 'Linear'
  })[model_type]

  out_data = []
  selected_rows = data[data['n_events'] > events_threshold]

  X_cols = combinations(selected_components, selected_attributes)

  X = selected_rows[X_cols].values
  y = selected_rows['ground_truth'].values

  kf = KFold(n_splits=5)

  all_rmse = []

  ch_y_true = []
  ch_y_pred = []

  for i, (train_index, test_index) in enumerate(kf.split(X)):
    print("=== Fold %d (training on %d samples, testing on %d) ===" % (i, len(train_index), len(test_index)))
    X_train, y_train = X[train_index], y[train_index]
    X_test, y_test = X[test_index], y[test_index]

    model = create_model()
    model.fit(X_train, y_train)

    y_pred = model.predict(X_test)

    rmse = mean_squared_error(y_test, y_pred, squared=False)

    print("RMSE for fold %d: %.6f" % (i, rmse))

    all_rmse.append(rmse)

    ch_y_true.extend(y_test)
    ch_y_pred.extend(y_pred)

  average_rmse = np.mean(all_rmse)
  print("=> Average RMSE: %.6f" % average_rmse)

  check_acceptance_criterion('ground_truth', average_rmse)

  plt.suptitle("[%s] W: real vs estimativa" % model_name)
  plt.title("RMSE = %.5f" % average_rmse)
  plt.plot(ch_y_true, 'bo', label="Real")
  plt.plot(ch_y_pred, 'rx', label="Estimativa")
  plt.grid()
  plt.legend()
  plt.show()


def analyze_features():
  features = combinations(selected_components, selected_attributes)

  for feature in features:
    ground_truth_vals = data['ground_truth'].values
    feat_vals = data[feature].values

    r, p = pearsonr(ground_truth_vals, feat_vals)

    feat_vals = feat_vals / np.max(np.abs(feat_vals))
    feat_vals = feat_vals * np.max(np.abs(ground_truth_vals))
    
    plt.title("%s vs ground truth (r = %.4f)" % (feature, r))
    plt.plot(ground_truth_vals, 'b.', label="Ground truth")
    plt.plot(feat_vals, 'r.', label=feature)
    plt.legend()
    plt.savefig("eeg/erp/img/features/%s vs ground truth.png" % feature)
    plt.close()

create_target_column('ground_truth')
data = data.sort_values(by=['ground_truth'])

model_type = 'svm'

#analyze_by_channels()
analyze_combined()

model_type = 'linear'

#analyze_by_channels()
analyze_combined()