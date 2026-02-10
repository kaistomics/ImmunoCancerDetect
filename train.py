# to generate input matrices for model training, combine all previously calculated features.

import tensorflow as tf
from tensorflow.keras.layers import Dense, Dropout
from keras.models import Sequential
from keras.callbacks import EarlyStopping
from sklearn.model_selection import train_test_split, KFold
import numpy as np
from numpy import loadtxt

# import input matrices
dat1 = loadtxt('case_train_input', delimiter=',')
dat2 = loadtxt('control_train_input', delimiter=',')
x_train = np.concatenate([dat1, dat2])
y_train = np.concatenate([np.ones(shape=(len(dat1))), np.zeros(shape=(len(dat2)))])

dat1 = loadtxt('case_test_input', delimiter=',')
dat2 = loadtxt('control_test_input', delimiter=',')
x_test = np.concatenate([dat1, dat2])
y_test = np.concatenate([np.ones(shape=(len(dat1))), np.zeros(shape=(len(dat2)))])

# hyperparameters
num_folds = 5
num_epoch = 10
num_batchsize = 128
early_stopping = EarlyStopping(monitor='val_loss', patience=5, restore_best_weights=True)

# define the K-fold Cross Validator
kfold = KFold(n_splits=num_folds, shuffle=True)

# define per-fold score containers
acc_per_fold = []
loss_per_fold = []
y_prob = []

# K-fold Cross Validation model evaluation
fold_no = 1
for train, test in kfold.split(x_train, y_train):
    
    # construct the model
    model = Sequential()
    model.add(Dense(32, input_shape=(124,), activation='relu'))
    model.add(Dense(8, activation='relu'))
    # model.add(Dropout(0.1))
    model.add(Dense(1, activation='sigmoid'))
    adam = tf.keras.optimizers.Adam(lr=0.001, beta_1=0.8, beta_2=0.999, epsilon=None, decay=0.0, amsgrad=False)
   
    model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])

    # generate a print
    print('------------------------------------------------------------------------')
    print(f'Training for fold {fold_no} ...')
    
    # fit data to the model
    model.fit(x_train[train], y_train[train], epochs=num_epoch, batch_size=num_batchsize, callbacks = [early_stopping], validation_split=0.2)
    
    # evaluate the model
    scores = model.evaluate(x_train[test], y_train[test], verbose=0)
    print(f'Score for fold {fold_no}: {model.metrics_names[0]} of {scores[0]}; {model.metrics_names[1]} of {scores[1]*100}%')
    acc_per_fold.append(scores[1] * 100)
    loss_per_fold.append(scores[0])
    y_prob.append(model.predict(x_test))
    
    # increase fold number
    fold_no = fold_no + 1

# provide average scores
print('------------------------------------------------------------------------')
print('Score per fold')
for i in range(0, len(acc_per_fold)):
  print('------------------------------------------------------------------------')
  print(f'> Fold {i+1} - Loss: {loss_per_fold[i]} - Accuracy: {acc_per_fold[i]}%')
print('------------------------------------------------------------------------')
print('Average scores for all folds:')
print(f'> Accuracy: {np.mean(acc_per_fold)} (+- {np.std(acc_per_fold)})')
print(f'> Loss: {np.mean(loss_per_fold)}')
print('------------------------------------------------------------------------')

from sklearn.metrics import roc_curve,roc_auc_score
import matplotlib.pyplot as plt

best_fold = np.argmax(acc_per_fold)

# predict with model of best performance
fpr, tpr, thresholds = roc_curve (y_test, y_prob[best_fold])

def plot_roc_curve(fper, tper):
    plt.plot(fper, tper, color='red', label='ROC')
    plt.plot([0, 1], [0, 1], color='green', linestyle='--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend()
    plt.show()

# plot the roc curve
plot_roc_curve(fpr,tpr) 
print('------------------------------------------------------------------------')
print(f'Final AUC: {roc_auc_score(y_test,y_prob[1])}')
