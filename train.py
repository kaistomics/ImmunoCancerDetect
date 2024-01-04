# To generate an input matrix for model training, combine all previously calculated features.

import tensorflow as tf
from keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout
from sklearn.model_selection import train_test_split
import numpy as np
from numpy import loadtxt

dat1 = loadtxt('case_train_input', delimiter=',')
dat2 = loadtxt('control_train_input', delimiter=',')
x_train = np.concatenate([dat1, dat2])
y_train = np.concatenate([np.ones(shape=(len(dat1))), np.zeros(shape=(len(dat2)))])

dat1 = loadtxt('case_test_input', delimiter=',')
dat2 = loadtxt('case_test_input', delimiter=',')
x_test = np.concatenate([dat1, dat2])
y_test = np.concatenate([np.ones(shape=(len(dat1))), np.zeros(shape=(len(dat2)))])

model = Sequential()
model.add(Dense(32, input_shape=(124,), activation='relu'))
model.add(Dense(8, activation='relu'))
# model.add(Dropout(0.1))
model.add(Dense(1, activation='sigmoid'))
adam = tf.keras.optimizers.Adam(lr=0.001, beta_1=0.8, beta_2=0.999, epsilon=None, decay=0.0, amsgrad=False)
model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])

model.fit(x_train, y_train, epochs=10, batch_size=128)

model.evaluate(x_test, y_test, verbose=2)

y_val_cat_prob=model.predict(x_test)

from sklearn.metrics import roc_curve,roc_auc_score
import matplotlib.pyplot as plt

fpr, tpr, thresholds = roc_curve (y_test, y_val_cat_prob)

def plot_roc_curve(fper, tper):
    plt.plot(fper, tper, color='red', label='ROC')
    plt.plot([0, 1], [0, 1], color='green', linestyle='--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend()
    plt.show()

plot_roc_curve (fpr,tpr) 
auc_score=roc_auc_score(y_test,y_val_cat_prob)
print(auc_score)
