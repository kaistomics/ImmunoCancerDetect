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
