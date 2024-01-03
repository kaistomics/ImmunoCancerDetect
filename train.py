# To generate an input matrix for model training, combine all previously calculated features.

import tensorflow as tf
from keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout
from sklearn.model_selection import train_test_split
import numpy as np
from numpy import loadtxt
