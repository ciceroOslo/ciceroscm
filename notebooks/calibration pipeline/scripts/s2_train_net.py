# %%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# %%
store = pd.HDFStore('../data/data.h5')
targ = store['targ']
parammat = store['parammat']
store.close()

# %%
from sklearn.model_selection import train_test_split
XTraining, XValidation, YTraining, YValidation = train_test_split(parammat,targ,test_size=0.1, shuffle=False) # before model building


# %%
import tensorflow as tf

from tensorflow import keras
from tensorflow.keras import layers

print(tf.__version__)
print(tf.config.list_physical_devices())

# %%
def build_and_compile_model(norm):
  model = keras.Sequential([
      norm,
      layers.Dense(512, activation='relu'),
      layers.Dense(32, activation='relu'),

      layers.Dense(nflds)
  ])

  model.compile(loss='mae',
                optimizer=tf.keras.optimizers.AdamW(0.001))
  return model

# %%
fflds=targ.columns.get_level_values(0)
fdates=targ.columns.get_level_values(1)
flds=targ.columns.get_level_values(0).unique()
dates=targ.columns.get_level_values(1).unique()
nflds=targ.shape[1]
nflds

# %%
X_train = tf.convert_to_tensor(XTraining, dtype=tf.float32)
Y_train = tf.convert_to_tensor(YTraining, dtype=tf.float32)
X_dev = tf.convert_to_tensor(XValidation, dtype=tf.float32)
Y_dev = tf.convert_to_tensor(YValidation, dtype=tf.float32)

# %%
normalizer = tf.keras.layers.Normalization(input_shape=[parammat.shape[1],], axis=-1)
normalizer.adapt(tf.convert_to_tensor((parammat)))
dnn_model = build_and_compile_model(normalizer)


# %%
from tqdm.keras import TqdmCallback

# %%
if 1:
    history=dnn_model.fit(X_train, Y_train, epochs=30000, batch_size=256, validation_data=(X_dev, Y_dev),verbose=0,callbacks=[TqdmCallback(verbose=0)])
    dnn_model.save_weights('../data/cicero_weights.h5')
else:
  dnn_model.load_weights('../data/cicero_weights.h5')

# %%
dnn_model.save('../data/dnn_cicero.keras')


# %%
pred = dnn_model.predict(X_dev)
pred_t = dnn_model.predict(X_train)

# %%
fig, ax = plt.subplots( len(flds), len(dates),figsize=(15, 10))

ax=ax.flatten()
for i in np.arange(nflds):
  ax[i].plot(pred_t[:,i],Y_train[:,i],'k.')
  ax[i].plot(pred[:,i],Y_dev[:,i],'r.')
  ax[i].set_title(fflds[i]+str(fdates[i]))
  ax[i].set_xlabel('predicted')
  ax[i].set_ylabel('observed')
plt.tight_layout()

# %%



