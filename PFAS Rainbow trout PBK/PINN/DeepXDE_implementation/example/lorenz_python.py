# -*- coding: utf-8 -*-
"""DeppXDE_Lorenz.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/150BmiF-xbNUGJQHRKytLEn3UFL-Dy8tS

Inverse problem for the Lorenz system with exogenous input

Example Source: https://deepxde.readthedocs.io/en/latest/demos/pinn_inverse/lorenz.inverse.forced.html
"""

import deepxde as dde
import numpy as np
import tensorflow as tf

# We also want to define our three unknown variables, 𝜎, 𝜌 and 𝛽
# which will now be called C1, C2, and C3, respectivly. These variables are given an initial guess of 1.0.

C1 = dde.Variable(1.0)
C2 = dde.Variable(1.0)
C3 = dde.Variable(1.0)

# Now we can begin by creating a TimeDomain class.

geom = dde.geometry.TimeDomain(0, 3)

def Lorenz_system(x, y):
    y1, y2, y3 = y[:, 0:1], y[:, 1:2], y[:, 2:]
    dy1_x = dde.grad.jacobian(y, x, i=0)
    dy2_x = dde.grad.jacobian(y, x, i=1)
    dy3_x = dde.grad.jacobian(y, x, i=2)
    return [
        dy1_x - C1 * (y2 - y1),
        dy2_x - y1 * (C2 - y3) + y2,
        dy3_x - y1 * y2 + C3 * y3,
    ]

def boundary(_, on_initial):
    return on_initial

ic1 = dde.icbc.IC(geom, lambda X: -8, boundary, component=0)
ic2 = dde.icbc.IC(geom, lambda X: 7, boundary, component=1)
ic3 = dde.icbc.IC(geom, lambda X: 27, boundary, component=2)

def gen_traindata():
    data = np.load("/Users/vassilis/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/PINN/DeepXDE_implementation/example/Lorenz.npz")
    return data["t"], data["y"]

observe_t, ob_y = gen_traindata()
observe_y0 = dde.icbc.PointSetBC(observe_t, ob_y[:, 0:1], component=0)
observe_y1 = dde.icbc.PointSetBC(observe_t, ob_y[:, 1:2], component=1)
observe_y2 = dde.icbc.PointSetBC(observe_t, ob_y[:, 2:3], component=2)

observe_t

# Now that the problem is fully setup, we define the PDE as:

data = dde.data.PDE(
    geom,
    Lorenz_system,
    [ic1, ic2, ic3, observe_y0, observe_y1, observe_y2],
    num_domain=400,
    num_boundary=2,
    anchors=observe_t,
)

# Where num_domain is the number of points inside the domain, and num_boundary is
# the number of points on the boundary. anchors are extra points beyond num_domain
# and num_boundary used for training. auxiliary_var_function is the interpolation
# function of 𝑓(𝑡)
# we defined above.

# Next, we choose the network. Here, we use a fully connected neural network of
# depth 4 (i.e., 3 hidden layers) and width 40:
net = dde.nn.FNN([1] + [40] * 3 + [3], "tanh", "Glorot uniform")
model = dde.Model(data, net)

external_trainable_variables =[C1, C2, C3]
# train adam
model.compile(
    "adam", lr=0.001, external_trainable_variables=external_trainable_variables
)

variable = dde.callbacks.VariableValue(
    external_trainable_variables, period=600, filename="variables.dat"
)
checkpointer = dde.callbacks.ModelCheckpoint(
    filepath= '/Users/vassilis/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/PINN/DeepXDE_implementation/example/model.ckpt',
    verbose=1, save_better_only=True, period=1, monitor='train loss')

losshistory, train_state = model.train(iterations=10000, callbacks=[variable, checkpointer])

dde.saveplot(losshistory, train_state, issave=True, isplot=True)

model.state_dict()

net = dde.nn.FNN([1] + [40] * 3 + [3], "tanh", "Glorot uniform")
model = dde.Model(data, net)

external_trainable_variables =[C1, C2, C3]

# train adam
model.compile(
    "adam", lr=0.001, external_trainable_variables=external_trainable_variables
)
model.restore('/Users/vassilis/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/PINN/DeepXDE_implementation/example/model.ckpt-10000.ckpt')

model.restore('/Users/vassilis/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/PINN/DeepXDE_implementation/example/model.ckpt-10000.ckpt', verbose=1)

model.state_dict()