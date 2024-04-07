#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

t = 1
U_PATH = "build/"
sample = 81
u_FTCS = np.loadtxt(U_PATH + "1d_FTCS_u.dat")
u_BTCS = np.loadtxt(U_PATH + "1d_BTCS_u.dat")
x = np.linspace(start=-1, stop=1, num=sample)

u_ana = np.empty_like(x)
u_diff_FTCS = np.empty_like(x)
u_diff_BTCS = np.empty_like(x)

for i in range(sample):
    u_ana[i] = -np.power(np.e, -t) * np.sin(np.pi * x[i])
    u_diff_FTCS[i] = np.abs(u_FTCS[i] - u_ana[i])
    u_diff_BTCS[i] = np.abs(u_BTCS[i] - u_ana[i])

plt.figure(figsize=(12, 5))
plt.subplot(121)
plt.plot(x, u_FTCS, "k--", label="FTCS simulation")  # simulation
plt.plot(x, u_BTCS, "b--", label="BTCS simulation")  # simulation
plt.plot(x, u_ana, "g-", label="analytical solution")
plt.xlabel("$x$")
plt.ylabel("$u$")
plt.legend()
plt.title("Solution Field")

plt.subplot(122)
plt.plot(x, u_diff_FTCS, label="FTCS simulation")
plt.plot(x, u_diff_BTCS, label="BTCS simulation")
plt.xlabel("$x$")
plt.ylabel("Difference $\epsilon$")
plt.title("Discretization error")
plt.legend()
plt.savefig("readu.pdf")
