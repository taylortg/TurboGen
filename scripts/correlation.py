import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

data = pd.read_csv('.//..//tmp//correlation.csv')

Mu2 = data[data['Category']=='Mu2']['Data'].item()
gamma = data[data['Category']=='gamma']['Data'].item()
eta = data[data['Category']=='eta_p']['Data'].item()

PR = np.linspace(1, 3, 201)
lamb = np.linspace(0.5, 0.8, 101)

workCoeff = (PR**((gamma-1)/(eta*gamma))-1)/((gamma-1)*Mu2**2)

plt.plot(PR, workCoeff, label="$\lambda$ = f(PR)")
plt.grid(True)
plt.xlabel("Pressure ratio, PR")
plt.ylabel("Work Coefficient, $\lambda$")
plt.legend()
plt.savefig("../out/work_coeff.png")

# clear figure
plt.clf()