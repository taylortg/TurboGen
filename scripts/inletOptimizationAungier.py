import matplotlib.pyplot as plt
import pandas as pd

url = './/..//tmp//optVelocities.csv'
df = pd.read_csv(url)
r1s = list = df['r1s'].tolist()
Cm = list = df['Cm'].tolist()
W1s = list = df['W1s'].tolist()
# Cm = pd.read_csv(url, index_col="Cm").tolist()
# W1s = pd.read_csv(url, index_col="W1s").tolist()


# plot shroud velocity triangle
fig, ax = plt.subplots()
plt.plot(Cm, W1s)
plt.xlabel("Cm (m/s)")
plt.ylabel("W1s (m/s)")
ymin = (round(min(W1s) / 10) * 10) - 10
ymax = (round(max(W1s) / 10) * 10) + 10
plt.ylim(ymin, ymax)
plt.xlim(40, 260)
ax.grid(which='major', linestyle='-', linewidth='0.5', color='black')
ax.minorticks_on()
ax.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
# plt.legend(loc='upper right')
plt.savefig("../out/Cm_W1s.png")
