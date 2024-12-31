import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv('.\\..\\tmp\\velocities.csv')
shroud_data = data[data['Category']=='Shroud'][['X','Y']]
hub_data = data[data['Category']=='Hub'][['X','Y']]

shroud_x = shroud_data['X'].tolist()
shroud_y = shroud_data['Y'].tolist()

hub_x = hub_data['X'].tolist()
hub_y = hub_data['Y'].tolist()

plt.plot(shroud_x, shroud_y, label='Inlet shroud')
plt.plot(hub_x, hub_y, label='Inlet hub')

plt.legend()
plt.savefig("../out/vel1.png")
