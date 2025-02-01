import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv('.//..//tmp//velocities.csv')
shroud_data = data[data['Category']=='Shroud'][['X','Y']]
hub_data = data[data['Category']=='Hub'][['X','Y']]
outlet_data = data[data['Category']=='Outlet'][['X','Y']]

shroud_x = shroud_data['X'].tolist()
shroud_y = shroud_data['Y'].tolist()

hub_x = hub_data['X'].tolist()
hub_y = hub_data['Y'].tolist()

outlet_x = outlet_data['X'].tolist()
outlet_y = outlet_data['Y'].tolist()

# plot shroud velocity triangle
fig, ax = plt.subplots()
ax.quiver(shroud_x[0], shroud_y[0], shroud_x[1], shroud_y[1], angles='xy', scale_units='xy', scale=1, label='C', color='black')
ax.quiver(shroud_x[2], shroud_y[2], shroud_x[3], 0, angles='xy', scale_units='xy', scale=1, label='U', color='green')
ax.quiver(shroud_x[4], shroud_y[4], shroud_x[5], shroud_y[5], angles='xy', scale_units='xy', scale=1, label='W', color='blue')
plt.xlabel("(m/s)")
plt.ylabel("(m/s)")
plt.ylim(0,round(shroud_y[1])+10)
plt.xlim(-10,round(shroud_x[3])+10)
ax.grid(which='major', linestyle='-', linewidth='0.5', color='black')
ax.minorticks_on()
ax.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
plt.legend(loc='lower right')
plt.savefig("../out/vel_inlet_shroud.png")

# clear figure
plt.clf()

# plot hub velocity triangle
fig, ax = plt.subplots()
ax.quiver(hub_x[0], hub_y[0], hub_x[1], hub_y[1], angles='xy', scale_units='xy', scale=1, label='C', color='black')
ax.quiver(hub_x[2], hub_y[2], hub_x[3], 0, angles='xy', scale_units='xy', scale=1, label='U', color='green')
ax.quiver(hub_x[4], hub_y[4], hub_x[5], hub_y[5], angles='xy', scale_units='xy', scale=1, label='W', color='blue')
plt.xlabel("(m/s)")
plt.ylabel("(m/s)")
plt.ylim(0,round(hub_y[1])+10)
plt.xlim(-10,round(hub_x[3])+10)
ax.grid(which='major', linestyle='-', linewidth='0.5', color='black')
ax.minorticks_on()
ax.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
plt.legend(loc='lower right')
plt.savefig("../out/vel_inlet_hub.png")

# clear figure
plt.clf()

# plot outlet velocity triangle
fig, ax = plt.subplots()
ax.quiver(outlet_x[0], outlet_y[0], -outlet_x[1], outlet_y[1], angles='xy', scale_units='xy', scale=1, label='W', color='blue')
ax.quiver(-outlet_x[4], outlet_y[4], outlet_x[5], outlet_y[5], angles='xy', scale_units='xy', scale=1, label='Wθ', color='blue')
ax.quiver(outlet_x[2], outlet_y[2], outlet_x[3], outlet_y[3], angles='xy', scale_units='xy', scale=1, label='Cm/Wm', color='green')
ax.quiver(outlet_x[6], outlet_y[6], outlet_x[7], outlet_y[7], angles='xy', scale_units='xy', scale=1, label='C', color='black')
ax.quiver(outlet_x[8], outlet_y[8], outlet_x[9], outlet_y[9], angles='xy', scale_units='xy', scale=1, label='Cθ', color='black')
plt.xlabel("(m/s)")
plt.ylabel("(m/s)")
plt.ylim(0,round(outlet_y[1])+10)
plt.xlim(round(-outlet_x[1])-10, round(outlet_x[9])+10)
ax.grid(which='major', linestyle='-', linewidth='0.5', color='black')
ax.minorticks_on()
ax.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
plt.legend(loc='lower right')
plt.savefig("../out/vel_outlet.png")
