import math

# beta2b = 0.0
# r2rms = 0.2
# r1h = 0.045
# r1s = 0.14
# D2 = 0.4
# r2 = D2/2

# D1rms = math.sqrt((2*r1h)**2 + (2*r1s)**2)
# r2rms = math.sqrt((r1h)**2 + (r1s)**2)

# LOD = 0.5*((1.0-(D1rms/0.3048))/math.cos(beta2b / 180 * math.pi))
# print(f"Galvas\nLOD: {LOD:.4f}\tLh: {LOD*D2:.4f}")

# LOD = ((1.0-((2*r2rms)/0.3048))/math.cos(beta2b / 180 * math.pi))
# print(f"Radcomp\nLOD: {LOD:.4f}\tLh: {LOD*r2:.4f}")

mstar_m4 = 0.05
beta_1t = -63.0
t3 = 2.11 / 1000.0
s3 = 2*math.pi*(140.0/1000.0)/20.0
m4 = 0.35 * 0.2
error = 1.0

while(error > 10e-6):
    old = mstar_m4
    beta_inf3 = 1.0 - 5*mstar_m4**2
    dbeta_inf = -10.0* beta_inf3*mstar_m4
    mstar = mstar_m4 * m4
    m = 2*mstar - t3
    beta_infa = 1.0 - 5*(m/m4)**2
    beta_infb = 0.5 * (beta_1t + beta_infa)
    mstar_m4 = ((0.5 * s3) / ((1.0/math.tan(math.radians(beta_infa)))+math.tan(math.radians(beta_infb)))) + 0.5*t3/m4
    error = abs(mstar_m4 - old) / old
    print(f"mstar_m4: {mstar_m4}\terror: {error}")