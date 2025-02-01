import math

beta2b = 0.0
r2rms = 0.2
r1h = 0.045
r1s = 0.14
D2 = 0.4
r2 = D2/2

D1rms = math.sqrt((2*r1h)**2 + (2*r1s)**2)
r2rms = math.sqrt((r1h)**2 + (r1s)**2)

LOD = 0.5*((1.0-(D1rms/0.3048))/math.cos(beta2b / 180 * math.pi))
print(f"Galvas\nLOD: {LOD:.4f}\tLh: {LOD*D2:.4f}")

LOD = ((1.0-((2*r2rms)/0.3048))/math.cos(beta2b / 180 * math.pi))
print(f"Radcomp\nLOD: {LOD:.4f}\tLh: {LOD*r2:.4f}")
