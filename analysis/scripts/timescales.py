import numpy as np
import constants as c

# See Guan & Gammie (2009)

# Calculate the Keplerian speed [cm s^-1]
# R - Radius [cm]
# M - Mass [Msun]
def vKeplerian(R, M):
    vK = np.sqrt(c.G * M * c.sun2g / R)
    return vK

# Calculate the diffusion timescale
def tdiffuse(R, M, alpha=0.1, H_R=0.01):
    vK = vKeplerian(R=R, M=M)
    t  = R / (alpha * H_R * vK)
    return t

# Binary orbital separation for GRO J1655-40

M     = 5.4  # [Msun]
Rg    = c.G * M * c.sun2g / c.c**2
R     = 100.0 * Rg  # [cm]
alpha = 0.1
H_R   = 0.01

tsec = tdiffuse(R=R, M=M, alpha=alpha, H_R=H_R)
print tsec
