import numpy as np

'''
PURPOSE:
--------
Generate various source-centred circular annulus region files to be used in the pile-up analysis.

RcutA <-- This tag means an annulus with Rcut = A pixels and Rout = 30 pixels
RA-B  <-- This tag means an annulus of width dR = 1 pixel with (Rin, Rout) = (A, B) pixels

The "RcutA" variety are used for pile-up analysis of grade distributions and spectral distortions.
The "RA-B" variety are used for pile-up analysis of radial intensity profiles.

Ultimately, we will use the "Rcut" variety for the scientific analysis.
'''

#====================================================================================================
# GRO J1655-40 (RA, DEC)
ra  = '16:54:00.137'
dec = '-39:50:44.900'

#====================================================================================================
# RcutA REGIONS

# Number of source regions to create with Rcut = 0,1,2,...,Nreg-1 pixels
Nreg = 26

# Specify the inner/outer annulus radii for the source region extraction
pix2arcsec  = 2.357                     # Conversion factor [arcsec / pixel]
Rcut_pixels = np.arange(Nreg)           # Inner radius of the extracted annulus [pixels]
Rcut_arcsec = Rcut_pixels * pix2arcsec  # Inner radius of the extracted annulus [arcsec]
Rout_pixels = 30                        # Outer radius of the extracted annulus [pixels]
Rout_arcsec = Rout_pixels * pix2arcsec  # Outer radius of the extracted annulus [arcsec]

# Region filenames
outdir   = '/Users/salvesen/research/fcolabs/data/regions/'
regfiles = []
for i in range(Nreg):
    regfiles.append(outdir + 'source_annulus_Rcut' + str(Rcut_pixels[i]) + '.reg')

# Loop through each circular annulus region to be extracted and write out the region file
for i in range(Nreg):

    # Region file header
    lines = []
    lines.append('# Region file format: DS9 version 4.1')
    lines.append('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1')
    lines.append('fk5')

    # Region
    Rcut_str = '{:.3f}'.format(Rcut_arcsec[i])
    Rout_str = '{:.3f}'.format(Rout_arcsec)
    region   = 'annulus(' + ra + ',' + dec + ',' + Rcut_str + '"' + ',' + Rout_str + '")'
    lines.append(region)

    # Write out the region file
    f = open(regfiles[i], 'w')
    for item in lines:
        f.write(item + '\n')
    f.close()

#====================================================================================================
# RA-B REGIONS

# Number of source regions to create with Rin = 0,1,2,...,Nreg-1 and Rout = 1,2,3,...,Nreg pixels
Nreg = 26

# Width of an annulus [pixels]
dR = 1.0

# Specify the inner/outer annulus radii for the source region extraction
pix2arcsec  = 2.357                     # Conversion factor [arcsec / pixel]
Rin_pixels  = np.arange(Nreg)           # Inner radius of the extracted annulus [pixels]
Rin_arcsec  = Rin_pixels * pix2arcsec   # Inner radius of the extracted annulus [arcsec]
Rout_pixels = np.arange(Nreg) + dR      # Outer radius of the extracted annulus [pixels]
Rout_arcsec = Rout_pixels * pix2arcsec  # Outer radius of the extracted annulus [arcsec]

# Region filenames
outdir   = '/Users/salvesen/research/fcolabs/data/regions/'
regfiles = []
for i in range(Nreg): regfiles.append(outdir + 'source_annulus_R' + str(Rin_pixels[i]) + '-' + str(Rout_pixels[i]) + '.reg')

# Loop through each circular annulus region to be extracted and write out the region file
for i in range(Nreg):

    # Region file header
    lines = []
    lines.append('# Region file format: DS9 version 4.1')
    lines.append('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1')
    lines.append('fk5')

    # Region
    Rin_str  = '{:.3f}'.format(Rin_arcsec[i])
    Rout_str = '{:.3f}'.format(Rout_arcsec[i])
    region   = 'annulus(' + ra + ',' + dec + ',' + Rin_str + '"' + ',' + Rout_str + '")'
    lines.append(region)

    # Write out the region file
    f = open(regfiles[i], 'w')
    for item in lines:
        f.write(item + '\n')
    f.close()

#====================================================================================================
