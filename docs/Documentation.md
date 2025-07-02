# Documentation

## Input file


### Header
```toml
title = "Text case"
version = "1.0.0"

```

$$y=a\,b$$
*a*

[link](http:\\www.google.com)

[unit]
scale = 1e-3

[unit]
Definition of the scale for all geometrical inputs. If [`unit`] label is not given, scale = 1 is assumed, that is all geometrical data are assumed to be in metres.

```toml
<Dir>\<File>

```
---
- list1
- list2
  - list2.1




[pm]
files = ['\\filesrv\HPC_GAUSS\A4IM\Magnets\RawFiles_v21\ring0.txt',
         '\\filesrv\HPC_GAUSS\A4IM\Magnets\RawFiles_v21\ring1.txt',
         '\\filesrv\HPC_GAUSS\A4IM\Magnets\RawFiles_v21\ring2.txt',
         '\\filesrv\HPC_GAUSS\A4IM\Magnets\RawFiles_v21\ring3.txt',
         '\\filesrv\HPC_GAUSS\A4IM\Magnets\RawFiles_v21\ring4.txt',
         '\\filesrv\HPC_GAUSS\A4IM\Magnets\RawFiles_v21\ring5.txt',
         '\\filesrv\HPC_GAUSS\A4IM\Magnets\RawFiles_v21\ring6.txt',
         '\\filesrv\HPC_GAUSS\A4IM\Magnets\RawFiles_v21\ring7.txt',
         '\\filesrv\HPC_GAUSS\A4IM\Magnets\RawFiles_v21\ring8.txt',
         '\\filesrv\HPC_GAUSS\A4IM\Magnets\RawFiles_v21\mirror\ring0.txt',
         '\\filesrv\HPC_GAUSS\A4IM\Magnets\RawFiles_v21\mirror\ring1.txt',
         '\\filesrv\HPC_GAUSS\A4IM\Magnets\RawFiles_v21\mirror\ring2.txt',
         '\\filesrv\HPC_GAUSS\A4IM\Magnets\RawFiles_v21\mirror\ring3.txt',
         '\\filesrv\HPC_GAUSS\A4IM\Magnets\RawFiles_v21\mirror\ring4.txt',
         '\\filesrv\HPC_GAUSS\A4IM\Magnets\RawFiles_v21\mirror\ring5.txt',
         '\\filesrv\HPC_GAUSS\A4IM\Magnets\RawFiles_v21\mirror\ring6.txt',
         '\\filesrv\HPC_GAUSS\A4IM\Magnets\RawFiles_v21\mirror\ring7.txt',
         '\\filesrv\HPC_GAUSS\A4IM\Magnets\RawFiles_v21\mirror\ring8.txt']
mat_code = [1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,2]

[pm]
Label which defines the input data of permanent magnets (PMs).
files: is the array of filenames which contain the PM data. Each file represents a group of PM, whose code is defined with the variable mat_code.
E.g. all magnets described in file "\\filesrv\HPC_GAUSS\A4IM\Magnets\RawFiles_v21\ring8.txt" have a group index equal to 2.
The description of the file data is given below.



[pm.modgeo]
active = true
code = 3
zpos = [9.875,29.875,49.875,69.625,89.875,109.626,129.875,149.625,187,-10.125,-30.125,-50.125,-69.875,-89.625,-109.625,-129.625,-149.625,-187]
file = '\\filesrv\HPC_GAUSS\A4IM\Torque\Momento_torcente_Dipole.mat'
angle = 1.1

[pm.modgeo]
Label which activate modification of the geometrical data of PMs.
If variable active = true means that the modification is active.
code: define the type of modification
code = 1: only modification of the position of PM groups along z-axis
code = 2: only modification of the orientation of PMs based on torque
code = 3: both previous modifications are active

if code = 1 or code = 3
zpos: new positions along z-axis of PMs group

if code = 2 or code = 3
file: name lf the matlab file containing the torque values for the entire list or PMs
angle: variation of PM orientation in degrees with respect to the original orientation (+/- depending on the sign of the torque)









[material]
JHcurve = true
MuR = [0.0223,0.0168]
Jresidua = [1.4311,1.4616]
files = ['T:\A4IM\Script\N52_estratto_HJ_23.txt',
         'T:\A4IM\Script\N52_long_estratto_HJ.txt']
tref = [18,18]
jvt = [-0.001263,-0.001263]
hvt1 = [-1e-2,-1e-2]
hvt2 = [3.8e-5,3.8e-5]



[material]
Label with the description of the PM material magnetic properties
JHcurve: flag which determines if linear or nonlinear magnetic properties are given.

if JHcurve = true (nonlinear behaviour)
   files: list of file names with the HJ curves of the PM materials
   tref: reference temperatures (in celsius degrees) corresponding to the HJ curves described in the previous files
   jvt: coefficients for linear rescaling values of J 
   hvt1,hvt2: linear and quadratic coefficients for rescaling values of H

if JHcurve = false (linear behaviour)
   MuR: relative permeability of the PM materials  (B = Mu0 H + J = Mu0 H + Jr + MuR Mu0 H)
   Jresidua: residual polarization (in tesla) of PM materials




[points]
file = 'T:\A4IM\Script\DSV_MeasurementPoints.txt'

[points]
Label describing the coordinates of computational points.
file: name of the text file with the list of x,y,z coordinates


[temperature]
type = 1
value = 23.7
file = 'T:\A4IM\Thermal_effects\Magnet_core_GC100W_BCTemp20\Export\3d\temperature.csv'

[temperature]
Label with the information about the PM temperature values
type: code describing the temperature status of PMs
type = 0: all PMs are at the reference temperature (non modification of the material data)
type = 1: all PMs are at the same temperature (given by value). The HJ curve are modified accordingly to the jvt and hvt1/hvt2 coefficients
type = 2: each PM has a different temperature (as listed in the .csv file). The HJ curve of each PM is modified accordingly to the jvt and hvt1/hvt2 coefficients


[simul]
type = 'det'
component = 'x'
outputfile = 'pippo.mat'
reaction = true
NL = true

[simul]
Label defining the simulation data
type: type of simulation
type = 'det' deterministic, that is one simulation using the reference input data
type = 'MC' stochastic, that is Montecarlo simulation using the variability of input data provided in the label [montecarlo]
component: cartesian component of the B field to be considered ('x', 'y' or 'z'). Also the magnitude ('mag') can be selected.
outputfile: name of the matlab output file where results of simulations will be stored
reaction: flag to activate or deactivate the PM reactions
NL: flag to activate or deactivate the nonlinear simulation


[torque]
compute=true
outputfile='papero'

[torque]
Label to activate the computation of torque acting on PMs.
Data written on output matlab file: outputfile











[montecarlo]
extractions=1000

[montecarlo.JHmaterial]
active=true
value=0.0075
distribution='Normal'

[montecarlo.Jmagnet]
active=true
value=0.028
distribution='Normal'

[montecarlo.Angle]
active=true
value=1
distribution='Uniform'

[montecarlo.Zring]
active=true
value=0.5e-3
distribution='Uniform'

[montecarlo.XMagnets]
active=true
value=0.2e-3
distribution='Uniform'

[montecarlo.YMagnets]
active=true
value=0.2e-3
distribution='Uniform'

[montecarlo.ZMagnets]
active=true
value=0.2e-3
distribution='Uniform'



[montecarlo]
Label to provide data for Montecarlo (MC) simulations.

extractions: Numbers of MC extractions



[montecarlo.JHmaterial]
Label to activate variability of HJ material curves
active: flag for activation
distribution: type of statistical distribution ('Normal' or 'Uniform')
value: Standard deviation of the relative variation of HJ curve (for 'Normal') or limits (+/-) of uniform distribution centered on 1

[montecarlo.Jmagnet]
Label to activate variability of magnetization of each single PM
active: flag for activation
distribution: type of statistical distribution ('Normal' or 'Uniform')
value: Standard deviation of the relative variation of magnetization (for 'Normal') or limits (+/-) of uniform distribution centered on 1


[montecarlo.Angle]
Label to activate variability of orientation angle of each single PM
active: flag for activation
distribution: type of statistical distribution ('Normal' or 'Uniform')
value: Standard deviation of the angle variation in degrees (for 'Normal') or limits (+/-) of uniform distribution centered




