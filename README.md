# PyKMOD
A Kurucz model atmosphere interpolator for Python! (Now with PHOENIX atmospheres!)

This is essentially a re-writing of the original tool, [KMOD](http://hebe.as.utexas.edu/stools/), which was written in IDL, but modified such that the output files are immediately readable by the popular spectroscopic analysis program [MOOG](http://www.as.utexas.edu/~chris/moog.html).
Python being a more flexible language, PyKMOD can be run either by importing it into a Python script and using the function pykmod(), or it can be run directly from the command line via your Python interpreter.

Included for convenience are ATLAS9 and PHOENIX BT-Settl model atmospheres which are ready to be used with PyKMOD. The ATLAS9 atmospheres were downloaded from [this page](http://research.iac.es/proyecto/ATLAS-APOGEE/) on the ATLAS-APOGEE website ([Meszaros et al. 2012](https://ui.adsabs.harvard.edu/abs/2012AJ....144..120M/abstract)). The PHOENIX atmospheres were extracted from the PHOENIX files [here](https://phoenix.ens-lyon.fr/Grids/BT-Settl/GNS93/STRUCTURES/), described in [Allard 2016](https://ui.adsabs.harvard.edu/abs/2016sf2a.conf..223A/abstract). PyKMOD can be adapted for reading any model atmosphere file by modifying its rd_kmod function.

(PyKMOD should be run with the latest versions of SciPy and Numpy)

### Installation
- Download the files provided and unzip them to a directory of your choice. 
- Extract the ATLAS and PHOENIX archives so that they are folders in the same directory as the Python scripts.
- If you want to import it as a python module, you can also copy the fully extracted files to your Python lib folder to make this easier.

### Using PyKMOD
#### Command Line
The command line syntax is as follows:
```
python pykmod.py (Teff (K)) (logg ([cm/s^2])) (vmicro (km/s)) ([Fe/H] (dex)) [output file path]
```

```
python phxmod.py (Teff (K)) (logg ([cm/s^2])) (vmicro (km/s)) ([Fe/H] (dex)) [output file path]
```
The output file path is optional. If none is provided, it will default to a file named "modelatmosphere.txt" in the PyKMOD directory.

#### As Python Module
```python
from pyKMOD.pykmod import pykmod

pykmod(Teff,logg,vmicro,FeH,outfile)
```

```python
from pyKMOD.phxmod import phxmod

phxmod(Teff,logg,vmicro,FeH,outfile)
```
In this case an output file path is required

