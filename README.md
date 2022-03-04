# SpectralLineReduction

Spectral Line Reduction software for SEQUOIA OTF mapping.

# LMTOY

Working within the LMTOY environment is the preferred method of using
"lmtslr", together with "dreampy3" and other LMT spectral line
reduction tools.  This module can of course be used standalone, but
tools like the pipeline and common tools are maintained in
LMTOY. There are also installation instructions that incoorporate
those listed below here.

See https://github.com/astroumd/lmtoy where we also maintain all [issues](https://github.com/astroumd/lmtoy/issues)

# Normal startup after cloning the distribution

$ virtualenv -p python3 venv
$ source venv/bin/activate
$ pip install -r requirements.txt

# Ipython shell

$ ipython profile create lmtslr
$ ipython --profile=lmtslr

# A test PS reduction

$ ipython --profile=lmtslr

```python
from lmtslr.reduction.line_reduction import read_obsnum_ps, LineData
import matplotlib.pyplot as pl
obsnum = 82480
bank = 0
tsys = 180
list_of_pixels = [0, 2]
use_calibration = True
data_path = '/home/gopal/engineering/lmt/spectral_line_software/SpectralLineReduction/example_data/'
I, S = read_obsnum_ps(obsnum, list_of_pixels, bank, use_calibration, tsys=tsys, path=data_path)
LD_0 = LineData(I, bank, S.nchan, S.bandwidth, S.roach[0].ps_spectrum)
LD_0.x_vsrc()
pl.plot(LD_0.xarray,LD_0.yarray)
pl.axis([-20, 20, -1, 1])
pl.xlabel(LD_0.xname)
```

# Windows Installation and Setup

After cloning the distribution, ensure you have virtualenv installed.  If not, install it via "pip install virtualenv".  Afterwards, change directory to the local installation of the distribution and initialize the virtualenv using virtualenv -p 3 venv.  

Change directory to venv/Scripts and run activate.bat. 

Return to the original SpectralLineReduction directory and run the command "pip install -r requirements.txt"

To test that setup was successful, run the test PS reduction above.

