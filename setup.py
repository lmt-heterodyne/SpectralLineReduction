from setuptools import setup, find_packages, Extension
import re

modcontents = open('lmtslr/__init__.py').read()
version = re.search(r"__version__ = '([^']*)'",modcontents).group(1)

setup(
    name = "SpectralLineReduction",
    version = version,
    description = "SpectralLineReduction: Package for spectral line reduction for LMT",
    author = "Gopal Narayanan",
    author_email = "gopal@astro.umass.edu",
    packages = find_packages(),
    scripts = ['bin/grid_data.py', 'bin/multiprocess_otf_map.py', 'bin/process_map.py',
               'bin/process_otf_map.py', 'bin/process_otf_map2.py',
               'bin/process_ps.py', 'bin/view_cube.py',
               'bin/lmtinfo.py', 'bin/make_spec_fits.py',
               'bin/view_spec_file.py',  'bin/view_spec_point.py']
    )
