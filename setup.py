from setuptools import find_packages
from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration
def configuration(parent_package='', top_path=None):
    config = Configuration('fibos', parent_package, top_path)
    config.add_extension('main75', sources=['fibos/fortran/main75.pyf','fibos/fortran/main75.f'])
    config.add_extension('ds75', sources=['fibos/fortran/ds75.pyf', 'fibos/fortran/ds75.f'])
    config.add_extension('occsurf75', sources=['fibos/fortran/occsurf75.pyf', 'fibos/fortran/occsurf75.f'])
    config.add_extension('renum75', sources=['fibos/fortran/renum75.pyf', 'fibos/fortran/renum75.f'])
    config.add_extension('surfcal76', sources=['fibos/fortran/surfcal76.pyf', 'fibos/fortran/surfcal76.f'])
    config.add_extension('respak75', sources=['fibos/fortran/respak75.pyf','fibos/fortran/respak75.f'])
    return config

if __name__ == '__main__':
    setup(
        name = 'fibos',
        version = '1.0.1',
        install_requires = [
            'setuptools>=67.8.0',
            'numpy>=1.21.5',
            'biopython>=1.81',
            'testresources>=1.8.0',
            'pandas>=2.0.0',
        ],
        description = 'Package to calcule Occluded Surfaces',
        author = 'Herson Hebert Mendes Soares',
        packages=find_packages(),
        data_files=[('fibos', ['fibos/radii'])],
        ext_modules=[ext for ext in configuration(top_path='').todict()['ext_modules']],
    )
