from setuptools import setup, find_packages

setup(
    name="ECCO",
    description="Python functionality of the ECCO project: ",
    long_description="Effects of Climate Change on Boreal Lake Ecosystems:"+/
    "Productivity and Community Responses Python software, to extract"+/
    "data from CORDEX climate models and output timeseries per lake in hdf5.",
    url="https://github.com/benlaken/ECCO",
    download_url="https://github.com/benlaken/ECCO",
    classifiers=['UiO', 'ECCO'],
    author="Benjamin Laken",
    author_email="benlaken@gmail.com",
    version="0.1",
    license="CC-BY-4.0",
    packages=find_packages(exclude=['*test']),
    install_requires=['argparse', 'time', 'NetCDF4', 'basemap',
                      'geos', 'simplejson', 'json', 'h5py', 'pandas']
)
