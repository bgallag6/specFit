from setuptools import setup, find_packages

# TODO: Relax allowed versions
install_requires = ["astropy==2.0.3",
                    "matplotlib==1.5.1",
                    "numpy==1.11.2",
                    "scipy==0.18.1",
                    "sunpy==0.8.4",
                    "PyYAML==5.4",
                    "Pillow"]

setup(
    name='specFit',
    version='0.0.1',
    author='Brendan Gallager',
    author_email='bgallag6@masonlive.gmu.edu',
    packages=find_packages(),
    url='http://pypi.python.org/pypi/specFit/',
    license='LICENSE.txt',
    description='A Python library for fitting spectral models to temporal image sequences and visualizing the results',
    long_description=open('README.rst').read(),
    install_requires=install_requires,
    include_package_data=True
)
