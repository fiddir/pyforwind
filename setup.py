import setuptools

with open("README.rst", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyforwind",
    version="0.1.2",
    description="Synthetic IEC-conform wind fields with extended turbulence characteristics",
    long_description= long_description,
    maintainer = "Jan Friedrich",
    license="LGPL-3.0",
    keywords="synthetic wind fields, inflow turbulence, wind energy",
    url="https://github.com/fiddir/pyforwind",
    packages=setuptools.find_packages(),
    install_requires=['requests','scipy','numpy','matplotlib', 'pandas'],
    platforms=['any'],
    classifiers=[],
    entry_points={"console_scripts": ["pyforwind = pyforwind.pyforwind:main"]},
)