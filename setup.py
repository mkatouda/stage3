from glob import glob
from os.path import basename, splitext
from setuptools import setup

setup(
    name="stage3",
    version="0.1.0",
    install_requires=[],
    entry_points={
        'console_scripts': [
            'stage3=stage3.stage3:main',
            'charmm2gromacs=stage3.charmm2gromacs:main',
            'xvg2yavg=stage3.xvg2yavg:main'
        ],
    },
    include_package_data=True,
    author="Michio Katouda",
    author_email="katouda@rist.or.jp",
    description="Automatic GROMACS Topology Generation tool of organic molecules using the GAFF, OPLS-AA, and CGenFF force fields",
    url="https://github.com/mkatouda/stage3",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL-3.0 license",
    ],
    python_requires='>=3.7',
)
