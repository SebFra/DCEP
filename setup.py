from setuptools import setup, find_packages

setup(
    name='DCEP',
    version='0.1',
    packages=['DcepGenerator'],
    keywords='generator dcep ctw2018',
    url='https://github.com/SebFra/DCEP.git',
    scripts=['DcepGenerator/ProblemGenerator.py', 'DcepGenerator/UnitigsGenerator.py'],
    license='GPLv3',
    author='Sebastien François',
    author_email='sebastien.francois@inria.fr',
    description='Simulator of instances for the DCEP problem. (CTW2018)',
    platforms='ALL',
    install_requires= ['matplotlib>=2.0.0', 'networkx>=1.11']
)
