from setuptools import setup, find_packages

setup(
    name='litho2strat',
    version='1.0.0',
    description='Stratigraphy estimation from drillhole lithology data',
    url='https://github.com/Loop3D/litho2strat',
    author='Vitaliy Ogarko',
    author_email='vogarko@gmail.com',
    packages=find_packages(),
    py_modules=['litho2strat'],
    install_requires=['numpy', 'matplotlib', 'networkx'],
)