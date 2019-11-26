from setuptools import setup
import os

def read(filename):
    return open(os.path.join(os.path.dirname(__file__), filename)).read()

setup(
    name = 'COMETSC',
    version = '0.1.12',
    long_description=read('README.md'),
    packages = ['Comet'],
    install_requires=[
        'adjustText==0.7.3' ,
        'atomicwrites==1.1.5' ,
        'attrs==18.1.0',
        'chardet==3.0.4',
        'cycler==0.10.0',
        'Cython==0.28.4',
        'decorator==4.3.0',
        'future==0.16.0',
        'idna==2.7',
        'ipython-genutils==0.2.0',
        'jsonschema==2.6.0',
        'jupyter-core==4.4.0',
        'kiwisolver==1.0.1',
        'matplotlib==3.0.0',
        'more-itertools==4.2.0',
        'nbformat==4.4.0',
        'numpy==1.14.5',
        'pandas==0.23.3',
        'plotly==2.7.0',
        'pluggy==0.6.0',
        'py==1.5.4',
        'pyparsing==2.2.0',
        'pytest==3.6.3',
        'python-dateutil==2.7.3',
        'pytz==2018.5',
        'requests==2.20.1',
        'scipy==1.1.0',
        'scikit-learn==0.21.0',
        'six==1.11.0',
        'traitlets==4.3.2',
        'urllib3==1.24.3',
        'xlmhg==2.4.9'],
    include_package_data=True,
    entry_points = {
        'console_scripts': [
            'Comet = Comet.__main__:main'
            ]
        }
    
)
