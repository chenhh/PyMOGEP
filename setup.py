#-*-coding:utf-8-*-
'''
@author: Hung-Hsin Chen
@mail: chenhh@par.cse.nsysu.edu.tw
'''
from distutils.core import setup 

setup(
    name = 'PyMOGEP', 
    version = '0.1.0',
    license = 'GPLv2',
    py_modules = ['PyMOGEP'], 
    author = 'Hung-Hsin Chen',
    author_email = 'chenhh@par.cse.nsysu.edu.tw', 
    url = 'https://github.com/chenhh/PyMOGEP', 
    description = 'multi-objective gene expression programming',
    install_requires = ['numpy', 'pandas'],
    
    package_dir = {'': 'src'},
)
