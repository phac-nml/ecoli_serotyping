try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'Ecoli Serotyping',
    'author': 'Camille La Rose',
    'url': '',
    'download_url': '',
    'author_email': 'claro100@uottawa.ca',
    'version': '0.1',
    'install_requires': ['nose'],
    'packages': ['main'],
    'scripts': [],
    'name': 'ecoli_serotyping'
}

setup(**config)