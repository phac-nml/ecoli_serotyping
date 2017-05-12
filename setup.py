try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'E. coli serotyping and virulence finding',
    'author': 'Camille La Rose, Chad Laing',
    'url': 'https://github.com/phac-nml/ecoli_serotyping',
    'download_url': '',
    'author_email': 'claro100@uottawa.ca, chad.laing@canada.ca',
    'version': '0.1',
    'install_requires': ['nose'],
    'packages': ['main','src'],
    'scripts': [],
    'name': 'ecoli_serotyping'
}

setup(**config)