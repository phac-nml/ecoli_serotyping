from setuptools import setup
from ectyper import __version__

setup(
    name='ectyper',
    version=__version__,
    description='E. coli serotyping',
    url='https://github.com/phac-nml/ecoli_serotyping',
    author='Chad Laing,Kyrylo Bessonov, Sam Sung, Camille La Rose, ',
    author_email='chad.laing@canada.ca, kyrylo.bessonov@canada.ca, sam.sung@canada.ca, claro100@uottawa.ca',
    license='Apache 2',
    scripts=['bin/ectyper'],
    packages=['ectyper'],
    package_data={'ectyper': ['Data/*']},
    zip_safe=False,
    test_suite='py.test'
)