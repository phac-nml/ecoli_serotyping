from setuptools import setup
from setuptools_scm import get_version

def theversion():
    get_version(root='.', relative_to=__file__)


setup(
    name='ectyper',
    version=theversion(),
    description='E. coli serotyping',
    url='https://github.com/phac-nml/ecoli_serotyping',
    author='Chad Laing, Sam Sung, Camille La Rose',
    author_email='chad.laing@canada.ca, sam.sung@canada.ca, claro100@uottawa.ca',
    license='Apache 2',
    scripts=['bin/ectyper'],
    packages=['ectyper'],
    package_data={'ectyper': ['Data/*']},
    zip_safe=False,
    test_suite='py.test'
)