from setuptools import setup
from ectyper import __version__


setup(
    name='ectyper',
    version=__version__,
    description='Escherichia coli fast serotyping using both raw reads and assemblies with automatic species identification',
    url='https://github.com/phac-nml/ecoli_serotyping',
    author='Chad Laing, Kyrylo Bessonov, Sam Sung, Camille La Rose, ',
    author_email='chad.laing@canada.ca, kyrylo.bessonov@canada.ca, sam.sung@canada.ca, claro100@uottawa.ca',
    license='Apache 2',
    scripts=['bin/ectyper'],
    packages=['ectyper'],
    install_requires=['requests','biopython','pandas'],
    package_data={'ectyper': ['Data/*.json', 'Data/*.py']},
    zip_safe=False,
    test_suite='py.test'
)

