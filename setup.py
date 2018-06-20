from setuptools import setup

setup(
    name='ectyper',
    version='0.3.1',
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