from setuptools import setup

setup(
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    name='ectyper',
    description='E. coli serotyping',
    url='https://github.com/phac-nml/ecoli_serotyping',
    author='Camille La Rose, Chad Laing, Sam Sung',
    author_email='claro100@uottawa.ca, chad.laing@canada.ca, sam.sung@canada.ca',
    license='Apache 2',
    scripts=['bin/ectyper'],
    packages=['ectyper'],
    package_data={'ectyper': ['Data/*']},
    zip_safe=False,
    test_suite='py.test'
)