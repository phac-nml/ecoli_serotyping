from setuptools import setup

setup(
    use_scm_version=myversion(),
    setup_requires=['setuptools_scm', 'setuptools_scm_git_archive'],
    name='ectyper',
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