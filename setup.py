from setuptools import setup


def myversion():
    from setuptools_scm.version import get_local_dirty_tag
    def clean_scheme(version):
        return get_local_dirty_tag(version) if version.dirty else '+clean'

    return {'local_scheme': clean_scheme}


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