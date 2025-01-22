from setuptools import setup
from ectyper import __version__


setup(
    name='ectyper',
    version=__version__,
    description='Escherichia coli fast serotyping using both raw reads and assemblies with automatic species identification',
    url='https://github.com/phac-nml/ecoli_serotyping',
    author='Kyrylo Bessonov, Chad Laing',
    author_email='chad.laing@canada.ca, kyrylo.bessonov@canada.ca',
    license='Apache 2',
    scripts=['bin/ectyper'],
    packages=['ectyper'],
    install_requires=['requests','biopython<1.85','pandas<3'],
    package_data={'ectyper': ['Data/*.json', 'Data/*.py']},
    zip_safe=False,
    test_suite='py.test',
    entry_points={
        'console_scripts': [
            'ectyper_init=ectyper.init:main',
            'ectyper=ectyper.ectyper:run_program'
        ],
    }
)

