from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='ectyper',
    version='0.0.9',
    description='E. coli serotyping',
    url='https://github.com/phac-nml/ecoli_serotyping',
    author='Camille La Rose, Chad Laing, Sam Sung',
    author_email='claro100@uottawa.ca, chad.laing@canada.ca, sam.sung@canada.ca',
    license='MIT',
    scripts=['bin/ectyper'],
    packages=['ectyper'],
    package_data={'ectyper': ['Data/*']},
    install_requires=[
          'future',
      ],
    zip_safe=False,
    test_suite='nose.collector'
)
