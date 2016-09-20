[2016/09/12 - Camille La Rose]
- finished writing tests for getListGenomes and checkFiles methods
- added the EcOH.fasta DB
- wrote the method to generate the EcOH.fasta DB
- minor directory fixes

[2016/09/13 - Camille La Rose]
- wrote the methods runBlastQuery and parseResults
- started the method getGenomeName to identify the genomes in the GENOMES dictionary

[2016/09/14 - Camille La Rose]
- documentation adds + fixes
- adjusted the methods runBlastQuery and parseResults (lists to dictionaries)
- wrote the tests for initializeDB and getGenomeName

[2016/09/15 - Camille La Rose]
- finished writing tests for parseResults and  findPerfectMatches
- created  new module ecprediction.py for prediction methods

[2016/09/16 - Camille La Rose]
- created new module ecvalidatingfiles.py for validation methods
- fixed tests for parseResults and findPerfectMatches
- fixed methods parseResults and findPerfectMatches so that it returns the GENOME dict.
- added filterPredictions method

[2016/09/19 - Camille La Rose]
- created a new test module (ecprediction_tests) to separate the tests by module
- created new method findTopMatch
- finished the tests for findPerfectMatches and filterPredictions

[2016/09/20 - Camille La Rose]
- removed findPerfectMatches and its test method
- finished findTopMatch