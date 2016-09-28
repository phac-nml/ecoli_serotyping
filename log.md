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

[2016/09/21 - Camille La Rose]
- fixed findTopMatch and findTopMatches
- started tests for the above, as well as for searchType and getProductPercentage

[2016/09/22 - Camille La Rose]
- fixed filterPredictions to include NAs
- fixed test for filterPredictions
- still writing tests for findTopMatch, findTopMatches, searchType and getProductPercentage

[2016/09/23 - Camille La Rose]
- finished tests for findTopMatch, findTopMatches, searchType and getProductPercentage
- added documentation for the above methods
- added logging
- will add the following to the DB:
      H:
        fliC_54_AJ605766_H17
        fliC_50_AJ515904_H17
        fliC_888881_AVRH01000047_H52
        fliC_274_AY337470_H42
        fliC_311_AGSG01000116_H25
        fliC_111197_JH965342_H29
        fliC_286_AY337485_H29
        fliC_109_AM231154_H27
        fliC_2_AIEY01000041_H6
        fliC_84_AIHL01000060_H21
        fliC_82_AIFA01000047_H6
        fliC_285_AY337483_H26
        fliC_314_CP000247_H31
        fliC_92_JH954529_H16
        fliC_93_JH953794_H16
        fliC_300_AB128919_H16
        fliC_21_AB028474_H7
        fliC_95_JH694260_H7
        fliC_97_KB000721_H7
        fliC_99_KB007180_H7
        fliC_100_AOES01000098_H7
        fliC_108_KB006714_H7
        fliC_101_JH946604_H7
        fliC_105_AAJT02000052_H28
        fliC_83_AIFX01000055_H12
        fllA_1_AB269771_H55
        fllA_2_AB269770_H44
      O:
        wzx_36_DQ462205-FJ539194_O28ac/O42
        wzx_86_289152760_O83
        wzx_103_KB021482_O104
        wzx_105_AFPS01000083_O104
        wzx_106_AFOB02000091_O104
        wzx_126_DQ676933_O123
        wzx_153_AAJT02000037_O148
        wzx_168_EU296420_O164
        wzx_178_DQ008592_O174
        wzx_194_NC:013353_O103
        wzx_199_AKMA01000036_O157
        wzx_200_AKLI01000048_O157
        wzy_121_AB812063_O153/O178
        wzy_152_AIGX01000028_O45
        wzy_166_AJ426423_O6
        wzy_171_CP003034_O7
        wzy_193_JH964427_O157
        wzy_195_JH970567_O157
        wzy_199_AKLQ01000042_O157
        wzy_200_AKLV01000031_O157
        wzy_201_JH953200_O157
        wzy_202_EF027120_O103
        wzy_205_AP010958_O103
        wzy_206_EF027115_O103
        wzy_215_AY863412_O145
        wzx_406_56-54Cigleris_O128ab
        wzy_406_E68_O141ab
        wzy_407_RVC2907_O141ac
        wzx_211_AB812082_O186
        wzy_224_AB972416_O77
        wzt_3_AB010293_O9
        wzt_7_D43637_O9
        wzm_3_AB010293_O9
        wzm_8_D43637_O9
        wzx_400_KP835694_O125ab
        wzx_405_5564-64_O128ac
        wzy_400_KP835694_O125ab
        wzy_405_5564-64_O128ac

[2016/09/26 - Camille La Rose]
- created script to add the above sequences to the database (as well as format their titles)
- documentation

[2016/09/27 - Camille La Rose]
- added test script to compare with the Danish database