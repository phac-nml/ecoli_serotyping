import hashlib
import os.path
import pytest
import random

from ectyper.src.ecprediction import *
from ectyper.src.ecvalidatingfiles import *

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"
REL_DIR = SCRIPT_DIRECTORY + '../temp/'

expected_dict = {'AAJT0200':{
   'gnl|BL_ORD_ID|175 8__wzx__wzx-O148__176 DQ167407.1;O antigen flippase;O148' : ['aee07ddae92b95074a93148a06624a29', 0.998320738875, 0.998320738875, 1191, 'wzx', 'O148'],    'gnl|BL_ORD_ID|165 8__wzx__wzx-O140__166 AB812060.1;O antigen flippase;O140' : ['27d4fd2e1144f746437c3fc4bfcf7ec3', 0.0357686453577, 0.959183673469, 1314, 'wzx', 'O140'],
   'gnl|BL_ORD_ID|28 1__fliC__fliC-H28__29 AY250010.1;flagellin;H28' : ['6dfdab8081e7d834f1b3243985472b2d', 0.987931034483, 0.987931034483, 1740, 'fliC', 'H28'],               'gnl|BL_ORD_ID|31 1__fliC__fliC-H31__32 AY250013.1;flagellin;H31' : ['0ef43266b0331ed986db7371dabdfc98', 0.337529976019, 0.980836236934, 1668, 'fliC', 'H31'],
   'gnl|BL_ORD_ID|10 1__fliC__fliC-H14__11 AY249998.1;flagellin;H14' : ['d4fb30f63599ca22f5ce3c6b9f8028d5', 0.338777979431, 0.962199312715, 1653, 'fliC', 'H14'],               'gnl|BL_ORD_ID|7 1__fliC__fliC-H12__8 AY337471.1;flagellin;H12' : ['b8c9c7ae76ae13eea5287cdd654d3907', 0.292505592841, 0.984934086629, 1788, 'fliC', 'H12'],
   'gnl|BL_ORD_ID|8 1__fliC__fliC-H12__9 AY337474.1;flagellin;H12' : ['ac09ce2077f2225059bc306c8c83318d', 0.291946308725, 0.983050847458, 1788, 'fliC', 'H12'],                 'gnl|BL_ORD_ID|62 1__fliC__fliC-H6__63 AY249991.1;flagellin;H6' : ['a52719aff71d6a30d58d0e6a2023a995', 0.312082574378, 0.975332068311, 1647, 'fliC', 'H6'],
   'gnl|BL_ORD_ID|35 1__fliC__fliC-H34__36 AF34585;flagellin;H34' : ['eadeb19641cb38052520d6c92037790b', 0.31684981685, 0.970093457944, 1638, 'fliC', 'H34'],                   'gnl|BL_ORD_ID|47 1__fliC__fliC-H41__48 AY250020.1;flagellin;H41' : ['1f3cef59dafb2c38024d508ac27439fd', 0.314285714286, 0.963503649635, 1680, 'fliC', 'H41'],
   'gnl|BL_ORD_ID|34 1__fliC__fliC-H34__35 AY250016.1;flagellin;H34' : ['eadeb19641cb38052520d6c92037790b', 0.31684981685, 0.970093457944, 1638, 'fliC', 'H34'],                'gnl|BL_ORD_ID|9 1__fliC__fliC-H12__10 AY249997.1;flagellin;H12' : ['779f01882e3901826062c3d55b8d319a', 0.288590604027, 0.971751412429, 1788, 'fliC', 'H12'],
   'gnl|BL_ORD_ID|0 1__fliC__fliC-H1__1 AB028471.1;flagellin;H1' : ['6143271ae1219405eab366766d4b7593', 0.287472035794, 0.967984934087, 1788, 'fliC', 'H1'],                    'gnl|BL_ORD_ID|50 1__fliC__fliC-H45__51 AY250023.1;flagellin;H45' : ['ad67533120333388945a7712ea6a53b2', 0.300527240773, 0.967924528302, 1707, 'fliC', 'H45'],
   'gnl|BL_ORD_ID|1 1__fliC__fliC-H1__2 L07387.1;flagellin;H1' : ['61dc9c31e66819e245cb7a9ca2e70f1d', 0.286912751678, 0.966101694915, 1788, 'fliC', 'H1'],                      'gnl|BL_ORD_ID|55 1__fliC__fliC-H49__56 AY250026.1;flagellin;H49' : ['8a06559dd8a39a175ed0aaaf58e83c6f', 0.298525073746, 0.971209213052, 1695, 'fliC', 'H49'],
   'gnl|BL_ORD_ID|54 1__fliC__fliC-H49__55 AB028480.1;flagellin;H49' : ['8a06559dd8a39a175ed0aaaf58e83c6f', 0.297997644287, 0.971209213052, 1698, 'fliC', 'H49'],               'gnl|BL_ORD_ID|52 1__fliC__fliC-H46__53 AY250024.1;flagellin;H46' : ['780b60fa959c1550fd8bdb0f39b5f579', 0.499709133217, 0.814990512334, 1719, 'fliC', 'H46'],
   'gnl|BL_ORD_ID|66 1__fliC__fliC-H7__67 AF228494.1;flagellin;H7' : ['3355cec6378330c6530aedc2700af69e', 0.287258248009, 0.951035781544, 1758, 'fliC', 'H7'],                  'gnl|BL_ORD_ID|64 1__fliC__fliC-H7__65 AF228491.1;flagellin;H7' : ['3355cec6378330c6530aedc2700af69e', 0.287258248009, 0.951035781544, 1758, 'fliC', 'H7'],
   'gnl|BL_ORD_ID|75 1__fliC__fliC-H9__76 AY249994.1;flagellin;H9' : ['1489feeac76fb141ac77e2163ecbf1a9', 0.260307998013, 0.935714285714, 2013, 'fliC', 'H9'],                  'gnl|BL_ORD_ID|51 1__fliC__fliC-H46__52 AB028478.1;flagellin;H46' : ['bbcc42bbde8a226400ae3494b900880a', 0.500291885581, 0.813092979127, 1713, 'fliC', 'H46'],
   'gnl|BL_ORD_ID|72 1__fliC__fliC-H7__73 AY249992.1;flagellin;H7' : ['86d2704ec42a318a7b26528261465ff5', 0.285551763367, 0.94538606403, 1758, 'fliC', 'H7'],                   'gnl|BL_ORD_ID|65 1__fliC__fliC-H7__66 AF228492.1;flagellin;H7' : ['86d2704ec42a318a7b26528261465ff5', 0.285551763367, 0.94538606403, 1758, 'fliC', 'H7'],
   'gnl|BL_ORD_ID|11 1__fliC__fliC-H15__12 AY249999.1;flagellin;H15' : ['cb457fe7fb890888e83ca1432af75528', 0.355239786856, 0.884955752212, 1689, 'fliC', 'H15'],               'gnl|BL_ORD_ID|68 1__fliC__fliC-H7__69 AF228496.1;flagellin;H7' : ['c8fd3b49aaadbef6c114577eb35ba872', 0.28441410694, 0.941619585687, 1758, 'fliC', 'H7'],
   'gnl|BL_ORD_ID|60 1__fliC__fliC-H52__61 AY250028.1;flagellin;H52' : ['a04d1a965df7932b0b42fa10764ec3d0', 0.35193452381, 0.9633401222, 1344, 'fliC', 'H52'],                  'gnl|BL_ORD_ID|71 1__fliC__fliC-H7__72 AB334575.1;flagellin;H7' : ['3fb9b7893b507725c49c4be19002b831', 0.283276450512, 0.937853107345, 1758, 'fliC', 'H7'],
   'gnl|BL_ORD_ID|70 1__fliC__fliC-H7__71 AB334574.1;flagellin;H7' : ['3fb9b7893b507725c49c4be19002b831', 0.283276450512, 0.937853107345, 1758, 'fliC', 'H7'],                  'gnl|BL_ORD_ID|69 1__fliC__fliC-H7__70 AY337468.1;flagellin;H7' : ['3fb9b7893b507725c49c4be19002b831', 0.283276450512, 0.937853107345, 1758, 'fliC', 'H7'],
   'gnl|BL_ORD_ID|67 1__fliC__fliC-H7__68 AF228495.1;flagellin;H7' : ['3fb9b7893b507725c49c4be19002b831', 0.283276450512, 0.937853107345, 1758, 'fliC', 'H7'],                  'gnl|BL_ORD_ID|63 1__fliC__fliC-H7__64 AF228487.1;flagellin;H7' : ['3fb9b7893b507725c49c4be19002b831', 0.283276450512, 0.937853107345, 1758, 'fliC', 'H7'],
   'gnl|BL_ORD_ID|26 1__fliC__fliC-H26__27 AY250008.1;flagellin;H26' : ['11bec9ac7863945e285044b85fa328bc', 0.297491039427, 0.936090225564, 1674, 'fliC', 'H26'],               'gnl|BL_ORD_ID|21 1__fliC__fliC-H20__22 AY250003.1;flagellin;H20' : ['c236f038a83412f0a51e6776f6ca2b9b', 0.30849220104, 0.911262798635, 1731, 'fliC', 'H20'],
   'gnl|BL_ORD_ID|16 1__fliC__fliC-H18__17 AY250001.1;flagellin;H18' : ['a16b6e1c69abbb85dcc3c43d5cc28552', 0.297897897898, 0.935849056604, 1665, 'fliC', 'H18'],               'gnl|BL_ORD_ID|18 1__fliC__fliC-H19__19 AY250002.1;flagellin;H19' : ['32eee96c2e861bfd9c1876647b3d5d3e', 0.27904451683, 0.914590747331, 1842, 'fliC', 'H19'],
   'gnl|BL_ORD_ID|17 1__fliC__fliC-H19__18 AY337479.1;flagellin;H19' : ['32eee96c2e861bfd9c1876647b3d5d3e', 0.28041462084, 0.914590747331, 1833, 'fliC', 'H19'],                'gnl|BL_ORD_ID|36 1__fliC__fliC-H37__37 AY250017.1;flagellin;H37' : ['bd80deb2fdbfd9d46ea82408be015218', 0.291221826809, 0.92641509434, 1686, 'fliC', 'H37'],
   'gnl|BL_ORD_ID|23 1__fliC__fliC-H23__24 AY250005.1;flagellin;H23' : ['eff1f961556b4a88b8072371521275c8', 0.276174306735, 0.919020715631, 1767, 'fliC', 'H23'],               'gnl|BL_ORD_ID|30 1__fliC__fliC-H30__31 AY250011.1;flagellin;H30' : ['e02e1735623030f5952d26e0f6d3b90d', 0.270286047869, 0.937246963563, 1713, 'fliC', 'H30'],
   'gnl|BL_ORD_ID|32 1__fliC__fliC-H32__33 AY250014.1;flagellin;H32' : ['ed23dd977c2e0f3cbbacd066e1a3e39c', 0.269118505546, 0.933198380567, 1713, 'fliC', 'H32'],               'gnl|BL_ORD_ID|53 1__fliC__fliC-H48__54 AY250025.1;flagellin;H48' : ['aa7a0638df03c38e5db1dcb09d7fa3a3', 0.299933199733, 0.914460285132, 1497, 'fliC', 'H48'],
   'gnl|BL_ORD_ID|49 1__fliC__fliC-H43__50 AY250022.1;flagellin;H43' : ['7c0cab8d80f814a75e241cd856aa22d7', 0.296812749004, 0.904858299595, 1506, 'fliC', 'H43'],               'gnl|BL_ORD_ID|59 1__fliC__fliC-H51__60 AY250027.1;flagellin;H51' : ['1dddbe105cd8f01595d077fb3e96ab1d', 0.259075907591, 0.887005649718, 1818, 'fliC', 'H51'],
   'gnl|BL_ORD_ID|58 1__fliC__fliC-H51__59 AB028481.1;flagellin;H51' : ['8568a7d4f80bbe9671256f44107a713c', 0.258525852585, 0.885122410546, 1818, 'fliC', 'H51'],               'gnl|BL_ORD_ID|81 3__fllA__fllA-H44__82 #VALUE!;flagellin;H44' : ['50034b78fa6139f4049455b3bf35d759', 0.285010940919, 0.856907894737, 1828, 'fllA', 'H44'],
   'gnl|BL_ORD_ID|3 1__fliC__fliC-H10__4 AY249995.1;flagellin;H10' : ['46d509f1c57c5cb74e9edf1a41760238', 0.349960411718, 0.892929292929, 1263, 'fliC', 'H10'],                 'gnl|BL_ORD_ID|2 1__fliC__fliC-H10__3 AY337482.1;flagellin;H10' : ['46d509f1c57c5cb74e9edf1a41760238', 0.349960411718, 0.892929292929, 1263, 'fliC', 'H10'],
   'gnl|BL_ORD_ID|82 3__fllA__fllA-H55__83 AB269771;flagellin;H55' : ['915d2ce61cdfab6a05054b48c479f306', 0.273694880892, 0.830769230769, 1973, 'fllA', 'H55'],                 'gnl|BL_ORD_ID|57 1__fliC__fliC-H5__58 AY249990.1;flagellin;H5' : ['2097d002ad1b83c687e03d7004d846d6', 0.331045003814, 0.878542510121, 1311, 'fliC', 'H5'],
   'gnl|BL_ORD_ID|61 1__fliC__fliC-H56__62 AY250029.1;flagellin;H56' : ['56cab81dfe4a766531738f8747bb0caa', 0.330282227307, 0.876518218623, 1311, 'fliC', 'H56'],               'gnl|BL_ORD_ID|37 1__fliC__fliC-H38__38 AY337466.1;flagellin;H38' : ['fa132af24eb3d853f0b92d4da688a981', 0.32068452381, 0.877800407332, 1344, 'fliC', 'H38'],
   'gnl|BL_ORD_ID|39 1__fliC__fliC-H38__40 AY250018.1;flagellin;H38' : ['94ab1c21417d8878ae8e3b0de44c88f5', 0.31994047619, 0.875763747454, 1344, 'fliC', 'H38'],                'gnl|BL_ORD_ID|38 1__fliC__fliC-H38__39 AY337473.1;flagellin;H38' : ['94ab1c21417d8878ae8e3b0de44c88f5', 0.31994047619, 0.875763747454, 1344, 'fliC', 'H38'],
   'gnl|BL_ORD_ID|44 1__fliC__fliC-H4__45 AY249989.1;flagellin;H4' : ['cd3369f667072a2d8d54f8562ed29cc3', 0.395238095238, 0.884861407249, 1050, 'fliC', 'H4'],                  'gnl|BL_ORD_ID|42 1__fliC__fliC-H4__43 AJ605764.1;flagellin;H4' : ['cd3369f667072a2d8d54f8562ed29cc3', 0.395238095238, 0.884861407249, 1050, 'fliC', 'H4'],
   'gnl|BL_ORD_ID|41 1__fliC__fliC-H4__42 AJ536600.1;flagellin;H4' : ['cd3369f667072a2d8d54f8562ed29cc3', 0.395238095238, 0.884861407249, 1050, 'fliC', 'H4'],                  'gnl|BL_ORD_ID|56 1__fliC__fliC-H5__57 AY337469.1;flagellin;H5' : ['9e28653f90c9e14104892dd3e94c5456', 0.327993897788, 0.87044534413, 1311, 'fliC', 'H5'],
   'gnl|BL_ORD_ID|43 1__fliC__fliC-H4__44 AJ605765.1;flagellin;H4' : ['f11f18caebb5c00b2a743df77b9c30a6', 0.384761904762, 0.887912087912, 1050, 'fliC', 'H4'],                  'gnl|BL_ORD_ID|33 1__fliC__fliC-H33__34 AY250015.1;flagellin;H33' : ['c2d5dbb952d454e42878c16a49569a4d', 0.341103341103, 0.860784313725, 1287, 'fliC', 'H33'],
   'gnl|BL_ORD_ID|48 1__fliC__fliC-H42__49 AY250021.1;flagellin;H42' : ['d877438aa5b2ce1b83a29ee661c33849', 0.334113973458, 0.866396761134, 1281, 'fliC', 'H42'],               'gnl|BL_ORD_ID|29 1__fliC__fliC-H29__30 AY250012.1;flagellin;H29' : ['d877438aa5b2ce1b83a29ee661c33849', 0.321321321321, 0.866396761134, 1332, 'fliC', 'H29'],
   'gnl|BL_ORD_ID|84 5__flnA__flnA-H17__85 CP002291;flagellin;H17' : ['0ba4eb36a2fcd2d9333cae1936431eec', 0.293307086614, 0.851428571429, 1524, 'flnA', 'H17'],                 'gnl|BL_ORD_ID|40 1__fliC__fliC-H39__41 AY250019.1;flagellin;H39' : ['a45424bf04fc3aa95036e92b25b7a276', 0.337182448037, 0.85546875, 1299, 'fliC', 'H39'],
   'gnl|BL_ORD_ID|45 1__fliC__fliC-H40__46 AJ865464.1;flagellin;H40' : ['d1a780bc0db40b1a4c90c7cb8ee32001', 0.267748478702, 0.796780684105, 1479, 'fliC', 'H40'],               'gnl|BL_ORD_ID|24 1__fliC__fliC-H24__25 AY250006.1;flagellin;H24' : ['5282b0acfb371dc93c2e41fee4a785f1', 0.1431670282, 0.970588235294, 1383, 'fliC', 'H24'],
   'gnl|BL_ORD_ID|73 1__fliC__fliC-H8__74 AJ865465.1;flagellin;H8' : ['9754f3d51f8129ec865f5c6828ab209e', 0.26707234618, 0.79476861167, 1479, 'fliC', 'H8'],                    'gnl|BL_ORD_ID|5 1__fliC__fliC-H11__6 AY337472.1;flagellin;H11' : ['8caff7bb15fd57658b9edebdd18d6c2e', 0.269256987048, 0.79476861167, 1467, 'fliC', 'H11'],
   'gnl|BL_ORD_ID|4 1__fliC__fliC-H11__5 AY337465.1;flagellin;H11' : ['8caff7bb15fd57658b9edebdd18d6c2e', 0.269256987048, 0.79476861167, 1467, 'fliC', 'H11'],                  'gnl|BL_ORD_ID|74 1__fliC__fliC-H8__75 AJ884569.1;flagellin;H8' : ['86b48ff8c3149dc361923c8b89a23563', 0.266396213658, 0.792756539235, 1479, 'fliC', 'H8'],
   'gnl|BL_ORD_ID|46 1__fliC__fliC-H40__47 AJ884568.1;flagellin;H40' : ['655ca62dbbbdaac23f26bc54c1f51659', 0.25016903313, 0.802603036876, 1479, 'fliC', 'H40'],                'gnl|BL_ORD_ID|20 1__fliC__fliC-H2__21 HQ116826;flagellin;H2' : ['8237e98abe46896f7bb6d8ddeb3275e1', 0.269745649264, 0.787109375, 1494, 'fliC', 'H2'],
   'gnl|BL_ORD_ID|19 1__fliC__fliC-H2__20 AF543692;flagellin;H2' : ['c4966dee16a320fc55a3f9619b74dc12', 0.266533066132, 0.785433070866, 1497, 'fliC', 'H2'],                    'gnl|BL_ORD_ID|27 1__fliC__fliC-H27__28 AY250009.1;flagellin;H27' : ['f3939ff5d0e1231ae3d89ed15e3e9d25', 0.242490842491, 0.805352798054, 1365, 'fliC', 'H27'],
   'gnl|BL_ORD_ID|22 1__fliC__fliC-H21__23 AY250004.1;flagellin;H21' : ['8b4b3b29ab37a007ca04a57923520032', 0.234782608696, 0.803970223325, 1380, 'fliC', 'H21'],               'gnl|BL_ORD_ID|6 1__fliC__fliC-H11__7 AY249996.1;flagellin;H11' : ['8b4b3b29ab37a007ca04a57923520032', 0.236842105263, 0.803970223325, 1368, 'fliC', 'H11'],
   'gnl|BL_ORD_ID|14 1__fliC__fliC-H16__15 AY337477.1;flagellin;H16' : ['09064bc24643dc9bb5ffa42624ec78d5', 0.227503168568, 0.785557986871, 1578, 'fliC', 'H16'],               'gnl|BL_ORD_ID|13 1__fliC__fliC-H16__14 AY337476.1;flagellin;H16' : ['09064bc24643dc9bb5ffa42624ec78d5', 0.227503168568, 0.785557986871, 1578, 'fliC', 'H16'] ,
   'gnl|BL_ORD_ID|12 1__fliC__fliC-H16__13 AY337475.1;flagellin;H16' : ['09064bc24643dc9bb5ffa42624ec78d5', 0.227503168568, 0.785557986871, 1578, 'fliC', 'H16'],               'gnl|BL_ORD_ID|15 1__fliC__fliC-H16__16 AY250000.1;flagellin;H16' : ['51f3c3538e86bd56fd16761505da4630', 0.191235059761, 0.793388429752, 1506, 'fliC', 'H16'],
   'gnl|BL_ORD_ID|25 1__fliC__fliC-H25__26 AY250007.1;flagellin;H25' : ['94f0f51796d70c5e0a890beb107786e0', 0.131996658312, 0.934911242604, 1197, 'fliC', 'H25'],               'gnl|BL_ORD_ID|377 9__wzy__wzy-O148__378 DQ167407.1;O antigen polyermase;O148' : ['c10c9586ec7845c902fbbd55f7df991f', 1.0, 1.0, 1143, 'wzy', 'O148'],
   'gnl|BL_ORD_ID|190 8__wzx__wzx-O159__191 EU294176.1;O antigen flippase;O159' : ['8cab52763f6262bcecaa41bece2ae41f', 0.0304347826087, 1.0, 1380, 'wzx', 'O159'] }}


def test_filterPredictions():

     test_prediction = parseResults([REL_DIR + 'NC_011749.1.xml',REL_DIR + 'AAJT02.1.xml'])
     i = 0

     while i<4:
        id_num = random.random() * 100
        length_num = random.random() * 100
        exp_dict = {}
        length_exp_dict = 0

        for expected_genome, expected_value in expected_dict.iteritems():
            if isinstance(expected_value, dict):
                expected_temp = {}
                for expected_title, expected_list in expected_value.iteritems():
                    expected_id_num = expected_list[2] * 100
                    expected_length_num = expected_list[1] * 100
                    if expected_length_num >= length_num:
                        if expected_id_num >=id_num:
                            expected_temp[expected_title] = [expected_list[0], expected_list[3]]
                            exp_dict[expected_genome]= expected_temp
                            length_exp_dict +=1
                        else:
                            expected_temp[expected_title] = [expected_list[0], 'NA']
                            exp_dict[expected_genome] = expected_temp
                            length_exp_dict +=1
                    else:
                        expected_temp[expected_title] = [expected_list[0], 'NA']
                        exp_dict[expected_genome] = expected_temp
                        length_exp_dict +=1
            else:
                exp_dict[expected_genome] = 'NA'
                length_exp_dict +=1

        test_pred = filterPredictions(test_prediction, id_num, length_num)

        length_test_pred = 0

        for test_genome, test_value in test_pred.iteritems():
            if test_value != 'NA':
                for test_title, test_hsp in test_value.iteritems():
                    hash_hsp = hashlib.md5()
                    hash_hsp.update(str(test_hsp[0]))
                    if exp_dict[test_genome][test_title][0] != str(hash_hsp.hexdigest()):
                        pytest.fail("The dictionaries aren't the same. Test failed. \nCAUSE \nAlignment title: " + str(test_title) + "\nHSP:\n" + str(test_hsp[0]))
                    if exp_dict[test_genome][test_title][1] != test_hsp[1]:
                        pytest.fail("The dictionaries aren't the same. Test failed. \nCAUSE \nAlignment title: " + str(test_title) + "\nLength:\n" + str(test_hsp[1]))

                    length_test_pred += 1
            else:
                if exp_dict[test_genome] != 'NA':
                    pytest.fail("The dictionaries aren't the same. Test failed. \nCAUSE \nGenome: \n" + test_genome)
                    length_test_pred +=1

        if length_exp_dict != length_test_pred:
            pytest.fail("The resulting dictionary is not the same length as the expected dictionary.")
        i+=1


def test_getProductPercentage():

    test_prediction = parseResults([REL_DIR + 'NC_011749.1.xml',REL_DIR + 'AAJT02.1.xml'])
    exp_dict = {}

    for expected_genome, expected_value in expected_dict.iteritems():
        expected_temp = {}
        for expected_title, expected_list in expected_value.iteritems():
            expected_id_num = expected_list[2]
            expected_length_num = expected_list[1]
            exp_product = expected_id_num * expected_length_num
            expected_temp[expected_title] = "{0:.4f}".format(exp_product)
        exp_dict[expected_genome]  = expected_temp

    for test_genome, test_value in test_prediction.iteritems():
        if test_value != 'NA':
            for test_title, test_info in test_value.iteritems():
                test_length = test_info.keys()[0]
                test_hsp = test_info[test_length]
                test_product = getProductPercentage(test_length, test_hsp)
                test_product = "{0:.4f}".format(test_product)

                assert getProductPercentage(1, test_hsp) == 0

                assert getProductPercentage(0, test_hsp) == 0

                assert test_product == exp_dict[test_genome][test_title]


def test_sortMatches():

    test_sortmatches = sortMatches(filterPredictions(parseResults([REL_DIR + 'AAJT02.1.xml']), 0, 0))

    for test_genome, test_typedict in test_sortmatches.iteritems():
        for test_type, test_infolist in test_typedict.iteritems():
            for test_info in test_infolist:
                test_title = test_info.get('title')
                assert expected_dict[test_genome][test_title][4] == test_type

    assert sortMatches({'NC_011749': 'NA'}) == {'NC_011749': 'NA'}


def test_searchType():

    for test_genome, test_value in expected_dict.iteritems():
        for test_title, test_info in test_value.iteritems():
            assert searchType(test_title, test_info[4]) == test_info[5]

    assert searchType('fli A456', 'wzx') == None
    assert searchType(' ', 'fliC') == None
    assert searchType('fliC-H15', ' ') == None


def test_findTopMatches():

    exp_top_match1 = {'AAJT0200':{
        'gnl|BL_ORD_ID|377 9__wzy__wzy-O148__378 DQ167407.1;O antigen polyermase;O148' : ['c10c9586ec7845c902fbbd55f7df991f', 1.0, 1.0, 1143, 'wzy', 'O148'],
        'gnl|BL_ORD_ID|28 1__fliC__fliC-H28__29 AY250010.1;flagellin;H28' : ['6dfdab8081e7d834f1b3243985472b2d', 0.987931034483, 0.987931034483, 1740, 'fliC', 'H28']
        }}

    exp_top_match_2 = {'GENOME_1': 'NA', 'GENOME_2':  {'wzy': [{'title': 'gnl|BL_ORD_ID|377 9__wzy__wzy-O148__378 DQ167407.1;O antigen polyermase;O148', 'hsp': 'example_hsp', 'length': '1166', 'perc': 0.964}]},
                       'GENOME_3': {'wzx':[{'title': 'gnl|BL_ORD_ID|377 9__wzx__wzx-O148__378 DQ167407.1;O antigen polyermase;O148', 'hsp': 'example_hsp', 'length': 'NA', 'perc': 0}]}}

    test_prediction = (filterPredictions(parseResults([REL_DIR + 'AAJT02.1.xml']), 0, 0))
    test_top_matches = findTopMatches(sortMatches(test_prediction))

    for test_genome, test_value in test_top_matches.iteritems():
        test_otype = test_value.get('otype')
        test_htype = test_value.get('htype')
        otype_title = test_otype['title']
        htype_title = test_htype['title']

        assert exp_top_match1[test_genome][otype_title][5] == 'O148'
        assert exp_top_match1[test_genome][htype_title][5] == 'H28'

    assert findTopMatches(exp_top_match_2) == {
        'GENOME_1': {'otype': 'NA', 'htype': 'NA', 'predictionstrength': 'NA'},
        'GENOME_2': {'otype': {'title': 'gnl|BL_ORD_ID|377 9__wzy__wzy-O148__378 DQ167407.1;O antigen polyermase;O148', 'hsp': 'example_hsp', 'length': '1166', 'perc': 0.964}, 'htype': 'NA', 'predictionstrength': 'Top match'},
        'GENOME_3': {'otype': 'NA', 'htype': 'NA', 'predictionstrength': 'Top match'}
    }