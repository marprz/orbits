function p = sat1d2
p=zeros(288, 3 );
p = [ 17457.5 -19827.1 3459.82;
18026.7 -19459.4 2519.71;
18561.6 -19054.7 1574.81;
19061.3 -18614 626.893;
19524.8 -18137.9 -322.232;
19951.2 -17627.6 -1270.77;
20339.8 -17083.8 -2216.91;
20689.7 -16507.7 -3158.87;
21000.5 -15900.4 -4094.87;
21271.4 -15263 -5023.13;
21502 -14596.7 -5941.88;
21691.9 -13902.8 -6849.41;
21840.6 -13182.6 -7743.97;
21947.9 -12437.4 -8623.88;
22013.6 -11668.8 -9487.48;
22037.6 -10878.1 -10333.1;
22019.9 -10066.8 -11159.2;
21960.4 -9236.56 -11964.2;
21859.3 -8388.83 -12746.5;
21716.8 -7525.26 -13504.7;
21533.2 -6647.46 -14237.4;
21308.7 -5757.12 -14943.1;
21043.9 -4855.91 -15620.5;
20739.1 -3945.53 -16268.3;
20395.1 -3027.72 -16885.3;
20012.4 -2104.2 -17470.3;
19591.7 -1176.72 -18022.1;
19133.9 -247.052 -18539.9;
18639.7 683.056 -19022.5;
18110.3 1611.84 -19469;
17546.4 2537.53 -19878.5;
16949.3 3458.38 -20250.4;
16319.9 4372.64 -20583.8;
15659.6 5278.57 -20878.1;
14969.5 6174.45 -21132.8;
14251 7058.58 -21347.3;
13505.4 7929.26 -21521.2;
12734.1 8784.85 -21654.3;
11938.6 9623.71 -21746.1;
11120.4 10444.2 -21796.5;
10281 11244.9 -21805.4;
9422.04 12024.1 -21772.8;
8545.14 12780.4 -21698.7;
7651.96 13512.3 -21583.3;
6744.21 14218.4 -21426.6;
5823.6 14897.4 -21229.1;
4891.88 15548 -20991.1;
3950.85 16168.8 -20712.9;
3002.28 16758.8 -20395.2;
2047.99 17316.7 -20038.5;
1089.8 17841.6 -19643.4;
129.548 18332.2 -19210.8;
-830.93 18787.8 -18741.3;
-1789.8 19207.5 -18236;
-2745.21 19590.3 -17695.7;
-3695.34 19935.6 -17121.5;
-4638.35 20242.6 -16514.4;
-5572.44 20510.9 -15875.6;
-6495.81 20739.7 -15206.3;
-7406.68 20928.8 -14507.8;
-8303.3 21077.7 -13781.4;
-9183.92 21186 -13028.5;
-10046.9 21253.6 -12250.6;
-10890.4 21280.4 -11449.1;
-11713 21266.2 -10625.5;
-12513.1 21211.1 -9781.39;
-13288.9 21115.2 -8918.47;
-14039.2 20978.6 -8038.35;
-14762.4 20801.5 -7142.72;
-15457 20584.4 -6233.3;
-16121.9 20327.6 -5311.85;
-16755.6 20031.6 -4380.13;
-17356.9 19697 -3439.94;
-17924.7 19324.4 -2493.1;
-18457.9 18914.4 -1541.43;
-18955.4 18468 -586.766;
-19416.2 17985.9 369.042;
-19839.5 17469.1 1324.15;
-20224.5 16918.5 2276.71;
-20570.3 16335.3 3224.89;
-20876.3 15720.5 4166.84;
-21142 15075.3 5100.74;
-21366.7 14401 6024.8;
-21550.2 13698.9 6937.21;
-21691.9 12970.3 7836.22;
-21791.6 12216.6 8720.08;
-21849.2 11439.3 9587.09;
-21864.5 10640 10435.6;
-21837.6 9820.05 11263.9;
-21768.3 8981.16 12070.4;
-21657 8124.92 12853.6;
-21503.8 7252.99 13612;
-21309 6367.06 14344;
-21073 5468.84 15048.4;
-20796.3 4560.07 15723.6;
-20479.3 3642.51 16368.5;
-20122.8 2717.92 16981.8;
-19727.4 1788.1 17562.3;
-19294 854.843 18108.9;
-18823.2 -80.0499 18620.5;
-18316.1 -1014.77 19096.2;
-17773.7 -1947.53 19535.1;
-17197 -2876.51 19936.3;
-16587.1 -3799.94 20299;
-15945.3 -4716.03 20622.7;
-15272.7 -5623.02 20906.6;
-14570.6 -6519.17 21150.2;
-13840.5 -7402.76 21353.2;
-13083.8 -8272.09 21515.1;
-12301.9 -9125.49 21635.5;
-11496.3 -9961.34 21714.4;
-10668.6 -10778 21751.6;
-9820.37 -11574 21746.9;
-8953.28 -12347.7 21700.5;
-8068.99 -13097.8 21612.5;
-7169.2 -13822.6 21483;
-6255.65 -14521 21312.3;
-5330.09 -15191.5 21100.8;
-4394.29 -15832.9 20848.8;
-3450.06 -16444 20556.9;
-2499.21 -17023.5 20225.7;
-1543.55 -17570.5 19855.8;
-584.919 -18083.9 19448;
374.854 -18562.7 19002.9;
1333.93 -19006 18521.6;
2290.49 -19413 18005;
3242.71 -19782.9 17453.9;
4188.76 -20115.1 16869.6;
5126.85 -20408.9 16253.1;
6055.2 -20663.8 15605.7;
6972.04 -20879.3 14928.5;
7875.63 -21055 14222.9;
8764.26 -21190.5 13490.3;
9636.25 -21285.8 12731.9;
10489.9 -21340.5 11949.4;
11323.7 -21354.7 11144.2;
12136 -21328.2 10317.8;
12925.3 -21261.2 9471.89;
13690 -21153.7 8607.95;
14428.8 -21006.1 7727.69;
15140.3 -20818.6 6832.77;
15823 -20591.6 5924.91;
16475.8 -20325.4 5005.84;
17097.4 -20020.7 4077.28;
17686.6 -19678 3141.02;
18242.3 -19298 2198.83;
18763.5 -18881.4 1252.48;
19249.2 -18428.9 303.787;
19698.5 -17941.6 -645.464;
20110.5 -17420.2 -1593.47;
20484.4 -16865.8 -2538.44;
20819.7 -16279.5 -3478.59;
21115.5 -15662.3 -4412.13;
21371.4 -15015.4 -5337.3;
21586.9 -14340 -6252.35;
21761.6 -13637.5 -7155.55;
21895.2 -12909.1 -8045.18;
21987.3 -12156.3 -8919.57;
22037.9 -11380.5 -9777.06;
22046.8 -10583 -10616;
22014 -9765.51 -11434.9;
21939.6 -8929.48 -12232.1;
21823.7 -8076.52 -13006.1;
21666.5 -7208.23 -13755.5;
21468.3 -6326.26 -14478.9;
21229.5 -5432.27 -15174.8;
20950.5 -4527.97 -15841.9;
20631.9 -3615.05 -16479;
20274.1 -2695.25 -17084.9;
19878 -1770.3 -17658.4;
19444.3 -841.972 -18198.4;
18973.6 87.9866 -18703.9;
18467.1 1017.81 -19173.9;
17925.5 1945.73 -19607.5;
17349.9 2869.99 -20003.9;
16741.3 3788.83 -20362.4;
16101 4700.5 -20682.1;
15430.1 5603.27 -20962.5;
14729.9 6495.41 -21203.1;
14001.7 7375.23 -21403.4;
13246.9 8241.04 -21563;
12466.9 9091.2 -21681.5;
11663.1 9924.07 -21758.7;
10837.1 10738.1 -21794.5;
9990.47 11531.6 -21788.8;
9124.78 12303.2 -21741.5;
8241.68 13051.3 -21652.8;
7342.85 13774.5 -21522.7;
6429.99 14471.5 -21351.5;
5504.84 15140.8 -21139.6;
4569.15 15781.2 -20887.3;
3624.71 16391.4 -20595;
2673.32 16970.3 -20263.3;
1716.8 17516.7 -19892.9;
756.967 18029.5 -19484.3;
-204.333 18507.9 -19038.5;
-1165.26 18950.7 -18556.1;
-2123.98 19357.2 -18038.2;
-3078.64 19726.5 -17485.6;
-4027.42 20058 -16899.5;
-4968.48 20350.9 -16280.9;
-5900.02 20604.7 -15631;
-6820.25 20818.9 -14951.1;
-7727.37 20993 -14242.4;
-8619.65 21126.6 -13506.3;
-9495.36 21219.6 -12744.2;
-10352.8 21271.7 -11957.6;
-11190.3 21282.8 -11147.9;
-12006.3 21252.8 -10316.7;
-12799.1 21181.9 -9465.65;
-13567.3 21070 -8596.33;
-14309.3 20917.4 -7710.42;
-15023.6 20724.4 -6809.63;
-15709 20491.3 -5895.7;
-16364.1 20218.6 -4970.38;
-16987.5 19906.8 -4035.46;
-17578.1 19556.5 -3092.75;
-18134.7 19168.3 -2144.06;
-18656.3 18743 -1191.24;
-19141.7 18281.4 -236.13;
-19590.2 17784.4 719.429;
-20000.7 17252.9 1673.58;
-20372.5 16688 2624.48;
-20704.8 16090.7 3570.28;
-20997.1 15462.3 4509.14;
-21248.6 14803.8 5439.26;
-21459.1 14116.6 6358.82;
-21627.9 13402.1 7266.04;
-21754.9 12661.6 8159.16;
-21839.7 11896.5 9036.45;
-21882.2 11108.3 9896.21;
-21882.3 10298.5 10736.8;
-21839.9 9468.81 11556.5;
-21755.3 8620.72 12353.8;
-21628.6 7755.9 13127.2;
-21459.9 6876.04 13875.1;
-21249.7 5982.84 14596.1;
-20998.3 5078.03 15288.8;
-20706.2 4163.36 15951.9;
-20374.1 3240.62 16584;
-20002.5 2311.59 17184.1;
-19592.3 1378.07 17750.8;
-19144.1 441.863 18283.2;
-18658.9 -495.209 18780.2;
-18137.6 -1431.34 19240.8;
-17581.2 -2364.71 19664.3;
-16990.9 -3293.53 20049.6;
-16367.8 -4215.99 20396.3;
-15713 -5130.33 20703.5;
-15028 -6034.77 20970.7;
-14313.9 -6927.57 21197.4;
-13572.2 -7807.02 21383.2;
-12804.4 -8671.42 21527.7;
-12011.9 -9519.12 21630.7;
-11196.3 -10348.5 21692;
-10359.1 -11157.9 21711.5;
-9501.96 -11945.9 21689.2;
-8626.57 -12710.9 21625.1;
-7734.61 -13451.5 21519.4;
-6827.79 -14166.2 21372.3;
-5907.87 -14853.7 21184.1;
-4976.62 -15512.6 20955.2;
-4035.83 -16141.9 20686;
-3087.32 -16740.1 20377.1;
-2132.91 -17306.3 20029.2;
-1174.43 -17839.3 19642.8;
-213.73 -18338.1 19218.7;
747.361 -18801.8 18757.8;
1707 -19229.6 18261;
2663.34 -19620.6 17729.2;
3614.57 -19974.1 17163.5;
4558.86 -20289.4 16565;
5494.41 -20566 15934.7;
6419.44 -20803.3 15274.1;
7332.19 -21001 14584.3;
8230.93 -21158.6 13866.6;
9113.93 -21275.9 13122.4;
9979.54 -21352.7 12353.2;
10826.1 -21388.8 11560.5;
11652 -21384.2 10745.7;
12455.7 -21338.9 9910.47;
13235.7 -21253.1 9056.33;
13990.4 -21126.8 8184.94;
14718.5 -20960.5 7297.97;
15418.6 -20754.3 6397.1;
16089.4 -20508.7 5484.06;
16729.6 -20224.2 4560.58;
17338 -19901.3 3628.41];
