function p = positionsGJ
p=zeros(1, 1009, 3 );
p(1,:,:) = [ 26707 -164.544 -164.544;
26707.4 -123.409 -123.409;
26707.8 -82.2729 -82.2729;
26707.9 -41.1365 -41.1365;
26708 -2.76265e-23 1.89024e-23;
26707.8 41.1363 41.1363;
26707.5 82.2723 82.2723;
26707.1 123.408 123.408;
26706.5 164.542 164.542;
26705.9 205.677 205.677;
26704.9 246.808 246.808;
26703.3 287.934 287.934;
26701.2 329.056 329.056;
26698.6 370.172 370.171;
26695.6 411.28 411.28;
26692 452.382 452.382;
26687.9 493.475 493.475;
26683.3 534.559 534.559;
26678.2 575.634 575.634;
26672.6 616.698 616.698;
26666.5 657.75 657.75;
26659.9 698.79 698.79;
26652.8 739.818 739.818;
26645.2 780.832 780.832;
26637.1 821.831 821.831;
26628.4 862.815 862.815;
26619.3 903.783 903.783;
26609.7 944.734 944.734;
26599.6 985.668 985.667;
26588.9 1026.58 1026.58;
26577.8 1067.48 1067.48;
26566.1 1108.35 1108.35;
26554 1149.21 1149.21;
26541.3 1190.04 1190.04;
26528.2 1230.85 1230.85;
26514.5 1271.64 1271.64;
26500.3 1312.41 1312.41;
26485.6 1353.15 1353.14;
26470.4 1393.86 1393.86;
26454.7 1434.55 1434.54;
26438.5 1475.21 1475.2;
26421.8 1515.84 1515.84;
26404.6 1556.44 1556.44;
26386.8 1597.01 1597.01;
26368.6 1637.56 1637.55;
26349.8 1678.07 1678.07;
26330.6 1718.55 1718.54;
26310.8 1758.99 1758.99;
26290.5 1799.4 1799.4;
26269.7 1839.78 1839.78;
26248.4 1880.12 1880.12;
26226.6 1920.43 1920.43;
26204.2 1960.7 1960.69;
26181.4 2000.93 2000.92;
26158 2041.12 2041.12;
26134.2 2081.27 2081.27;
26109.8 2121.38 2121.38;
26084.8 2161.45 2161.45;
26059.4 2201.47 2201.47;
26033.5 2241.46 2241.45;
26007 2281.4 2281.39;
25980 2321.29 2321.29;
25952.5 2361.14 2361.13;
25924.5 2400.94 2400.93;
25896 2440.69 2440.69;
25866.9 2480.4 2480.39;
25837.3 2520.05 2520.05;
25807.2 2559.65 2559.65;
25776.6 2599.21 2599.2;
25745.5 2638.71 2638.7;
25713.8 2678.15 2678.15;
25681.6 2717.54 2717.54;
25648.9 2756.88 2756.87;
25615.6 2796.16 2796.15;
25581.9 2835.38 2835.38;
25547.5 2874.55 2874.54;
25512.7 2913.65 2913.65;
25477.3 2952.7 2952.69;
25441.5 2991.68 2991.68;
25405 3030.61 3030.6;
25368.1 3069.47 3069.46;
25330.6 3108.26 3108.25;
25292.5 3146.99 3146.98;
25254 3185.65 3185.64;
25214.9 3224.25 3224.24;
25175.2 3262.78 3262.77;
25135.1 3301.23 3301.22;
25094.3 3339.62 3339.61;
25053.1 3377.94 3377.93;
25011.3 3416.18 3416.17;
24969 3454.35 3454.34;
24926.1 3492.45 3492.43;
24882.6 3530.46 3530.45;
24838.7 3568.41 3568.39;
24794.2 3606.27 3606.26;
24749.1 3644.05 3644.04;
24703.5 3681.76 3681.74;
24657.3 3719.38 3719.37;
24610.6 3756.92 3756.9;
24563.3 3794.37 3794.36;
24515.5 3831.75 3831.73;
24467.1 3869.03 3869.01;
24418.2 3906.22 3906.21;
24368.7 3943.33 3943.31;
24318.6 3980.35 3980.33;
24268 4017.28 4017.26;
24216.8 4054.11 4054.09;
24165.1 4090.85 4090.83;
24112.8 4127.49 4127.47;
24059.9 4164.04 4164.02;
24006.5 4200.49 4200.47;
23952.5 4236.84 4236.82;
23897.9 4273.09 4273.07;
23842.7 4309.24 4309.21;
23787 4345.28 4345.26;
23730.7 4381.22 4381.2;
23673.8 4417.06 4417.03;
23616.4 4452.79 4452.76;
23558.4 4488.4 4488.38;
23499.8 4523.91 4523.88;
23440.6 4559.31 4559.28;
23380.8 4594.59 4594.56;
23320.4 4629.76 4629.73;
23259.5 4664.81 4664.78;
23197.9 4699.75 4699.71;
23135.8 4734.56 4734.52;
23073 4769.25 4769.22;
23009.7 4803.82 4803.79;
22945.8 4838.27 4838.23;
22881.3 4872.59 4872.55;
22816.2 4906.79 4906.75;
22750.4 4940.85 4940.81;
22684.1 4974.78 4974.74;
22617.2 5008.58 5008.54;
22549.6 5042.25 5042.21;
22481.5 5075.78 5075.74;
22412.7 5109.17 5109.13;
22343.3 5142.43 5142.38;
22273.3 5175.54 5175.49;
22202.7 5208.51 5208.46;
22131.5 5241.33 5241.28;
22059.6 5274 5273.95;
21987.1 5306.53 5306.48;
21914 5338.91 5338.85;
21840.2 5371.13 5371.07;
21765.8 5403.2 5403.14;
21690.8 5435.11 5435.05;
21615.2 5466.86 5466.8;
21538.9 5498.45 5498.38;
21461.9 5529.87 5529.81;
21384.4 5561.13 5561.07;
21306.1 5592.23 5592.16;
21227.2 5623.15 5623.08;
21147.7 5653.9 5653.83;
21067.5 5684.47 5684.4;
20986.7 5714.87 5714.8;
20905.2 5745.09 5745.02;
20823 5775.13 5775.05;
20740.2 5804.98 5804.9;
20656.6 5834.65 5834.57;
20572.5 5864.13 5864.04;
20487.6 5893.41 5893.33;
20402.1 5922.5 5922.41;
20315.9 5951.4 5951.31;
20229 5980.09 5980;
20141.4 6008.58 6008.49;
20053.1 6036.87 6036.77;
19964.1 6064.95 6064.85;
19874.5 6092.81 6092.71;
19784.1 6120.47 6120.36;
19693 6147.9 6147.8;
19601.2 6175.12 6175.01;
19508.8 6202.11 6202;
19415.6 6228.88 6228.77;
19321.6 6255.42 6255.3;
19227 6281.72 6281.61;
19131.6 6307.79 6307.67;
19035.5 6333.62 6333.5;
18938.7 6359.21 6359.08;
18841.1 6384.55 6384.42;
18742.8 6409.64 6409.51;
18643.8 6434.48 6434.35;
18544 6459.06 6458.92;
18443.5 6483.38 6483.24;
18342.2 6507.43 6507.29;
18240.1 6531.22 6531.07;
18137.3 6554.73 6554.59;
18033.7 6577.97 6577.82;
17929.3 6600.92 6600.77;
17824.2 6623.59 6623.44;
17718.3 6645.98 6645.82;
17611.5 6668.06 6667.9;
17504 6689.85 6689.69;
17395.7 6711.34 6711.17;
17286.7 6732.52 6732.34;
17176.8 6753.38 6753.21;
17066 6773.93 6773.75;
16954.5 6794.16 6793.97;
16842.2 6814.06 6813.87;
16729 6833.62 6833.43;
16615 6852.85 6852.66;
16500.2 6871.74 6871.54;
16384.5 6890.27 6890.07;
16268 6908.45 6908.25;
16150.6 6926.28 6926.06;
16032.3 6943.73 6943.51;
15913.2 6960.81 6960.59;
15793.3 6977.51 6977.29;
15672.4 6993.83 6993.6;
15550.7 7009.76 7009.52;
15428.1 7025.29 7025.05;
15304.6 7040.41 7040.16;
15180.2 7055.12 7054.87;
15054.9 7069.41 7069.15;
14928.7 7083.27 7083.01;
14801.5 7096.7 7096.43;
14673.5 7109.69 7109.42;
14544.5 7122.23 7121.95;
14414.5 7134.3 7134.02;
14283.6 7145.92 7145.63;
14151.8 7157.05 7156.76;
14019 7167.7 7167.4;
13885.2 7177.86 7177.55;
13750.5 7187.52 7187.2;
13614.7 7196.66 7196.34;
13478 7205.29 7204.96;
13340.3 7213.38 7213.04;
13201.6 7220.92 7220.58;
13061.8 7227.92 7227.57;
12921.1 7234.35 7233.99;
12779.3 7240.2 7239.84;
12636.4 7245.47 7245.1;
12492.5 7250.14 7249.76;
12347.6 7254.2 7253.81;
12201.6 7257.63 7257.23;
12054.5 7260.43 7260.02;
11906.3 7262.57 7262.16;
11757 7264.05 7263.63;
11606.6 7264.85 7264.42;
11455.1 7264.96 7264.52;
11302.5 7264.36 7263.9;
11148.8 7263.03 7262.56;
10993.9 7260.96 7260.48;
10837.8 7258.13 7257.65;
10680.6 7254.52 7254.03;
10522.2 7250.13 7249.62;
10362.7 7244.91 7244.39;
10201.9 7238.87 7238.34;
10039.9 7231.97 7231.43;
9876.69 7224.2 7223.65;
9712.26 7215.53 7214.97;
9546.59 7205.95 7205.37;
9379.66 7195.43 7194.83;
9211.47 7183.94 7183.33;
9042 7171.46 7170.84;
8871.24 7157.97 7157.33;
8699.18 7143.43 7142.78;
8525.81 7127.82 7127.15;
8351.12 7111.11 7110.43;
8175.08 7093.28 7092.57;
7997.7 7074.28 7073.56;
7818.96 7054.08 7053.34;
7638.85 7032.66 7031.9;
7457.35 7009.97 7009.19;
7274.45 6985.97 6985.18;
7090.15 6960.64 6959.82;
6904.42 6933.92 6933.08;
6717.26 6905.78 6904.91;
6528.67 6876.16 6875.28;
6338.61 6845.03 6844.12;
6147.1 6812.33 6811.4;
5954.11 6778.02 6777.06;
5759.63 6742.03 6741.04;
5563.67 6704.32 6703.3;
5366.21 6664.81 6663.77;
5167.24 6623.46 6622.39;
4966.76 6580.2 6579.09;
4764.76 6534.94 6533.8;
4561.24 6487.63 6486.46;
4356.2 6438.19 6436.98;
4149.65 6386.53 6385.28;
3941.57 6332.58 6331.29;
3731.99 6276.23 6274.9;
3520.9 6217.4 6216.02;
3308.32 6155.98 6154.56;
3094.27 6091.87 6090.4;
2878.76 6024.95 6023.43;
2661.82 5955.11 5953.54;
2443.48 5882.22 5880.59;
2223.79 5806.14 5804.45;
2002.78 5726.74 5724.99;
1780.51 5643.86 5642.04;
1557.05 5557.34 5555.44;
1332.47 5467.01 5465.04;
1106.87 5372.69 5370.64;
880.352 5274.19 5272.05;
653.039 5171.31 5169.07;
425.077 5063.82 5061.49;
196.635 4951.5 4949.07;
-32.0879 4834.12 4831.57;
-260.865 4711.4 4708.73;
-489.43 4583.1 4580.3;
-717.477 4448.91 4445.97;
-944.65 4308.56 4305.47;
-1170.54 4161.74 4158.49;
-1394.67 4008.13 4004.7;
-1616.5 3847.42 3843.8;
-1835.41 3679.29 3675.46;
-2050.68 3503.41 3499.37;
-2261.49 3319.51 3315.22;
-2466.93 3127.29 3122.74;
-2665.97 2926.51 2921.68;
-2857.43 2717 2711.87;
-3040.04 2498.64 2493.19;
-3212.42 2271.42 2265.62;
-3373.06 2035.43 2029.27;
-3520.39 1790.95 1784.41;
-3652.82 1538.39 1531.45;
-3768.74 1278.37 1271.03;
-3866.64 1011.72 1003.97;
-3945.18 739.46 731.319;
-4003.21 462.809 454.287;
-4039.93 183.132 174.256;
-4054.86 -98.1017 -107.295;
-4047.93 -379.379 -388.839;
-4019.45 -659.208 -668.877;
-3970.12 -936.182 -945.993;
-3900.93 -1209.03 -1218.91;
-3813.15 -1476.66 -1486.53;
-3708.18 -1738.16 -1747.95;
-3587.56 -1992.84 -2002.47;
-3452.83 -2240.17 -2249.58;
-3305.53 -2479.83 -2488.94;
-3147.13 -2711.62 -2720.36;
-2979.02 -2935.48 -2943.8;
-2802.46 -3151.46 -3159.3;
-2618.62 -3359.66 -3366.98;
-2428.54 -3560.27 -3567.02;
-2233.16 -3753.52 -3759.65;
-2033.3 -3939.64 -3945.13;
-1829.69 -4118.91 -4123.73;
-1622.98 -4291.6 -4295.73;
-1413.73 -4458 -4461.4;
-1202.43 -4618.38 -4621.04;
-989.502 -4773.02 -4774.92;
-775.326 -4922.17 -4923.31;
-560.227 -5066.11 -5066.46;
-344.486 -5205.06 -5204.62;
-128.348 -5339.27 -5338.02;
87.9722 -5468.95 -5466.89;
304.291 -5594.32 -5591.45;
520.446 -5715.57 -5711.88;
736.3 -5832.89 -5828.38;
951.733 -5946.47 -5941.13;
1166.64 -6056.46 -6050.29;
1380.93 -6163.03 -6156.03;
1594.52 -6266.33 -6258.49;
1807.36 -6366.49 -6357.82;
2019.38 -6463.65 -6454.15;
2230.53 -6557.93 -6547.61;
2440.77 -6649.46 -6638.3;
2650.07 -6738.33 -6726.36;
2858.39 -6824.67 -6811.87;
3065.72 -6908.57 -6894.94;
3272.02 -6990.12 -6975.67;
3477.3 -7069.4 -7054.14;
3681.52 -7146.52 -7130.44;
3884.68 -7221.54 -7204.65;
4086.77 -7294.54 -7276.84;
4287.79 -7365.59 -7347.09;
4487.73 -7434.77 -7415.46;
4686.59 -7502.13 -7482.02;
4884.37 -7567.74 -7546.83;
5081.07 -7631.65 -7609.95;
5276.69 -7693.93 -7671.44;
5471.24 -7754.62 -7731.34;
5664.71 -7813.78 -7789.71;
5857.12 -7871.44 -7846.6;
6048.47 -7927.67 -7902.06;
6238.76 -7982.5 -7956.12;
6428 -8035.98 -8008.82;
6616.2 -8088.14 -8060.22;
6803.37 -8139.03 -8110.35;
6989.51 -8188.68 -8159.24;
7174.62 -8237.12 -8206.92;
7358.73 -8284.39 -8253.45;
7541.83 -8330.52 -8298.83;
7723.94 -8375.55 -8343.11;
7905.06 -8419.49 -8386.32;
8085.2 -8462.39 -8428.48;
8264.37 -8504.27 -8469.63;
8442.58 -8545.15 -8509.78;
8619.84 -8585.06 -8548.97;
8796.15 -8624.03 -8587.21;
8971.53 -8662.07 -8624.54;
9145.98 -8699.21 -8660.96;
9319.52 -8735.48 -8696.52;
9492.14 -8770.89 -8731.22;
9663.87 -8805.46 -8765.09;
9834.7 -8839.21 -8798.14;
10004.7 -8872.17 -8830.4;
10173.7 -8904.34 -8861.88;
10341.9 -8935.76 -8892.6;
10509.3 -8966.42 -8922.58;
10675.8 -8996.36 -8951.84;
10841.4 -9025.59 -8980.38;
11006.3 -9054.11 -9008.23;
11170.2 -9081.96 -9035.4;
11333.4 -9109.13 -9061.9;
11495.8 -9135.65 -9087.76;
11657.3 -9161.53 -9112.97;
11818.1 -9186.78 -9137.56;
11978.1 -9211.42 -9161.54;
12137.3 -9235.45 -9184.91;
12295.7 -9258.89 -9207.7;
12453.3 -9281.75 -9229.91;
12610.2 -9304.04 -9251.55;
12766.4 -9325.77 -9272.64;
12921.8 -9346.96 -9293.19;
13076.4 -9367.61 -9313.2;
13230.3 -9387.73 -9332.69;
13383.5 -9407.34 -9351.66;
13536 -9426.43 -9370.13;
13687.8 -9445.03 -9388.1;
13838.8 -9463.14 -9405.59;
13989.2 -9480.77 -9422.59;
14138.8 -9497.92 -9439.13;
14287.8 -9514.61 -9455.2;
14436.1 -9530.85 -9470.82;
14583.7 -9546.63 -9486;
14730.6 -9561.97 -9500.74;
14876.9 -9576.88 -9515.04;
15022.5 -9591.37 -9528.93;
15167.5 -9605.43 -9542.39;
15311.8 -9619.09 -9555.45;
15455.4 -9632.33 -9568.1;
15598.4 -9645.18 -9580.35;
15740.8 -9657.63 -9592.22;
15882.6 -9669.7 -9603.7;
16023.7 -9681.39 -9614.8;
16164.2 -9692.7 -9625.53;
16304.2 -9703.64 -9635.9;
16443.4 -9714.22 -9645.9;
16582.1 -9724.44 -9655.55;
16720.2 -9734.31 -9664.84;
16857.7 -9743.83 -9673.8;
16994.6 -9753.01 -9682.41;
17131 -9761.85 -9690.69;
17266.7 -9770.36 -9698.64;
17401.9 -9778.55 -9706.26;
17536.5 -9786.41 -9713.56;
17670.5 -9793.95 -9720.55;
17803.9 -9801.18 -9727.23;
17936.8 -9808.1 -9733.6;
18069.2 -9814.72 -9739.67;
18201 -9821.03 -9745.43;
18332.2 -9827.05 -9750.91;
18462.9 -9832.78 -9756.09;
18593.1 -9838.22 -9760.99;
18722.7 -9843.37 -9765.61;
18851.8 -9848.24 -9769.95;
18980.3 -9852.84 -9774.01;
19108.4 -9857.16 -9777.8;
19235.9 -9861.22 -9781.33;
19362.9 -9865 -9784.59;
19489.4 -9868.53 -9787.59;
19615.3 -9871.79 -9790.33;
19740.8 -9874.8 -9792.82;
19865.8 -9877.56 -9795.06;
19990.2 -9880.06 -9797.05;
20114.2 -9882.32 -9798.79;
20237.7 -9884.34 -9800.3;
20360.6 -9886.12 -9801.56;
20483.1 -9887.66 -9802.59;
20605.1 -9888.96 -9803.39;
20726.6 -9890.03 -9803.96;
20847.7 -9890.88 -9804.3;
20968.3 -9891.49 -9804.42;
21088.3 -9891.89 -9804.31;
21208 -9892.06 -9803.99;
21327.1 -9892.01 -9803.45;
21445.8 -9891.75 -9802.7;
21564.1 -9891.28 -9801.73;
21681.8 -9890.59 -9800.56;
21799.1 -9889.7 -9799.18;
21916 -9888.6 -9797.59;
22032.4 -9887.3 -9795.81;
22148.4 -9885.8 -9793.82;
22263.9 -9884.09 -9791.64;
22379 -9882.19 -9789.26;
22493.6 -9880.1 -9786.69;
22607.8 -9877.81 -9783.93;
22721.6 -9875.34 -9780.98;
22834.9 -9872.67 -9777.84;
22947.8 -9869.82 -9774.52;
23060.3 -9866.78 -9771.02;
23172.4 -9863.57 -9767.34;
23284 -9860.17 -9763.47;
23395.2 -9856.59 -9759.44;
23506 -9852.84 -9755.22;
23616.4 -9848.91 -9750.84;
23726.4 -9844.81 -9746.28;
23835.9 -9840.54 -9741.55;
23945 -9836.1 -9736.66;
24053.8 -9831.49 -9731.6;
24162.1 -9826.72 -9726.37;
24270 -9821.78 -9720.99;
24377.6 -9816.68 -9715.44;
24484.7 -9811.42 -9709.73;
24591.4 -9806 -9703.87;
24697.8 -9800.42 -9697.85;
24803.7 -9794.69 -9691.68;
24909.3 -9788.8 -9685.35;
25014.4 -9782.76 -9678.87;
25119.2 -9776.57 -9672.24;
25223.6 -9770.23 -9665.47;
25327.6 -9763.74 -9658.54;
25431.3 -9757.1 -9651.48;
25534.5 -9750.31 -9644.26;
25637.4 -9743.39 -9636.91;
25739.9 -9736.32 -9629.41;
25842 -9729.11 -9621.78;
25943.8 -9721.76 -9614;
26045.2 -9714.27 -9606.09;
26146.2 -9706.64 -9598.04;
26246.9 -9698.87 -9589.86;
26347.2 -9690.97 -9581.54;
26447.1 -9682.94 -9573.09;
26546.7 -9674.78 -9564.51;
26645.9 -9666.48 -9555.81;
26744.7 -9658.06 -9546.97;
26843.2 -9649.5 -9538;
26941.3 -9640.82 -9528.91;
27039.1 -9632.01 -9519.7;
27136.6 -9623.07 -9510.36;
27233.7 -9614.02 -9500.89;
27330.4 -9604.83 -9491.31;
27426.8 -9595.53 -9481.6;
27522.8 -9586.1 -9471.78;
27618.5 -9576.56 -9461.83;
27713.9 -9566.9 -9451.77;
27808.9 -9557.12 -9441.6;
27903.6 -9547.22 -9431.3;
27997.9 -9537.2 -9420.9;
28091.9 -9527.07 -9410.38;
28185.6 -9516.83 -9399.75;
28278.9 -9506.48 -9389;
28371.9 -9496.01 -9378.15;
28464.6 -9485.43 -9367.18;
28557 -9474.74 -9356.11;
28649 -9463.95 -9344.93;
28740.7 -9453.04 -9333.64;
28832 -9442.03 -9322.25;
28923.1 -9430.91 -9310.75;
29013.8 -9419.69 -9299.15;
29104.2 -9408.36 -9287.45;
29194.3 -9396.93 -9275.64;
29284 -9385.39 -9263.73;
29373.5 -9373.75 -9251.72;
29462.6 -9362.02 -9239.61;
29551.4 -9350.18 -9227.4;
29639.9 -9338.24 -9215.1;
29728.1 -9326.2 -9202.69;
29815.9 -9314.07 -9190.19;
29903.5 -9301.84 -9177.6;
29990.7 -9289.51 -9164.9;
30077.7 -9277.09 -9152.12;
30164.3 -9264.57 -9139.24;
30250.7 -9251.96 -9126.27;
30336.7 -9239.25 -9113.2;
30422.4 -9226.46 -9100.05;
30507.8 -9213.57 -9086.8;
30592.9 -9200.58 -9073.47;
30677.8 -9187.51 -9060.04;
30762.3 -9174.35 -9046.53;
30846.5 -9161.1 -9032.93;
30930.4 -9147.77 -9019.24;
31014.1 -9134.34 -9005.47;
31097.4 -9120.83 -8991.61;
31180.4 -9107.23 -8977.66;
31263.2 -9093.54 -8963.63;
31345.6 -9079.77 -8949.52;
31427.8 -9065.92 -8935.32;
31509.7 -9051.98 -8921.04;
31591.3 -9037.96 -8906.68;
31672.6 -9023.86 -8892.24;
31753.6 -9009.68 -8877.72;
31834.3 -8995.41 -8863.12;
31914.8 -8981.06 -8848.43;
31994.9 -8966.64 -8833.67;
32074.8 -8952.13 -8818.83;
32154.4 -8937.55 -8803.92;
32233.7 -8922.88 -8788.92;
32312.7 -8908.14 -8773.85;
32391.5 -8893.33 -8758.71;
32470 -8878.43 -8743.49;
32548.2 -8863.46 -8728.19;
32626.1 -8848.42 -8712.82;
32703.7 -8833.3 -8697.38;
32781.1 -8818.1 -8681.86;
32858.2 -8802.84 -8666.27;
32935 -8787.49 -8650.61;
33011.6 -8772.08 -8634.87;
33087.9 -8756.59 -8619.07;
33163.9 -8741.04 -8603.2;
33239.6 -8725.41 -8587.25;
33315.1 -8709.71 -8571.24;
33390.3 -8693.94 -8555.15;
33465.3 -8678.1 -8539;
33539.9 -8662.19 -8522.78;
33614.4 -8646.22 -8506.49;
33688.5 -8630.17 -8490.14;
33762.4 -8614.06 -8473.72;
33836 -8597.88 -8457.23;
33909.4 -8581.63 -8440.68;
33982.4 -8565.32 -8424.06;
34055.3 -8548.94 -8407.38;
34127.9 -8532.5 -8390.63;
34200.2 -8515.99 -8373.82;
34272.2 -8499.41 -8356.95;
34344 -8482.78 -8340.01;
34415.6 -8466.08 -8323.01;
34486.9 -8449.31 -8305.95;
34557.9 -8432.48 -8288.82;
34628.7 -8415.6 -8271.64;
34699.2 -8398.64 -8254.39;
34769.5 -8381.63 -8237.09;
34839.5 -8364.56 -8219.72;
34909.2 -8347.43 -8202.3;
34978.7 -8330.23 -8184.81;
35048 -8312.98 -8167.27;
35117 -8295.67 -8149.67;
35185.8 -8278.29 -8132.01;
35254.3 -8260.86 -8114.29;
35322.6 -8243.38 -8096.52;
35390.6 -8225.83 -8078.69;
35458.4 -8208.23 -8060.8;
35525.9 -8190.56 -8042.86;
35593.2 -8172.85 -8024.86;
35660.2 -8155.07 -8006.81;
35727 -8137.24 -7988.7;
35793.6 -8119.36 -7970.54;
35859.9 -8101.42 -7952.32;
35925.9 -8083.43 -7934.05;
35991.8 -8065.38 -7915.73;
36057.3 -8047.27 -7897.35;
36122.7 -8029.12 -7878.92;
36187.8 -8010.91 -7860.44;
36252.7 -7992.64 -7841.9;
36317.3 -7974.33 -7823.32;
36381.7 -7955.96 -7804.68;
36445.8 -7937.54 -7785.99;
36509.7 -7919.07 -7767.25;
36573.4 -7900.55 -7748.47;
36636.9 -7881.97 -7729.63;
36700.1 -7863.35 -7710.74;
36763.1 -7844.67 -7691.8;
36825.8 -7825.95 -7672.81;
36888.3 -7807.17 -7653.78;
36950.6 -7788.35 -7634.69;
37012.6 -7769.48 -7615.56;
37074.4 -7750.56 -7596.38;
37136 -7731.59 -7577.16;
37197.4 -7712.57 -7557.88;
37258.5 -7693.51 -7538.56;
37319.4 -7674.39 -7519.19;
37380 -7655.23 -7499.78;
37440.5 -7636.03 -7480.32;
37500.7 -7616.77 -7460.82;
37560.7 -7597.47 -7441.27;
37620.4 -7578.13 -7421.67;
37680 -7558.74 -7402.03;
37739.3 -7539.3 -7382.35;
37798.3 -7519.82 -7362.62;
37857.2 -7500.3 -7342.85;
37915.8 -7480.72 -7323.03;
37974.2 -7461.11 -7303.17;
38032.4 -7441.45 -7283.27;
38090.4 -7421.75 -7263.33;
38148.1 -7402 -7243.34;
38205.6 -7382.22 -7223.31;
38262.9 -7362.38 -7203.24;
38320 -7342.51 -7183.13;
38376.8 -7322.59 -7162.97;
38433.5 -7302.63 -7142.78;
38489.9 -7282.63 -7122.54;
38546.1 -7262.59 -7102.26;
38602 -7242.51 -7081.95;
38657.8 -7222.39 -7061.59;
38713.3 -7202.22 -7041.19;
38768.7 -7182.02 -7020.75;
38823.8 -7161.77 -7000.28;
38878.6 -7141.49 -6979.76;
38933.3 -7121.16 -6959.21;
38987.8 -7100.8 -6938.62;
39042 -7080.39 -6917.99;
39096.1 -7059.95 -6897.32;
39149.9 -7039.47 -6876.61;
39203.5 -7018.95 -6855.87;
39256.9 -6998.39 -6835.09;
39310 -6977.8 -6814.27;
39363 -6957.16 -6793.41;
39415.7 -6936.49 -6772.52;
39468.3 -6915.78 -6751.59;
39520.6 -6895.04 -6730.63;
39572.7 -6874.25 -6709.63;
39624.6 -6853.43 -6688.59;
39676.3 -6832.58 -6667.52;
39727.8 -6811.69 -6646.41;
39779.1 -6790.76 -6625.27;
39830.2 -6769.8 -6604.09;
39881 -6748.8 -6582.88;
39931.7 -6727.76 -6561.63;
39982.1 -6706.7 -6540.35;
40032.4 -6685.59 -6519.04;
40082.4 -6664.45 -6497.69;
40132.2 -6643.28 -6476.31;
40181.8 -6622.07 -6454.9;
40231.2 -6600.83 -6433.45;
40280.4 -6579.56 -6411.97;
40329.5 -6558.25 -6390.45;
40378.2 -6536.91 -6368.91;
40426.8 -6515.53 -6347.33;
40475.2 -6494.12 -6325.72;
40523.4 -6472.68 -6304.08;
40571.4 -6451.21 -6282.4;
40619.2 -6429.71 -6260.7;
40666.8 -6408.17 -6238.96;
40714.1 -6386.6 -6217.19;
40761.3 -6365 -6195.39;
40808.3 -6343.37 -6173.56;
40855 -6321.7 -6151.7;
40901.6 -6300.01 -6129.81;
40948 -6278.28 -6107.89;
40994.1 -6256.52 -6085.94;
41040.1 -6234.74 -6063.96;
41085.9 -6212.92 -6041.95;
41131.4 -6191.07 -6019.91;
41176.8 -6169.19 -5997.85;
41222 -6147.29 -5975.75;
41267 -6125.35 -5953.62;
41311.7 -6103.38 -5931.47;
41356.3 -6081.39 -5909.29;
41400.7 -6059.36 -5887.08;
41444.9 -6037.31 -5864.84;
41488.9 -6015.22 -5842.57;
41532.7 -5993.11 -5820.27;
41576.3 -5970.97 -5797.95;
41619.7 -5948.81 -5775.6;
41662.9 -5926.61 -5753.23;
41705.9 -5904.39 -5730.82;
41748.7 -5882.14 -5708.39;
41791.3 -5859.86 -5685.93;
41833.7 -5837.55 -5663.45;
41876 -5815.22 -5640.94;
41918 -5792.86 -5618.4;
41959.8 -5770.47 -5595.84;
42001.5 -5748.06 -5573.25;
42042.9 -5725.62 -5550.64;
42084.2 -5703.15 -5528;
42125.3 -5680.66 -5505.34;
42166.2 -5658.14 -5482.65;
42206.8 -5635.6 -5459.93;
42247.3 -5613.03 -5437.19;
42287.7 -5590.44 -5414.43;
42327.8 -5567.82 -5391.64;
42367.7 -5545.17 -5368.83;
42407.4 -5522.5 -5345.99;
42447 -5499.81 -5323.13;
42486.3 -5477.09 -5300.24;
42525.5 -5454.34 -5277.34;
42564.5 -5431.57 -5254.4;
42603.2 -5408.78 -5231.45;
42641.8 -5385.96 -5208.47;
42680.3 -5363.12 -5185.47;
42718.5 -5340.26 -5162.44;
42756.5 -5317.37 -5139.4;
42794.3 -5294.46 -5116.33;
42832 -5271.52 -5093.23;
42869.5 -5248.57 -5070.12;
42906.8 -5225.59 -5046.98;
42943.9 -5202.58 -5023.82;
42980.8 -5179.56 -5000.64;
43017.5 -5156.51 -4977.44;
43054 -5133.44 -4954.21;
43090.4 -5110.34 -4930.97;
43126.5 -5087.23 -4907.7;
43162.5 -5064.09 -4884.41;
43198.3 -5040.93 -4861.1;
43233.9 -5017.75 -4837.77;
43269.3 -4994.55 -4814.42;
43304.6 -4971.33 -4791.05;
43339.6 -4948.08 -4767.66;
43374.5 -4924.82 -4744.25;
43409.2 -4901.53 -4720.81;
43443.7 -4878.22 -4697.36;
43478 -4854.9 -4673.89;
43512.2 -4831.55 -4650.4;
43546.1 -4808.18 -4626.89;
43579.9 -4784.79 -4603.35;
43613.5 -4761.38 -4579.8;
43646.9 -4737.95 -4556.23;
43680.1 -4714.5 -4532.65;
43713.1 -4691.03 -4509.04;
43746 -4667.55 -4485.41;
43778.7 -4644.04 -4461.77;
43811.2 -4620.51 -4438.1;
43843.5 -4596.96 -4414.42;
43875.6 -4573.4 -4390.72;
43907.6 -4549.82 -4367;
43939.4 -4526.21 -4343.27;
43970.9 -4502.59 -4319.51;
44002.4 -4478.95 -4295.74;
44033.6 -4455.29 -4271.95;
44064.6 -4431.62 -4248.14;
44095.5 -4407.92 -4224.32;
44126.2 -4384.21 -4200.48;
44156.7 -4360.48 -4176.62;
44187.1 -4336.73 -4152.74;
44217.2 -4312.97 -4128.85;
44247.2 -4289.18 -4104.94;
44277 -4265.38 -4081.02;
44306.6 -4241.57 -4057.07;
44336.1 -4217.73 -4033.11;
44365.3 -4193.88 -4009.14;
44394.4 -4170.01 -3985.15;
44423.3 -4146.13 -3961.14;
44452.1 -4122.22 -3937.12;
44480.6 -4098.31 -3913.08;
44509 -4074.37 -3889.03;
44537.2 -4050.42 -3864.96;
44565.2 -4026.45 -3840.87;
44593.1 -4002.47 -3816.77;
44620.8 -3978.47 -3792.66;
44648.3 -3954.46 -3768.53;
44675.6 -3930.43 -3744.38;
44702.7 -3906.38 -3720.22;
44729.7 -3882.32 -3696.04;
44756.5 -3858.24 -3671.86;
44783.1 -3834.15 -3647.65;
44809.5 -3810.05 -3623.43;
44835.8 -3785.92 -3599.2;
44861.9 -3761.79 -3574.95;
44887.8 -3737.64 -3550.69;
44913.6 -3713.47 -3526.42;
44939.1 -3689.29 -3502.13;
44964.5 -3665.1 -3477.83;
44989.7 -3640.89 -3453.52;
45014.8 -3616.67 -3429.19;
45039.7 -3592.43 -3404.85;
45064.4 -3568.18 -3380.49;
45088.9 -3543.91 -3356.12;
45113.2 -3519.64 -3331.74;
45137.4 -3495.34 -3307.35;
45161.4 -3471.04 -3282.94;
45185.2 -3446.72 -3258.52;
45208.9 -3422.39 -3234.09;
45232.4 -3398.04 -3209.65;
45255.7 -3373.69 -3185.19;
45278.8 -3349.31 -3160.72;
45301.8 -3324.93 -3136.24;
45324.6 -3300.54 -3111.75;
45347.2 -3276.13 -3087.24;
45369.7 -3251.71 -3062.73;
45392 -3227.27 -3038.2;
45414.1 -3202.83 -3013.66;
45436 -3178.37 -2989.11;
45457.8 -3153.9 -2964.55;
45479.4 -3129.42 -2939.97;
45500.8 -3104.92 -2915.39;
45522.1 -3080.42 -2890.8;
45543.1 -3055.9 -2866.19;
45564 -3031.37 -2841.57;
45584.8 -3006.83 -2816.94;
45605.4 -2982.28 -2792.31;
45625.8 -2957.72 -2767.66;
45646 -2933.14 -2743;
45666 -2908.56 -2718.33;
45685.9 -2883.96 -2693.65;
45705.6 -2859.36 -2668.96;
45725.2 -2834.74 -2644.26;
45744.6 -2810.11 -2619.55;
45763.8 -2785.48 -2594.83;
45782.8 -2760.83 -2570.1;
45801.7 -2736.17 -2545.36;
45820.4 -2711.5 -2520.62;
45838.9 -2686.82 -2495.86;
45857.3 -2662.13 -2471.09;
45875.5 -2637.44 -2446.32;
45893.5 -2612.73 -2421.53;
45911.4 -2588.01 -2396.74;
45929 -2563.28 -2371.94;
45946.6 -2538.55 -2347.12;
45963.9 -2513.8 -2322.3;
45981.1 -2489.04 -2297.48;
45998.1 -2464.28 -2272.64;
46014.9 -2439.51 -2247.79;
46031.6 -2414.72 -2222.94;
46048.1 -2389.93 -2198.08;
46064.5 -2365.13 -2173.21;
46080.6 -2340.32 -2148.33;
46096.6 -2315.51 -2123.45;
46112.5 -2290.68 -2098.55;
46128.1 -2265.85 -2073.65;
46143.6 -2241 -2048.74;
46159 -2216.15 -2023.83;
46174.1 -2191.29 -1998.91;
46189.1 -2166.43 -1973.97;
46204 -2141.55 -1949.04;
46218.6 -2116.67 -1924.09;
46233.1 -2091.78 -1899.14;
46247.5 -2066.88 -1874.18;
46261.6 -2041.98 -1849.22;
46275.6 -2017.07 -1824.24;
46289.5 -1992.15 -1799.27;
46303.1 -1967.22 -1774.28;
46316.6 -1942.28 -1749.29;
46330 -1917.34 -1724.29;
46343.1 -1892.4 -1699.29;
46356.1 -1867.44 -1674.28;
46369 -1842.48 -1649.26;
46381.6 -1817.51 -1624.24;
46394.1 -1792.54 -1599.21;
46406.5 -1767.55 -1574.17;
46418.6 -1742.57 -1549.14;
46430.6 -1717.57 -1524.09;
46442.5 -1692.57 -1499.04;
46454.1 -1667.57 -1473.98;
46465.6 -1642.56 -1448.92;
46477 -1617.54 -1423.86;
46488.1 -1592.51 -1398.78;
46499.2 -1567.48 -1373.71;
46510 -1542.45 -1348.63;
46520.7 -1517.41 -1323.54;
46531.2 -1492.36 -1298.45;
46541.5 -1467.31 -1273.35;
46551.7 -1442.25 -1248.25;
46561.7 -1417.19 -1223.15;
46571.6 -1392.13 -1198.04;
46581.3 -1367.05 -1172.92;
46590.8 -1341.98 -1147.81;
46600.1 -1316.9 -1122.68;
46609.3 -1291.81 -1097.56;
46618.4 -1266.72 -1072.43;
46627.2 -1241.62 -1047.3;
46635.9 -1216.52 -1022.16;
46644.5 -1191.42 -997.017;
46652.8 -1166.31 -971.873;
46661 -1141.2 -946.725;
46669.1 -1116.08 -921.574;
46676.9 -1090.96 -896.419;
46684.6 -1065.84 -871.261;
46692.2 -1040.71 -846.1;
46699.6 -1015.58 -820.937;
46706.8 -990.444 -795.77;
46713.8 -965.305 -770.6;
46720.7 -940.163 -745.428;
46727.5 -915.017 -720.253;
46734 -889.869 -695.075;
46740.4 -864.716 -669.895;
46746.6 -839.561 -644.713;
46752.7 -814.403 -619.528;
46758.6 -789.242 -594.341;
46764.4 -764.079 -569.152;
46769.9 -738.912 -543.961;
46775.3 -713.743 -518.768;
46780.6 -688.571 -493.573;
46785.7 -663.397 -468.376;
46790.6 -638.221 -443.178;
46795.3 -613.042 -417.978;
46799.9 -587.861 -392.777;
46804.4 -562.678 -367.574;
46808.6 -537.493 -342.37;
46812.7 -512.306 -317.164;
46816.7 -487.118 -291.958;
46820.4 -461.927 -266.75;
46824.1 -436.735 -241.542;
46827.5 -411.542 -216.332;
46830.8 -386.346 -191.122;
46833.9 -361.15 -165.911;
46836.9 -335.952 -140.7;
46839.7 -310.753 -115.488;
46842.3 -285.553 -90.2753;
46844.7 -260.352 -65.0625;
46847 -235.15 -39.8494;
46849.2 -209.947 -14.6362;
46851.2 -184.743 10.5772;
46853 -159.539 35.7905;
46854.6 -134.334 61.0037;
46856.1 -109.128 86.2168;
46857.4 -83.9223 111.43;
46858.6 -58.7161 136.642;
46859.6 -33.5097 161.854;
46860.4 -8.30305 187.065;
46861.1 16.9036 212.276;
46861.6 42.1103 237.486;
46861.9 67.3168 262.696;
46862.1 92.5231 287.904;
46862.1 117.729 313.112;
46862 142.935 338.318;
46861.7 168.14 363.523;
46861.2 193.345 388.727;
46860.5 218.549 413.93;
46859.7 243.752 439.131;
46858.8 268.954 464.331;
46857.6 294.156 489.529;
46856.3 319.357 514.726;
46854.9 344.556 539.921;
46853.3 369.754 565.114];
