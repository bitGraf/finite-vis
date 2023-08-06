package main

import (
	"fmt"
	"strings"
)

var (
	parula_table = [...][3]float64{
		{0.2081, 0.1663, 0.5292},
		{0.2091, 0.1721, 0.5411},
		{0.2101, 0.1779, 0.553},
		{0.2109, 0.1837, 0.565},
		{0.2116, 0.1895, 0.5771},
		{0.2121, 0.1954, 0.5892},
		{0.2124, 0.2013, 0.6013},
		{0.2125, 0.2072, 0.6135},
		{0.2123, 0.2132, 0.6258},
		{0.2118, 0.2192, 0.6381},
		{0.2111, 0.2253, 0.6505},
		{0.2099, 0.2315, 0.6629},
		{0.2084, 0.2377, 0.6753},
		{0.2063, 0.244, 0.6878},
		{0.2038, 0.2503, 0.7003},
		{0.2006, 0.2568, 0.7129},
		{0.1968, 0.2632, 0.7255},
		{0.1921, 0.2698, 0.7381},
		{0.1867, 0.2764, 0.7507},
		{0.1802, 0.2832, 0.7634},
		{0.1728, 0.2902, 0.7762},
		{0.1641, 0.2975, 0.789},
		{0.1541, 0.3052, 0.8017},
		{0.1427, 0.3132, 0.8145},
		{0.1295, 0.3217, 0.8269},
		{0.1147, 0.3306, 0.8387},
		{0.0986, 0.3397, 0.8495},
		{0.0816, 0.3486, 0.8588},
		{0.0646, 0.3572, 0.8664},
		{0.0482, 0.3651, 0.8722},
		{0.0329, 0.3724, 0.8765},
		{0.0213, 0.3792, 0.8796},
		{0.0136, 0.3853, 0.8815},
		{0.0086, 0.3911, 0.8827},
		{0.006, 0.3965, 0.8833},
		{0.0051, 0.4017, 0.8834},
		{0.0054, 0.4066, 0.8831},
		{0.0067, 0.4113, 0.8825},
		{0.0089, 0.4159, 0.8816},
		{0.0116, 0.4203, 0.8805},
		{0.0148, 0.4246, 0.8793},
		{0.0184, 0.4288, 0.8779},
		{0.0223, 0.4329, 0.8763},
		{0.0264, 0.437, 0.8747},
		{0.0306, 0.441, 0.8729},
		{0.0349, 0.4449, 0.8711},
		{0.0394, 0.4488, 0.8692},
		{0.0437, 0.4526, 0.8672},
		{0.0477, 0.4564, 0.8652},
		{0.0514, 0.4602, 0.8632},
		{0.0549, 0.464, 0.8611},
		{0.0582, 0.4677, 0.8589},
		{0.0612, 0.4714, 0.8568},
		{0.064, 0.4751, 0.8546},
		{0.0666, 0.4788, 0.8525},
		{0.0689, 0.4825, 0.8503},
		{0.071, 0.4862, 0.8481},
		{0.0729, 0.4899, 0.846},
		{0.0746, 0.4937, 0.8439},
		{0.0761, 0.4974, 0.8418},
		{0.0773, 0.5012, 0.8398},
		{0.0782, 0.5051, 0.8378},
		{0.0789, 0.5089, 0.8359},
		{0.0794, 0.5129, 0.8341},
		{0.0795, 0.5169, 0.8324},
		{0.0793, 0.521, 0.8308},
		{0.0788, 0.5251, 0.8293},
		{0.0778, 0.5295, 0.828},
		{0.0764, 0.5339, 0.827},
		{0.0746, 0.5384, 0.8261},
		{0.0724, 0.5431, 0.8253},
		{0.0698, 0.5479, 0.8247},
		{0.0668, 0.5527, 0.8243},
		{0.0636, 0.5577, 0.8239},
		{0.06, 0.5627, 0.8237},
		{0.0562, 0.5677, 0.8234},
		{0.0523, 0.5727, 0.8231},
		{0.0484, 0.5777, 0.8228},
		{0.0445, 0.5826, 0.8223},
		{0.0408, 0.5874, 0.8217},
		{0.0372, 0.5922, 0.8209},
		{0.0342, 0.5968, 0.8198},
		{0.0317, 0.6012, 0.8186},
		{0.0296, 0.6055, 0.8171},
		{0.0279, 0.6097, 0.8154},
		{0.0265, 0.6137, 0.8135},
		{0.0255, 0.6176, 0.8114},
		{0.0248, 0.6214, 0.8091},
		{0.0243, 0.625, 0.8066},
		{0.0239, 0.6285, 0.8039},
		{0.0237, 0.6319, 0.801},
		{0.0235, 0.6352, 0.798},
		{0.0233, 0.6384, 0.7948},
		{0.0231, 0.6415, 0.7916},
		{0.023, 0.6445, 0.7881},
		{0.0229, 0.6474, 0.7846},
		{0.0227, 0.6503, 0.781},
		{0.0227, 0.6531, 0.7773},
		{0.0232, 0.6558, 0.7735},
		{0.0238, 0.6585, 0.7696},
		{0.0246, 0.6611, 0.7656},
		{0.0263, 0.6637, 0.7615},
		{0.0282, 0.6663, 0.7574},
		{0.0306, 0.6688, 0.7532},
		{0.0338, 0.6712, 0.749},
		{0.0373, 0.6737, 0.7446},
		{0.0418, 0.6761, 0.7402},
		{0.0467, 0.6784, 0.7358},
		{0.0516, 0.6808, 0.7313},
		{0.0574, 0.6831, 0.7267},
		{0.0629, 0.6854, 0.7221},
		{0.0692, 0.6877, 0.7173},
		{0.0755, 0.6899, 0.7126},
		{0.082, 0.6921, 0.7078},
		{0.0889, 0.6943, 0.7029},
		{0.0956, 0.6965, 0.6979},
		{0.1031, 0.6986, 0.6929},
		{0.1104, 0.7007, 0.6878},
		{0.118, 0.7028, 0.6827},
		{0.1258, 0.7049, 0.6775},
		{0.1335, 0.7069, 0.6723},
		{0.1418, 0.7089, 0.6669},
		{0.1499, 0.7109, 0.6616},
		{0.1585, 0.7129, 0.6561},
		{0.1671, 0.7148, 0.6507},
		{0.1758, 0.7168, 0.6451},
		{0.1849, 0.7186, 0.6395},
		{0.1938, 0.7205, 0.6338},
		{0.2033, 0.7223, 0.6281},
		{0.2128, 0.7241, 0.6223},
		{0.2224, 0.7259, 0.6165},
		{0.2324, 0.7275, 0.6107},
		{0.2423, 0.7292, 0.6048},
		{0.2527, 0.7308, 0.5988},
		{0.2631, 0.7324, 0.5929},
		{0.2735, 0.7339, 0.5869},
		{0.2845, 0.7354, 0.5809},
		{0.2953, 0.7368, 0.5749},
		{0.3064, 0.7381, 0.5689},
		{0.3177, 0.7394, 0.563},
		{0.3289, 0.7406, 0.557},
		{0.3405, 0.7417, 0.5512},
		{0.352, 0.7428, 0.5453},
		{0.3635, 0.7438, 0.5396},
		{0.3753, 0.7446, 0.5339},
		{0.3869, 0.7454, 0.5283},
		{0.3986, 0.7461, 0.5229},
		{0.4103, 0.7467, 0.5175},
		{0.4218, 0.7473, 0.5123},
		{0.4334, 0.7477, 0.5072},
		{0.4447, 0.7482, 0.5021},
		{0.4561, 0.7485, 0.4972},
		{0.4672, 0.7487, 0.4924},
		{0.4783, 0.7489, 0.4877},
		{0.4892, 0.7491, 0.4831},
		{0.5, 0.7491, 0.4786},
		{0.5106, 0.7492, 0.4741},
		{0.5212, 0.7492, 0.4698},
		{0.5315, 0.7491, 0.4655},
		{0.5418, 0.749, 0.4613},
		{0.5519, 0.7489, 0.4571},
		{0.5619, 0.7487, 0.4531},
		{0.5718, 0.7485, 0.449},
		{0.5816, 0.7482, 0.4451},
		{0.5913, 0.7479, 0.4412},
		{0.6009, 0.7476, 0.4374},
		{0.6103, 0.7473, 0.4335},
		{0.6197, 0.7469, 0.4298},
		{0.629, 0.7465, 0.4261},
		{0.6382, 0.746, 0.4224},
		{0.6473, 0.7456, 0.4188},
		{0.6564, 0.7451, 0.4152},
		{0.6653, 0.7446, 0.4116},
		{0.6742, 0.7441, 0.4081},
		{0.683, 0.7435, 0.4046},
		{0.6918, 0.743, 0.4011},
		{0.7004, 0.7424, 0.3976},
		{0.7091, 0.7418, 0.3942},
		{0.7176, 0.7412, 0.3908},
		{0.7261, 0.7405, 0.3874},
		{0.7346, 0.7399, 0.384},
		{0.743, 0.7392, 0.3806},
		{0.7513, 0.7385, 0.3773},
		{0.7596, 0.7378, 0.3739},
		{0.7679, 0.7372, 0.3706},
		{0.7761, 0.7364, 0.3673},
		{0.7843, 0.7357, 0.3639},
		{0.7924, 0.735, 0.3606},
		{0.8005, 0.7343, 0.3573},
		{0.8085, 0.7336, 0.3539},
		{0.8166, 0.7329, 0.3506},
		{0.8246, 0.7322, 0.3472},
		{0.8325, 0.7315, 0.3438},
		{0.8405, 0.7308, 0.3404},
		{0.8484, 0.7301, 0.337},
		{0.8563, 0.7294, 0.3336},
		{0.8642, 0.7288, 0.33},
		{0.872, 0.7282, 0.3265},
		{0.8798, 0.7276, 0.3229},
		{0.8877, 0.7271, 0.3193},
		{0.8954, 0.7266, 0.3156},
		{0.9032, 0.7262, 0.3117},
		{0.911, 0.7259, 0.3078},
		{0.9187, 0.7256, 0.3038},
		{0.9264, 0.7256, 0.2996},
		{0.9341, 0.7256, 0.2953},
		{0.9417, 0.7259, 0.2907},
		{0.9493, 0.7264, 0.2859},
		{0.9567, 0.7273, 0.2808},
		{0.9639, 0.7285, 0.2754},
		{0.9708, 0.7303, 0.2696},
		{0.9773, 0.7326, 0.2634},
		{0.9831, 0.7355, 0.257},
		{0.9882, 0.739, 0.2504},
		{0.9922, 0.7431, 0.2437},
		{0.9952, 0.7476, 0.2373},
		{0.9973, 0.7524, 0.231},
		{0.9986, 0.7573, 0.2251},
		{0.9991, 0.7624, 0.2195},
		{0.999, 0.7675, 0.2141},
		{0.9985, 0.7726, 0.209},
		{0.9976, 0.7778, 0.2042},
		{0.9964, 0.7829, 0.1995},
		{0.995, 0.788, 0.1949},
		{0.9933, 0.7931, 0.1905},
		{0.9914, 0.7981, 0.1863},
		{0.9894, 0.8032, 0.1821},
		{0.9873, 0.8083, 0.178},
		{0.9851, 0.8133, 0.174},
		{0.9828, 0.8184, 0.17},
		{0.9805, 0.8235, 0.1661},
		{0.9782, 0.8286, 0.1622},
		{0.9759, 0.8337, 0.1583},
		{0.9736, 0.8389, 0.1544},
		{0.9713, 0.8441, 0.1505},
		{0.9692, 0.8494, 0.1465},
		{0.9672, 0.8548, 0.1425},
		{0.9654, 0.8603, 0.1385},
		{0.9638, 0.8659, 0.1343},
		{0.9623, 0.8716, 0.1301},
		{0.9611, 0.8774, 0.1258},
		{0.96, 0.8834, 0.1215},
		{0.9593, 0.8895, 0.1171},
		{0.9588, 0.8958, 0.1126},
		{0.9586, 0.9022, 0.1082},
		{0.9587, 0.9088, 0.1036},
		{0.9591, 0.9155, 0.099},
		{0.9599, 0.9225, 0.0944},
		{0.961, 0.9296, 0.0897},
		{0.9624, 0.9368, 0.085},
		{0.9641, 0.9443, 0.0802},
		{0.9662, 0.9518, 0.0753},
		{0.9685, 0.9595, 0.0703},
		{0.971, 0.9673, 0.0651},
		{0.9736, 0.9752, 0.0597},
		{0.9763, 0.9831, 0.0538},
	}
	plasma_table = [...][3]float64{
		{0.050383, 0.029803, 0.527975},
		{0.063536, 0.028426, 0.533124},
		{0.075353, 0.027206, 0.538007},
		{0.086222, 0.026125, 0.542658},
		{0.096379, 0.025165, 0.547103},
		{0.105980, 0.024309, 0.551368},
		{0.115124, 0.023556, 0.555468},
		{0.123903, 0.022878, 0.559423},
		{0.132381, 0.022258, 0.563250},
		{0.140603, 0.021687, 0.566959},
		{0.148607, 0.021154, 0.570562},
		{0.156421, 0.020651, 0.574065},
		{0.164070, 0.020171, 0.577478},
		{0.171574, 0.019706, 0.580806},
		{0.178950, 0.019252, 0.584054},
		{0.186213, 0.018803, 0.587228},
		{0.193374, 0.018354, 0.590330},
		{0.200445, 0.017902, 0.593364},
		{0.207435, 0.017442, 0.596333},
		{0.214350, 0.016973, 0.599239},
		{0.221197, 0.016497, 0.602083},
		{0.227983, 0.016007, 0.604867},
		{0.234715, 0.015502, 0.607592},
		{0.241396, 0.014979, 0.610259},
		{0.248032, 0.014439, 0.612868},
		{0.254627, 0.013882, 0.615419},
		{0.261183, 0.013308, 0.617911},
		{0.267703, 0.012716, 0.620346},
		{0.274191, 0.012109, 0.622722},
		{0.280648, 0.011488, 0.625038},
		{0.287076, 0.010855, 0.627295},
		{0.293478, 0.010213, 0.629490},
		{0.299855, 0.009561, 0.631624},
		{0.306210, 0.008902, 0.633694},
		{0.312543, 0.008239, 0.635700},
		{0.318856, 0.007576, 0.637640},
		{0.325150, 0.006915, 0.639512},
		{0.331426, 0.006261, 0.641316},
		{0.337683, 0.005618, 0.643049},
		{0.343925, 0.004991, 0.644710},
		{0.350150, 0.004382, 0.646298},
		{0.356359, 0.003798, 0.647810},
		{0.362553, 0.003243, 0.649245},
		{0.368733, 0.002724, 0.650601},
		{0.374897, 0.002245, 0.651876},
		{0.381047, 0.001814, 0.653068},
		{0.387183, 0.001434, 0.654177},
		{0.393304, 0.001114, 0.655199},
		{0.399411, 0.000859, 0.656133},
		{0.405503, 0.000678, 0.656977},
		{0.411580, 0.000577, 0.657730},
		{0.417642, 0.000564, 0.658390},
		{0.423689, 0.000646, 0.658956},
		{0.429719, 0.000831, 0.659425},
		{0.435734, 0.001127, 0.659797},
		{0.441732, 0.001540, 0.660069},
		{0.447714, 0.002080, 0.660240},
		{0.453677, 0.002755, 0.660310},
		{0.459623, 0.003574, 0.660277},
		{0.465550, 0.004545, 0.660139},
		{0.471457, 0.005678, 0.659897},
		{0.477344, 0.006980, 0.659549},
		{0.483210, 0.008460, 0.659095},
		{0.489055, 0.010127, 0.658534},
		{0.494877, 0.011990, 0.657865},
		{0.500678, 0.014055, 0.657088},
		{0.506454, 0.016333, 0.656202},
		{0.512206, 0.018833, 0.655209},
		{0.517933, 0.021563, 0.654109},
		{0.523633, 0.024532, 0.652901},
		{0.529306, 0.027747, 0.651586},
		{0.534952, 0.031217, 0.650165},
		{0.540570, 0.034950, 0.648640},
		{0.546157, 0.038954, 0.647010},
		{0.551715, 0.043136, 0.645277},
		{0.557243, 0.047331, 0.643443},
		{0.562738, 0.051545, 0.641509},
		{0.568201, 0.055778, 0.639477},
		{0.573632, 0.060028, 0.637349},
		{0.579029, 0.064296, 0.635126},
		{0.584391, 0.068579, 0.632812},
		{0.589719, 0.072878, 0.630408},
		{0.595011, 0.077190, 0.627917},
		{0.600266, 0.081516, 0.625342},
		{0.605485, 0.085854, 0.622686},
		{0.610667, 0.090204, 0.619951},
		{0.615812, 0.094564, 0.617140},
		{0.620919, 0.098934, 0.614257},
		{0.625987, 0.103312, 0.611305},
		{0.631017, 0.107699, 0.608287},
		{0.636008, 0.112092, 0.605205},
		{0.640959, 0.116492, 0.602065},
		{0.645872, 0.120898, 0.598867},
		{0.650746, 0.125309, 0.595617},
		{0.655580, 0.129725, 0.592317},
		{0.660374, 0.134144, 0.588971},
		{0.665129, 0.138566, 0.585582},
		{0.669845, 0.142992, 0.582154},
		{0.674522, 0.147419, 0.578688},
		{0.679160, 0.151848, 0.575189},
		{0.683758, 0.156278, 0.571660},
		{0.688318, 0.160709, 0.568103},
		{0.692840, 0.165141, 0.564522},
		{0.697324, 0.169573, 0.560919},
		{0.701769, 0.174005, 0.557296},
		{0.706178, 0.178437, 0.553657},
		{0.710549, 0.182868, 0.550004},
		{0.714883, 0.187299, 0.546338},
		{0.719181, 0.191729, 0.542663},
		{0.723444, 0.196158, 0.538981},
		{0.727670, 0.200586, 0.535293},
		{0.731862, 0.205013, 0.531601},
		{0.736019, 0.209439, 0.527908},
		{0.740143, 0.213864, 0.524216},
		{0.744232, 0.218288, 0.520524},
		{0.748289, 0.222711, 0.516834},
		{0.752312, 0.227133, 0.513149},
		{0.756304, 0.231555, 0.509468},
		{0.760264, 0.235976, 0.505794},
		{0.764193, 0.240396, 0.502126},
		{0.768090, 0.244817, 0.498465},
		{0.771958, 0.249237, 0.494813},
		{0.775796, 0.253658, 0.491171},
		{0.779604, 0.258078, 0.487539},
		{0.783383, 0.262500, 0.483918},
		{0.787133, 0.266922, 0.480307},
		{0.790855, 0.271345, 0.476706},
		{0.794549, 0.275770, 0.473117},
		{0.798216, 0.280197, 0.469538},
		{0.801855, 0.284626, 0.465971},
		{0.805467, 0.289057, 0.462415},
		{0.809052, 0.293491, 0.458870},
		{0.812612, 0.297928, 0.455338},
		{0.816144, 0.302368, 0.451816},
		{0.819651, 0.306812, 0.448306},
		{0.823132, 0.311261, 0.444806},
		{0.826588, 0.315714, 0.441316},
		{0.830018, 0.320172, 0.437836},
		{0.833422, 0.324635, 0.434366},
		{0.836801, 0.329105, 0.430905},
		{0.840155, 0.333580, 0.427455},
		{0.843484, 0.338062, 0.424013},
		{0.846788, 0.342551, 0.420579},
		{0.850066, 0.347048, 0.417153},
		{0.853319, 0.351553, 0.413734},
		{0.856547, 0.356066, 0.410322},
		{0.859750, 0.360588, 0.406917},
		{0.862927, 0.365119, 0.403519},
		{0.866078, 0.369660, 0.400126},
		{0.869203, 0.374212, 0.396738},
		{0.872303, 0.378774, 0.393355},
		{0.875376, 0.383347, 0.389976},
		{0.878423, 0.387932, 0.386600},
		{0.881443, 0.392529, 0.383229},
		{0.884436, 0.397139, 0.379860},
		{0.887402, 0.401762, 0.376494},
		{0.890340, 0.406398, 0.373130},
		{0.893250, 0.411048, 0.369768},
		{0.896131, 0.415712, 0.366407},
		{0.898984, 0.420392, 0.363047},
		{0.901807, 0.425087, 0.359688},
		{0.904601, 0.429797, 0.356329},
		{0.907365, 0.434524, 0.352970},
		{0.910098, 0.439268, 0.349610},
		{0.912800, 0.444029, 0.346251},
		{0.915471, 0.448807, 0.342890},
		{0.918109, 0.453603, 0.339529},
		{0.920714, 0.458417, 0.336166},
		{0.923287, 0.463251, 0.332801},
		{0.925825, 0.468103, 0.329435},
		{0.928329, 0.472975, 0.326067},
		{0.930798, 0.477867, 0.322697},
		{0.933232, 0.482780, 0.319325},
		{0.935630, 0.487712, 0.315952},
		{0.937990, 0.492667, 0.312575},
		{0.940313, 0.497642, 0.309197},
		{0.942598, 0.502639, 0.305816},
		{0.944844, 0.507658, 0.302433},
		{0.947051, 0.512699, 0.299049},
		{0.949217, 0.517763, 0.295662},
		{0.951344, 0.522850, 0.292275},
		{0.953428, 0.527960, 0.288883},
		{0.955470, 0.533093, 0.285490},
		{0.957469, 0.538250, 0.282096},
		{0.959424, 0.543431, 0.278701},
		{0.961336, 0.548636, 0.275305},
		{0.963203, 0.553865, 0.271909},
		{0.965024, 0.559118, 0.268513},
		{0.966798, 0.564396, 0.265118},
		{0.968526, 0.569700, 0.261721},
		{0.970205, 0.575028, 0.258325},
		{0.971835, 0.580382, 0.254931},
		{0.973416, 0.585761, 0.251540},
		{0.974947, 0.591165, 0.248151},
		{0.976428, 0.596595, 0.244767},
		{0.977856, 0.602051, 0.241387},
		{0.979233, 0.607532, 0.238013},
		{0.980556, 0.613039, 0.234646},
		{0.981826, 0.618572, 0.231287},
		{0.983041, 0.624131, 0.227937},
		{0.984199, 0.629718, 0.224595},
		{0.985301, 0.635330, 0.221265},
		{0.986345, 0.640969, 0.217948},
		{0.987332, 0.646633, 0.214648},
		{0.988260, 0.652325, 0.211364},
		{0.989128, 0.658043, 0.208100},
		{0.989935, 0.663787, 0.204859},
		{0.990681, 0.669558, 0.201642},
		{0.991365, 0.675355, 0.198453},
		{0.991985, 0.681179, 0.195295},
		{0.992541, 0.687030, 0.192170},
		{0.993032, 0.692907, 0.189084},
		{0.993456, 0.698810, 0.186041},
		{0.993814, 0.704741, 0.183043},
		{0.994103, 0.710698, 0.180097},
		{0.994324, 0.716681, 0.177208},
		{0.994474, 0.722691, 0.174381},
		{0.994553, 0.728728, 0.171622},
		{0.994561, 0.734791, 0.168938},
		{0.994495, 0.740880, 0.166335},
		{0.994355, 0.746995, 0.163821},
		{0.994141, 0.753137, 0.161404},
		{0.993851, 0.759304, 0.159092},
		{0.993482, 0.765499, 0.156891},
		{0.993033, 0.771720, 0.154808},
		{0.992505, 0.777967, 0.152855},
		{0.991897, 0.784239, 0.151042},
		{0.991209, 0.790537, 0.149377},
		{0.990439, 0.796859, 0.147870},
		{0.989587, 0.803205, 0.146529},
		{0.988648, 0.809579, 0.145357},
		{0.987621, 0.815978, 0.144363},
		{0.986509, 0.822401, 0.143557},
		{0.985314, 0.828846, 0.142945},
		{0.984031, 0.835315, 0.142528},
		{0.982653, 0.841812, 0.142303},
		{0.981190, 0.848329, 0.142279},
		{0.979644, 0.854866, 0.142453},
		{0.977995, 0.861432, 0.142808},
		{0.976265, 0.868016, 0.143351},
		{0.974443, 0.874622, 0.144061},
		{0.972530, 0.881250, 0.144923},
		{0.970533, 0.887896, 0.145919},
		{0.968443, 0.894564, 0.147014},
		{0.966271, 0.901249, 0.148180},
		{0.964021, 0.907950, 0.149370},
		{0.961681, 0.914672, 0.150520},
		{0.959276, 0.921407, 0.151566},
		{0.956808, 0.928152, 0.152409},
		{0.954287, 0.934908, 0.152921},
		{0.951726, 0.941671, 0.152925},
		{0.949151, 0.948435, 0.152178},
		{0.946602, 0.955190, 0.150328},
		{0.944152, 0.961916, 0.146861},
		{0.941896, 0.968590, 0.140956},
		{0.940015, 0.975158, 0.131326},
	}
)

func get_index(x float64, N int) int {
	n := int(x * float64(N-1))

	if n < 0 {
		n = 0
	} else if n >= N {
		n = N - 1
	}

	return n
}

func parula(x float64) [3]float64 {
	N := len(parula_table)
	n := get_index(x, N)

	return parula_table[n]
}

func plasma(x float64) [3]float64 {
	N := len(parula_table)
	n := get_index(x, N)

	return plasma_table[n]
}

type colormap_type int

const (
	Parula colormap_type = iota
	Plasma

	MaxTypes
)

func colormap(x float64, which colormap_type) [3]float64 {
	switch which {
	case Parula:
		return parula(x)
	case Plasma:
		return plasma(x)
	default:
		return [3]float64{0, 0, 0}
	}
}

func make_bands(x float64, numBands int) float64 {
	if numBands == 0 {
		return x
	}

	n := int(x * float64(numBands))

	if n < 0 {
		n = 0
	} else if n > numBands {
		n = numBands
	}

	return float64(n) / float64(numBands)
}
func colormapN(x float64, which colormap_type, N int) [3]float64 {
	x = make_bands(x, N)
	switch which {
	case Parula:
		return parula(x)
	case Plasma:
		return plasma(x)
	default:
		return [3]float64{0, 0, 0}
	}
}

func (t colormap_type) ToString() string {
	switch t {
	case Parula:
		return "parula"
	case Plasma:
		return "plasma"
	default:
		return "..."
	}
}

func map_type_from_string(name *string) colormap_type {

	lower := strings.ToLower(*name)

	switch lower {
	case "parula":
		return Parula
	case "plasma":
		return Plasma
	default:
		fmt.Println(name, " is not a valid colormap")
		return Parula
	}
}

func print_map_types() {
	fmt.Printf("Supported colormaps: \n")
	for n := 0; n < int(MaxTypes); n++ {
		fmt.Printf(" - %v\n", colormap_type(n).ToString())
	}
}

func colormap_info() {
	fmt.Printf("Parula: %v points\n", len(parula_table))
	fmt.Printf("Plasma: %v points\n", len(plasma_table))
}