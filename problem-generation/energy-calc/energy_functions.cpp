/*

    This file contains the energy functions
    for the HP and Miyazawa-Jernigan models.

*/

#include <iostream>

#include "energy_functions.h"


double hp_energy(char aa1, char aa2) {

    if (aa1 == 'H' and aa2 == 'H') {
        return 1.0;
    } else {
        return 0.0;
    }

}

double mj_energy(char aa1, char aa2) {

    if ((aa1=='C' and aa2=='C') or (aa2=='C' and aa1=='C')) return 5.440000;
    if ((aa1=='M' and aa2=='C') or (aa2=='M' and aa1=='C')) return 4.990000;
    if ((aa1=='M' and aa2=='M') or (aa2=='M' and aa1=='M')) return 5.460000;
    if ((aa1=='F' and aa2=='C') or (aa2=='F' and aa1=='C')) return 5.800000;
    if ((aa1=='F' and aa2=='M') or (aa2=='F' and aa1=='M')) return 6.560000;
    if ((aa1=='F' and aa2=='F') or (aa2=='F' and aa1=='F')) return 7.260000;
    if ((aa1=='I' and aa2=='C') or (aa2=='I' and aa1=='C')) return 5.500000;
    if ((aa1=='I' and aa2=='M') or (aa2=='I' and aa1=='M')) return 6.020000;
    if ((aa1=='I' and aa2=='F') or (aa2=='I' and aa1=='F')) return 6.840000;
    if ((aa1=='I' and aa2=='I') or (aa2=='I' and aa1=='I')) return 6.540000;
    if ((aa1=='L' and aa2=='C') or (aa2=='L' and aa1=='C')) return 5.830000;
    if ((aa1=='L' and aa2=='M') or (aa2=='L' and aa1=='M')) return 6.410000;
    if ((aa1=='L' and aa2=='F') or (aa2=='L' and aa1=='F')) return 7.280000;
    if ((aa1=='L' and aa2=='I') or (aa2=='L' and aa1=='I')) return 7.040000;
    if ((aa1=='L' and aa2=='L') or (aa2=='L' and aa1=='L')) return 7.370000;
    if ((aa1=='V' and aa2=='C') or (aa2=='V' and aa1=='C')) return 4.960000;
    if ((aa1=='V' and aa2=='M') or (aa2=='V' and aa1=='M')) return 5.320000;
    if ((aa1=='V' and aa2=='F') or (aa2=='V' and aa1=='F')) return 6.290000;
    if ((aa1=='V' and aa2=='I') or (aa2=='V' and aa1=='I')) return 6.050000;
    if ((aa1=='V' and aa2=='L') or (aa2=='V' and aa1=='L')) return 6.480000;
    if ((aa1=='V' and aa2=='V') or (aa2=='V' and aa1=='V')) return 5.520000;
    if ((aa1=='W' and aa2=='C') or (aa2=='W' and aa1=='C')) return 4.950000;
    if ((aa1=='W' and aa2=='M') or (aa2=='W' and aa1=='M')) return 5.550000;
    if ((aa1=='W' and aa2=='F') or (aa2=='W' and aa1=='F')) return 6.160000;
    if ((aa1=='W' and aa2=='I') or (aa2=='W' and aa1=='I')) return 5.780000;
    if ((aa1=='W' and aa2=='L') or (aa2=='W' and aa1=='L')) return 6.140000;
    if ((aa1=='W' and aa2=='V') or (aa2=='W' and aa1=='V')) return 5.180000;
    if ((aa1=='W' and aa2=='W') or (aa2=='W' and aa1=='W')) return 5.060000;
    if ((aa1=='Y' and aa2=='C') or (aa2=='Y' and aa1=='C')) return 4.160000;
    if ((aa1=='Y' and aa2=='M') or (aa2=='Y' and aa1=='M')) return 4.910000;
    if ((aa1=='Y' and aa2=='F') or (aa2=='Y' and aa1=='F')) return 5.660000;
    if ((aa1=='Y' and aa2=='I') or (aa2=='Y' and aa1=='I')) return 5.250000;
    if ((aa1=='Y' and aa2=='L') or (aa2=='Y' and aa1=='L')) return 5.670000;
    if ((aa1=='Y' and aa2=='V') or (aa2=='Y' and aa1=='V')) return 4.620000;
    if ((aa1=='Y' and aa2=='W') or (aa2=='Y' and aa1=='W')) return 4.660000;
    if ((aa1=='Y' and aa2=='Y') or (aa2=='Y' and aa1=='Y')) return 4.170000;
    if ((aa1=='A' and aa2=='C') or (aa2=='A' and aa1=='C')) return 3.570000;
    if ((aa1=='A' and aa2=='M') or (aa2=='A' and aa1=='M')) return 3.940000;
    if ((aa1=='A' and aa2=='F') or (aa2=='A' and aa1=='F')) return 4.810000;
    if ((aa1=='A' and aa2=='I') or (aa2=='A' and aa1=='I')) return 4.580000;
    if ((aa1=='A' and aa2=='L') or (aa2=='A' and aa1=='L')) return 4.910000;
    if ((aa1=='A' and aa2=='V') or (aa2=='A' and aa1=='V')) return 4.040000;
    if ((aa1=='A' and aa2=='W') or (aa2=='A' and aa1=='W')) return 3.820000;
    if ((aa1=='A' and aa2=='Y') or (aa2=='A' and aa1=='Y')) return 3.360000;
    if ((aa1=='A' and aa2=='A') or (aa2=='A' and aa1=='A')) return 2.720000;
    if ((aa1=='G' and aa2=='C') or (aa2=='G' and aa1=='C')) return 3.160000;
    if ((aa1=='G' and aa2=='M') or (aa2=='G' and aa1=='M')) return 3.390000;
    if ((aa1=='G' and aa2=='F') or (aa2=='G' and aa1=='F')) return 4.130000;
    if ((aa1=='G' and aa2=='I') or (aa2=='G' and aa1=='I')) return 3.780000;
    if ((aa1=='G' and aa2=='L') or (aa2=='G' and aa1=='L')) return 4.160000;
    if ((aa1=='G' and aa2=='V') or (aa2=='G' and aa1=='V')) return 3.380000;
    if ((aa1=='G' and aa2=='W') or (aa2=='G' and aa1=='W')) return 3.420000;
    if ((aa1=='G' and aa2=='Y') or (aa2=='G' and aa1=='Y')) return 3.010000;
    if ((aa1=='G' and aa2=='A') or (aa2=='G' and aa1=='A')) return 2.310000;
    if ((aa1=='G' and aa2=='G') or (aa2=='G' and aa1=='G')) return 2.240000;
    if ((aa1=='T' and aa2=='C') or (aa2=='T' and aa1=='C')) return 3.110000;
    if ((aa1=='T' and aa2=='M') or (aa2=='T' and aa1=='M')) return 3.510000;
    if ((aa1=='T' and aa2=='F') or (aa2=='T' and aa1=='F')) return 4.280000;
    if ((aa1=='T' and aa2=='I') or (aa2=='T' and aa1=='I')) return 4.030000;
    if ((aa1=='T' and aa2=='L') or (aa2=='T' and aa1=='L')) return 4.340000;
    if ((aa1=='T' and aa2=='V') or (aa2=='T' and aa1=='V')) return 3.460000;
    if ((aa1=='T' and aa2=='W') or (aa2=='T' and aa1=='W')) return 3.220000;
    if ((aa1=='T' and aa2=='Y') or (aa2=='T' and aa1=='Y')) return 3.010000;
    if ((aa1=='T' and aa2=='A') or (aa2=='T' and aa1=='A')) return 2.320000;
    if ((aa1=='T' and aa2=='G') or (aa2=='T' and aa1=='G')) return 2.080000;
    if ((aa1=='T' and aa2=='T') or (aa2=='T' and aa1=='T')) return 2.120000;
    if ((aa1=='S' and aa2=='C') or (aa2=='S' and aa1=='C')) return 2.860000;
    if ((aa1=='S' and aa2=='M') or (aa2=='S' and aa1=='M')) return 3.030000;
    if ((aa1=='S' and aa2=='F') or (aa2=='S' and aa1=='F')) return 4.020000;
    if ((aa1=='S' and aa2=='I') or (aa2=='S' and aa1=='I')) return 3.520000;
    if ((aa1=='S' and aa2=='L') or (aa2=='S' and aa1=='L')) return 3.920000;
    if ((aa1=='S' and aa2=='V') or (aa2=='S' and aa1=='V')) return 3.050000;
    if ((aa1=='S' and aa2=='W') or (aa2=='S' and aa1=='W')) return 2.990000;
    if ((aa1=='S' and aa2=='Y') or (aa2=='S' and aa1=='Y')) return 2.780000;
    if ((aa1=='S' and aa2=='A') or (aa2=='S' and aa1=='A')) return 2.010000;
    if ((aa1=='S' and aa2=='G') or (aa2=='S' and aa1=='G')) return 1.820000;
    if ((aa1=='S' and aa2=='T') or (aa2=='S' and aa1=='T')) return 1.960000;
    if ((aa1=='S' and aa2=='S') or (aa2=='S' and aa1=='S')) return 1.670000;
    if ((aa1=='N' and aa2=='C') or (aa2=='N' and aa1=='C')) return 2.590000;
    if ((aa1=='N' and aa2=='M') or (aa2=='N' and aa1=='M')) return 2.950000;
    if ((aa1=='N' and aa2=='F') or (aa2=='N' and aa1=='F')) return 3.750000;
    if ((aa1=='N' and aa2=='I') or (aa2=='N' and aa1=='I')) return 3.240000;
    if ((aa1=='N' and aa2=='L') or (aa2=='N' and aa1=='L')) return 3.740000;
    if ((aa1=='N' and aa2=='V') or (aa2=='N' and aa1=='V')) return 2.830000;
    if ((aa1=='N' and aa2=='W') or (aa2=='N' and aa1=='W')) return 3.070000;
    if ((aa1=='N' and aa2=='Y') or (aa2=='N' and aa1=='Y')) return 2.760000;
    if ((aa1=='N' and aa2=='A') or (aa2=='N' and aa1=='A')) return 1.840000;
    if ((aa1=='N' and aa2=='G') or (aa2=='N' and aa1=='G')) return 1.740000;
    if ((aa1=='N' and aa2=='T') or (aa2=='N' and aa1=='T')) return 1.880000;
    if ((aa1=='N' and aa2=='S') or (aa2=='N' and aa1=='S')) return 1.580000;
    if ((aa1=='N' and aa2=='N') or (aa2=='N' and aa1=='N')) return 1.680000;
    if ((aa1=='Q' and aa2=='C') or (aa2=='Q' and aa1=='C')) return 2.850000;
    if ((aa1=='Q' and aa2=='M') or (aa2=='Q' and aa1=='M')) return 3.300000;
    if ((aa1=='Q' and aa2=='F') or (aa2=='Q' and aa1=='F')) return 4.100000;
    if ((aa1=='Q' and aa2=='I') or (aa2=='Q' and aa1=='I')) return 3.670000;
    if ((aa1=='Q' and aa2=='L') or (aa2=='Q' and aa1=='L')) return 4.040000;
    if ((aa1=='Q' and aa2=='V') or (aa2=='Q' and aa1=='V')) return 3.070000;
    if ((aa1=='Q' and aa2=='W') or (aa2=='Q' and aa1=='W')) return 3.110000;
    if ((aa1=='Q' and aa2=='Y') or (aa2=='Q' and aa1=='Y')) return 2.970000;
    if ((aa1=='Q' and aa2=='A') or (aa2=='Q' and aa1=='A')) return 1.890000;
    if ((aa1=='Q' and aa2=='G') or (aa2=='Q' and aa1=='G')) return 1.660000;
    if ((aa1=='Q' and aa2=='T') or (aa2=='Q' and aa1=='T')) return 1.900000;
    if ((aa1=='Q' and aa2=='S') or (aa2=='Q' and aa1=='S')) return 1.490000;
    if ((aa1=='Q' and aa2=='N') or (aa2=='Q' and aa1=='N')) return 1.710000;
    if ((aa1=='Q' and aa2=='Q') or (aa2=='Q' and aa1=='Q')) return 1.540000;
    if ((aa1=='D' and aa2=='C') or (aa2=='D' and aa1=='C')) return 2.410000;
    if ((aa1=='D' and aa2=='M') or (aa2=='D' and aa1=='M')) return 2.570000;
    if ((aa1=='D' and aa2=='F') or (aa2=='D' and aa1=='F')) return 3.480000;
    if ((aa1=='D' and aa2=='I') or (aa2=='D' and aa1=='I')) return 3.170000;
    if ((aa1=='D' and aa2=='L') or (aa2=='D' and aa1=='L')) return 3.400000;
    if ((aa1=='D' and aa2=='V') or (aa2=='D' and aa1=='V')) return 2.480000;
    if ((aa1=='D' and aa2=='W') or (aa2=='D' and aa1=='W')) return 2.840000;
    if ((aa1=='D' and aa2=='Y') or (aa2=='D' and aa1=='Y')) return 2.760000;
    if ((aa1=='D' and aa2=='A') or (aa2=='D' and aa1=='A')) return 1.700000;
    if ((aa1=='D' and aa2=='G') or (aa2=='D' and aa1=='G')) return 1.590000;
    if ((aa1=='D' and aa2=='T') or (aa2=='D' and aa1=='T')) return 1.800000;
    if ((aa1=='D' and aa2=='S') or (aa2=='D' and aa1=='S')) return 1.630000;
    if ((aa1=='D' and aa2=='N') or (aa2=='D' and aa1=='N')) return 1.680000;
    if ((aa1=='D' and aa2=='Q') or (aa2=='D' and aa1=='Q')) return 1.460000;
    if ((aa1=='D' and aa2=='D') or (aa2=='D' and aa1=='D')) return 1.210000;
    if ((aa1=='E' and aa2=='C') or (aa2=='E' and aa1=='C')) return 2.270000;
    if ((aa1=='E' and aa2=='M') or (aa2=='E' and aa1=='M')) return 2.890000;
    if ((aa1=='E' and aa2=='F') or (aa2=='E' and aa1=='F')) return 3.560000;
    if ((aa1=='E' and aa2=='I') or (aa2=='E' and aa1=='I')) return 3.270000;
    if ((aa1=='E' and aa2=='L') or (aa2=='E' and aa1=='L')) return 3.590000;
    if ((aa1=='E' and aa2=='V') or (aa2=='E' and aa1=='V')) return 2.670000;
    if ((aa1=='E' and aa2=='W') or (aa2=='E' and aa1=='W')) return 2.990000;
    if ((aa1=='E' and aa2=='Y') or (aa2=='E' and aa1=='Y')) return 2.790000;
    if ((aa1=='E' and aa2=='A') or (aa2=='E' and aa1=='A')) return 1.510000;
    if ((aa1=='E' and aa2=='G') or (aa2=='E' and aa1=='G')) return 1.220000;
    if ((aa1=='E' and aa2=='T') or (aa2=='E' and aa1=='T')) return 1.740000;
    if ((aa1=='E' and aa2=='S') or (aa2=='E' and aa1=='S')) return 1.480000;
    if ((aa1=='E' and aa2=='N') or (aa2=='E' and aa1=='N')) return 1.510000;
    if ((aa1=='E' and aa2=='Q') or (aa2=='E' and aa1=='Q')) return 1.420000;
    if ((aa1=='E' and aa2=='D') or (aa2=='E' and aa1=='D')) return 1.020000;
    if ((aa1=='E' and aa2=='E') or (aa2=='E' and aa1=='E')) return 0.910000;
    if ((aa1=='H' and aa2=='C') or (aa2=='H' and aa1=='C')) return 3.600000;
    if ((aa1=='H' and aa2=='M') or (aa2=='H' and aa1=='M')) return 3.980000;
    if ((aa1=='H' and aa2=='F') or (aa2=='H' and aa1=='F')) return 4.770000;
    if ((aa1=='H' and aa2=='I') or (aa2=='H' and aa1=='I')) return 4.140000;
    if ((aa1=='H' and aa2=='L') or (aa2=='H' and aa1=='L')) return 4.540000;
    if ((aa1=='H' and aa2=='V') or (aa2=='H' and aa1=='V')) return 3.580000;
    if ((aa1=='H' and aa2=='W') or (aa2=='H' and aa1=='W')) return 3.980000;
    if ((aa1=='H' and aa2=='Y') or (aa2=='H' and aa1=='Y')) return 3.520000;
    if ((aa1=='H' and aa2=='A') or (aa2=='H' and aa1=='A')) return 2.410000;
    if ((aa1=='H' and aa2=='G') or (aa2=='H' and aa1=='G')) return 2.150000;
    if ((aa1=='H' and aa2=='T') or (aa2=='H' and aa1=='T')) return 2.420000;
    if ((aa1=='H' and aa2=='S') or (aa2=='H' and aa1=='S')) return 2.110000;
    if ((aa1=='H' and aa2=='N') or (aa2=='H' and aa1=='N')) return 2.080000;
    if ((aa1=='H' and aa2=='Q') or (aa2=='H' and aa1=='Q')) return 1.980000;
    if ((aa1=='H' and aa2=='D') or (aa2=='H' and aa1=='D')) return 2.320000;
    if ((aa1=='H' and aa2=='E') or (aa2=='H' and aa1=='E')) return 2.150000;
    if ((aa1=='H' and aa2=='H') or (aa2=='H' and aa1=='H')) return 3.050000;
    if ((aa1=='R' and aa2=='C') or (aa2=='R' and aa1=='C')) return 2.570000;
    if ((aa1=='R' and aa2=='M') or (aa2=='R' and aa1=='M')) return 3.120000;
    if ((aa1=='R' and aa2=='F') or (aa2=='R' and aa1=='F')) return 3.980000;
    if ((aa1=='R' and aa2=='I') or (aa2=='R' and aa1=='I')) return 3.630000;
    if ((aa1=='R' and aa2=='L') or (aa2=='R' and aa1=='L')) return 4.030000;
    if ((aa1=='R' and aa2=='V') or (aa2=='R' and aa1=='V')) return 3.070000;
    if ((aa1=='R' and aa2=='W') or (aa2=='R' and aa1=='W')) return 3.410000;
    if ((aa1=='R' and aa2=='Y') or (aa2=='R' and aa1=='Y')) return 3.160000;
    if ((aa1=='R' and aa2=='A') or (aa2=='R' and aa1=='A')) return 1.830000;
    if ((aa1=='R' and aa2=='G') or (aa2=='R' and aa1=='G')) return 1.720000;
    if ((aa1=='R' and aa2=='T') or (aa2=='R' and aa1=='T')) return 1.900000;
    if ((aa1=='R' and aa2=='S') or (aa2=='R' and aa1=='S')) return 1.620000;
    if ((aa1=='R' and aa2=='N') or (aa2=='R' and aa1=='N')) return 1.640000;
    if ((aa1=='R' and aa2=='Q') or (aa2=='R' and aa1=='Q')) return 1.800000;
    if ((aa1=='R' and aa2=='D') or (aa2=='R' and aa1=='D')) return 2.290000;
    if ((aa1=='R' and aa2=='E') or (aa2=='R' and aa1=='E')) return 2.270000;
    if ((aa1=='R' and aa2=='H') or (aa2=='R' and aa1=='H')) return 2.160000;
    if ((aa1=='R' and aa2=='R') or (aa2=='R' and aa1=='R')) return 1.550000;
    if ((aa1=='K' and aa2=='C') or (aa2=='K' and aa1=='C')) return 1.950000;
    if ((aa1=='K' and aa2=='M') or (aa2=='K' and aa1=='M')) return 2.480000;
    if ((aa1=='K' and aa2=='F') or (aa2=='K' and aa1=='F')) return 3.360000;
    if ((aa1=='K' and aa2=='I') or (aa2=='K' and aa1=='I')) return 3.010000;
    if ((aa1=='K' and aa2=='L') or (aa2=='K' and aa1=='L')) return 3.370000;
    if ((aa1=='K' and aa2=='V') or (aa2=='K' and aa1=='V')) return 2.490000;
    if ((aa1=='K' and aa2=='W') or (aa2=='K' and aa1=='W')) return 2.690000;
    if ((aa1=='K' and aa2=='Y') or (aa2=='K' and aa1=='Y')) return 2.600000;
    if ((aa1=='K' and aa2=='A') or (aa2=='K' and aa1=='A')) return 1.310000;
    if ((aa1=='K' and aa2=='G') or (aa2=='K' and aa1=='G')) return 1.150000;
    if ((aa1=='K' and aa2=='T') or (aa2=='K' and aa1=='T')) return 1.310000;
    if ((aa1=='K' and aa2=='S') or (aa2=='K' and aa1=='S')) return 1.050000;
    if ((aa1=='K' and aa2=='N') or (aa2=='K' and aa1=='N')) return 1.210000;
    if ((aa1=='K' and aa2=='Q') or (aa2=='K' and aa1=='Q')) return 1.290000;
    if ((aa1=='K' and aa2=='D') or (aa2=='K' and aa1=='D')) return 1.680000;
    if ((aa1=='K' and aa2=='E') or (aa2=='K' and aa1=='E')) return 1.800000;
    if ((aa1=='K' and aa2=='H') or (aa2=='K' and aa1=='H')) return 1.350000;
    if ((aa1=='K' and aa2=='R') or (aa2=='K' and aa1=='R')) return 0.590000;
    if ((aa1=='K' and aa2=='K') or (aa2=='K' and aa1=='K')) return 0.120000;
    if ((aa1=='P' and aa2=='C') or (aa2=='P' and aa1=='C')) return 3.070000;
    if ((aa1=='P' and aa2=='M') or (aa2=='P' and aa1=='M')) return 3.450000;
    if ((aa1=='P' and aa2=='F') or (aa2=='P' and aa1=='F')) return 4.250000;
    if ((aa1=='P' and aa2=='I') or (aa2=='P' and aa1=='I')) return 3.760000;
    if ((aa1=='P' and aa2=='L') or (aa2=='P' and aa1=='L')) return 4.200000;
    if ((aa1=='P' and aa2=='V') or (aa2=='P' and aa1=='V')) return 3.320000;
    if ((aa1=='P' and aa2=='W') or (aa2=='P' and aa1=='W')) return 3.730000;
    if ((aa1=='P' and aa2=='Y') or (aa2=='P' and aa1=='Y')) return 3.190000;
    if ((aa1=='P' and aa2=='A') or (aa2=='P' and aa1=='A')) return 2.030000;
    if ((aa1=='P' and aa2=='G') or (aa2=='P' and aa1=='G')) return 1.870000;
    if ((aa1=='P' and aa2=='T') or (aa2=='P' and aa1=='T')) return 1.900000;
    if ((aa1=='P' and aa2=='S') or (aa2=='P' and aa1=='S')) return 1.570000;
    if ((aa1=='P' and aa2=='N') or (aa2=='P' and aa1=='N')) return 1.530000;
    if ((aa1=='P' and aa2=='Q') or (aa2=='P' and aa1=='Q')) return 1.730000;
    if ((aa1=='P' and aa2=='D') or (aa2=='P' and aa1=='D')) return 1.330000;
    if ((aa1=='P' and aa2=='E') or (aa2=='P' and aa1=='E')) return 1.260000;
    if ((aa1=='P' and aa2=='H') or (aa2=='P' and aa1=='H')) return 2.250000;
    if ((aa1=='P' and aa2=='R') or (aa2=='P' and aa1=='R')) return 1.700000;
    if ((aa1=='P' and aa2=='K') or (aa2=='P' and aa1=='K')) return 0.970000;
    if ((aa1=='P' and aa2=='P') or (aa2=='P' and aa1=='P')) return 1.750000;

    std::cout << "Warning: no value returned in MJ energy for " << aa1 << ", " << aa2 << ".\n";

    return 0.0;
 
}

