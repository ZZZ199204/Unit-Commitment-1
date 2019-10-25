global UnitStates;
UnitStates=zeros(8,24);
MIN_UP_DOWN_TIME_FLAG       = 1;
RAMP_UP_DOWN_FLAG           = 0;
N_PRED                      = 1;
COMPLETE_ENUMERATION_FLAG   = 1;
DETAIL_PRINT_FLAG           = 0;
DISPATCH_METHOD             = 3;
RESERVE_FLAG                = 0;
START_UP_COST_METHOD        = 1;
gen_data = [...
      1        25     80      10440           213.00       350            2.00          4             2           -5          150                4                50             75          NaN         NaN           NaN               0               NaN 
      2        60    250       9000           585.62       400            2.00          5             3           +8          170                5                80            120          NaN         NaN           NaN               0               NaN 
      3        75    300       8730           684.74      1100            2.00          5             4           +8          500                5               100            150          NaN         NaN           NaN               0               NaN 
      4        20     60      11900           252.00      0.02            2.00          1             1           -6            0                0                80            120          NaN         NaN           NaN               0               NaN 
      5        25     80      10440           213.00       350            2.00          4             2           -5          150                4                50             75          NaN         NaN           NaN               0               NaN 
      6        60    250       9000           585.62       400            2.00          5             3           +8          170                5                80            120          NaN         NaN           NaN               0               NaN 
      7        75    300       8730           684.74      1100            2.00          5             4           +8          500                5               100            150          NaN         NaN           NaN               0               NaN 
      8        20     60      11900           252.00      0.02            2.00          1             1           -6            0                0                80            120          NaN         NaN           NaN               0               NaN 
      ];
DEMAND = [900;1060;1200;1080;800;560;580;1000;900;1060;1200;1080;800;560;580;1000;900;1060;1200;1080;800;560;580;1000];
K_RES_UP = 0.00;
K_RES_DN = 0.00;