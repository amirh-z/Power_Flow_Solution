clc;
clear;
close all;

% bus data
%      ----bus----   Voltage   Angle    -Load-     -------Generator-------
%      bus_i  type   Mag.      Degree   Pd   Qd    PG    QG   Min Q  Max Q
bus = [
        1      1     1.04        0      0    0     72.3  0   -300    300 ;
        2      2     1.025       0      0    0     163   0   -300    300 ;
        3      2     1.025       0      0    0     85    0   -300    300 ;
        4      0     1           0      0    0     0     0    0      0   ;
        5      0     1           0      90   30    0     0    0      0   ;
        6      0     1           0      0    0     0     0    0      0   ;
        7      0     1           0      100  35    0     0    0      0   ;
        8      0     1           0      0    0     0     0    0      0   ;
        9      0     1           0      125  50    0     0    0      0   
      ];
  
  
% branch data
%          fbus    tbus       RL(pu)    XL(pu) 
branch = [
            1       4         0         0.0576 ;
            4       5         0.017     0.092  ;
            5       6         0.039     0.17   ;
            3       6         0         0.0586 ;
            6       7         0.0119    0.1008 ;
            7       8         0.0085    0.072  ;
            8       2         0         0.0625 ;
            8       9         0.032     0.161  ;
            9       4         0.01      0.085  
         ];
  

        