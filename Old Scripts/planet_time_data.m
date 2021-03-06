

target_planets = {'mercury' 'venus' 'jupiter' 'saturn' 'uranus' 'neptune'};

ti_min_planets = [2016 2016 2016 2016  2017 2020 ;  % year
                  1    6    6    9     1     1   ;  % month
                  1    1    1    1     1     1  ];  % day

ti_max_planets = [2017 2018 2018 2018 2019 2021 ;
                  1    3    6    10   1   10    ;
                  1    1    1    1    1   1    ];
    
tf_min_planets = [2016 2016 2018 2019 2021 2031 ;
                  4    12   6    4    4    1    ;
                  1    1    1    1    1    1   ];

tf_max_planets = [2017 2019 2022 2025 2035 2050 ;
                  3    1    3    3    3    6
                  1    1    1    1    1    1   ];           

v_inf_planets =   [7.0  3.0  9.1  11.5 12.1 12.5];




target_smallbodies = {'pluto' 'toutatis' 'adonis' 'hermes' 'icarus' 'geographos' 'castalia' 'eros' 'itokawa' 'apollo' 'nereus' 'oljato' 'orpheus'  'dyonisus' 'nyx' 'minos' 'florence' 'pan' '1996FG3' 'apophis' '1999KW4' 'aten' 'ganymed' 'quetzalcoatl' '2000SG344' '1989UQ' 'SG286' 'SG344 ' '2001CC21' 'SG3441999YB' '1994CN2'};

ti_min_smallbodies = [2019 2020 2020 2020 2023 2018 2025 2020 2020 2019 2019 2019 2020 2023 2021 2020 2024 2021 2021 2018 2024 2021 2023 2018 2025 2029 2019 2020 2020 2021 2021 ;
                      10   5    5    7    1    5    2    7    5    1    5    9    6    1    5    8    3    11   12   6    6    12   4    12   10   8    3    1    10   10   4    ;
                      1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1   ];

ti_max_smallbodies = [2021 2022 2023 2023 2025 2020 2027 2023 2022 2022 2022 2021 2022 2024 2024 2022 2027 2024 2024 2021 2027 2024 2026 2022 2030 2033 2021 2022 2023 2024 2025 ;
                      5    6    2    4    1    8    8    6    8    8    8    11   8    11   5    10   1    1    1    12   1    6    1    6    1    1    12   12   6    10   2    ;
                      1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1   ];
       
tf_min_smallbodies = [2031 2021 2021 2021 2023 2019 2025 2021 2021 2020 2020 2021 2020 2023 2022 2021 2024 2021 2021 2018 2024 2022 2024 2019 2026 2030 2020 2021 2021 2022 2021 ;
                      1    6    6    6    4    6    11   1    6    6    6    1    11   11   6    6    11   12   12   9    12   6    6    11   2    5    1    1    9    6    9    ;        
                      1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1   ];
tf_max_smallbodies = [2050 2028 2026 2025 2027 2021 2029 2026 2024 2024 2024 2025 2023 2027 2026 2024 2029 2025 2025 2022 2027 2025 2029 2025 2031 2034 2024 2025 2025 2026 2027 ;
                      7    9    2    12   6    12   2    2    8    1    8    1    8    2    5    2    6    1    8    12   5    12   4    1    3    3    1    1    7    10   3    ;
                      1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1   ];
                  
v_inf_smallbodies =   [15.9 6.8  7.1  6.8  7.8  5.2  6.2  4.1  3.9  6.0  4.2  8.0  2.6  6.3  4.5  4.1  5.8  6.1  3.1  2.3  7.8  5.7  6.4  6.8  0.5  2.5  3.7  2.1  2.1  2.9  4.8 ];

target = cell(1,78); 
ti_min = zeros(3,78);
ti_max = zeros(3,78);
tf_min = zeros(3,78);
tf_max = zeros(3,78);
v_inf = zeros(1,78);

positions = [1,2,5,6,7,8,9,14,15,16,18,19,36,70,71,12,13,17,20,21,22,23,35,37,40,55,73,74,76,78,60,38,26,31,39,41,42];
not_positions = [3 4 10 11 24 25 27:30 32:34 43:54 56:59 61:69 72 75 77];

target(1:37) = horzcat(target_planets,target_smallbodies);
v_inf(1:37) = [v_inf_planets v_inf_smallbodies];
ti_min(1:3,1:37) = [ti_min_planets ti_min_smallbodies];
ti_max(1:3,1:37) = [ti_max_planets ti_max_smallbodies];
tf_min(1:3,1:37) = [tf_min_planets tf_min_smallbodies];
tf_max(1:3,1:37) = [tf_max_planets tf_max_smallbodies];


target(positions) = target(1:37);
target(not_positions) = cell(1,1);
v_inf(positions) = v_inf(1:37);
v_inf(not_positions) = 0;
ti_min(:,positions) = ti_min(:,1:37);
ti_max(:,positions) = ti_max(:,1:37);
tf_min(:,positions) = tf_min(:,1:37);
tf_max(:,positions) = tf_max(:,1:37);

for i = 1:41
ti_min(1:3,not_positions(i)) = zeros(3,1);
ti_max(1:3,not_positions(i)) = zeros(3,1);
tf_min(1:3,not_positions(i)) = zeros(3,1);
tf_max(1:3,not_positions(i)) = zeros(3,1);
end

clear('ti_min_planets','ti_max_planets','tf_max_planets','tf_min_planets','v_inf_planets','target_smallbodies','target_planets','ti_min_smallbodies','ti_max_smallbodies','tf_min_smallbodies','tf_max_smallbodies','v_inf_smallbodies','positions','not_positions','i');



