function [gwf, dt] = safe_example_gwf()
% function [gwf, dt] = safe_example_gwf()
% Waveform with some frequency matching by Filip Szczepankiewicz.

%% STE
dt  = 1e-3; % ms

gwf = 0.08 * [ % T/m
     0         0         0
    -0.2005    0.9334    0.3029
    -0.2050    0.9324    0.3031
    -0.2146    0.9302    0.3032
    -0.2313    0.9263    0.3030
    -0.2589    0.9193    0.3019
    -0.3059    0.9060    0.2980
    -0.3892    0.8767    0.2883
    -0.3850    0.7147    0.3234
    -0.3687    0.5255    0.3653
    -0.3509    0.3241    0.4070
    -0.3323    0.1166    0.4457
    -0.3136   -0.0906    0.4783
    -0.2956   -0.2913    0.5019
    -0.2790   -0.4793    0.5139
    -0.2642   -0.6491    0.5118
    -0.2518   -0.7957    0.4939
    -0.2350   -0.8722    0.4329
    -0.2187   -0.9111    0.3541
    -0.2063   -0.9409    0.2747
    -0.1977   -0.9627    0.1933
    -0.1938   -0.9768    0.1080
    -0.1967   -0.9820    0.0159
    -0.2114   -0.9751   -0.0883
    -0.2292   -0.9219   -0.2150
    -0.2299   -0.8091   -0.3561
    -0.2290   -0.6748   -0.5011
    -0.2253   -0.5239   -0.6460
    -0.2178   -0.3620   -0.7868
    -0.2056   -0.1948   -0.9194
    -0.1391   -0.0473   -0.9908
    -0.0476    0.0607   -0.9987
     0.0215    0.1452   -0.9909
     0.0725    0.2136   -0.9759
     0.1114    0.2709   -0.9579
     0.1426    0.3204   -0.9383
     0.1690    0.3641   -0.9177
     0         0         0
     0         0         0
     0         0         0
     0         0         0
     0         0         0
     0         0         0
     0         0         0
    -0.3734   -0.1768    0.9125
    -0.3825   -0.2310    0.8965
    -0.3919   -0.2895    0.8752
    -0.4015   -0.3543    0.8465
    -0.4108   -0.4290    0.8065
    -0.4182   -0.5202    0.7469
    -0.4178   -0.6423    0.6451
    -0.3855   -0.8173    0.4321
    -0.3110   -0.9418    0.1401
    -0.2526   -0.9669   -0.0674
    -0.2100   -0.9541   -0.2213
    -0.1766   -0.9227   -0.3474
    -0.1491   -0.8788   -0.4570
    -0.1258   -0.8239   -0.5555
    -0.1056   -0.7583   -0.6459
    -0.0882   -0.6809   -0.7293
    -0.0734   -0.5900   -0.8061
    -0.0615   -0.4830   -0.8753
    -0.0533   -0.3556   -0.9349
    -0.0506   -0.2005   -0.9801
    -0.0575   -0.0019   -1.0000
    -0.0909    0.2976   -0.9521
    -0.3027    0.9509   -0.0860
    -0.2737    0.9610   -0.0692
    -0.2524    0.9675   -0.0596
    -0.2364    0.9719   -0.0533
    -0.2245    0.9749   -0.0490
    -0.2158    0.9770   -0.0459
    -0.2097    0.9785   -0.0439
    -0.2058    0.9794   -0.0426
    -0.2039    0.9798   -0.0420
     0         0         0
    ];


