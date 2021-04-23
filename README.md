## Requirements
CARLsim 4, ECJ-28, MATLAB 2020a

## ca1_sub_snn_final
Compile with ``make``, and use ``./launchCARLsimECJ.sh`` to launch evolutionary runs. Network fitness values and evolved hyper-parameters saved to ``out.stat``.

## ca1_sub_snn_test_final
Use ``cat params/params_ca1.csv | ./cogmap-snn`` to run evolved individual networks. Network activity and input information saved to disk.

## matlab_analysis_scripts
MATLAB scripts for data analysis and visualization. 
