clear all
close all
format shortG
%Run up and run down averaged data used. Remove '_avg' and ranges to use
%original data

Run1 = readtable('Messdaten(1).xlsx','Sheet','Average','Range','B2:N11');
Data.Run1 = removevars(Run1,{'Tt3__C_','Tt4__C_','T9__C_','Q_Br__ml_min_','Q_L__l_s_'});
Error1 = readtable('Messdaten(1).xlsx','Sheet','Average','Range','R2:AA11');
Data.Error1 = removevars(Error1,{'Q_Br__ml_min_','Q_L__l_s_'});

Run2 = readtable('Messdaten(2).xlsx','Sheet','Average','Range','B2:N11');
Data.Run2 = removevars(Run2,{'Tt3__C_','Tt4__C_','T9__C_','Q_Br__ml_min_','Q_L__l_s_'});
Error2 = readtable('Messdaten(1).xlsx','Sheet','Average','Range','R2:AA11');
Data.Error2 = removevars(Error2,{'Q_Br__ml_min_','Q_L__l_s_'});

save('Data_struct.mat','Data')
disp('Done')