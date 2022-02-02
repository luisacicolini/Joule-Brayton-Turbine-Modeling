clear all
close all
load Data_struct

                                 %   Temperatures
figure('Name','Temperature measurements')

errorbar(Data.Run1.N1_1_min_,Data.Run1.Tt3_K_,Data.Error1.Tt3_K_,'Color', [0 0.447 0.741])
hold on
errorbar(Data.Run1.N1_1_min_,Data.Run1.Tt4_K_,Data.Error1.Tt4_K_,'r')
errorbar(Data.Run1.N1_1_min_,Data.Run1.T9_K_,Data.Error1.T9_K_,'Color', [0.9290 0.6940 0.1250])

errorbar(Data.Run2.N1_1_min_,Data.Run2.Tt3_K_,Data.Error2.Tt3_K_,'Color', [0 0.447 0.741],'LineStyle','--')
errorbar(Data.Run2.N1_1_min_,Data.Run2.Tt4_K_,Data.Error2.Tt4_K_,'r','LineStyle','--')
errorbar(Data.Run2.N1_1_min_,Data.Run2.T9_K_,Data.Error2.T9_K_,'Color', [0.9290 0.6940 0.1250],'LineStyle','--')
title('Temperature measurement','FontSize',14)
legend('Tt3 (run 1)','Tt4 (run 1)','Tt9 (run 1)','Tt3 (run 2)','Tt4 (run 2)','Tt9 (run 2)','Location','northwest')
xlabel('Engine speed [1/min]');
ylabel('Temperature [K]');
xlim([0 110000])
ylim([0 1250])
ax = gca;
ax.XAxis.Exponent = 3;
grid on
grid minor



                                %   Volumetric flow - AIR
figure('Name','Volumetric flow measurements - air')

errorbar(Data.Run1.N1_1_min_,Data.Run1.Q_L__m_3_s_,Data.Error1.Q_L__m_3_s_,'Color', [0 0.447 0.741])
hold on
errorbar(Data.Run1.N1_1_min_,Data.Run2.Q_L__m_3_s_,Data.Error2.Q_L__m_3_s_,'Color', 'r')

title('Volumetric flow measurement','FontSize',14)
legend('Q air (run 1)','Q air (run 2)','Location','northwest')
xlabel('Engine speed [1/min]');
ylabel('Volumetric flow [m^3/s]');
xlim([0 110000])
%ylim([0 1250])
ax = gca;
ax.XAxis.Exponent = 3;
grid on
grid minor

                             %   Volumetric flow - FUEL
figure('Name','Volumetric flow measurements - fuel')

errorbar(Data.Run1.N1_1_min_,Data.Run1.Q_Br__m_3_s_,Data.Error1.Q_Br__m_3_s_,'Color', [0 0.447 0.741])
hold on
errorbar(Data.Run1.N1_1_min_,Data.Run2.Q_Br__m_3_s_,Data.Error2.Q_Br__m_3_s_,'Color', 'r')

title('Volumetric flow measurement','FontSize',14)
legend('Q fuel (run 1)','Q fuel (run 2)','Location','northwest')
xlabel('Engine speed [1/min]');
ylabel('Volumetric flow [m^3/s]');
xlim([0 110000])
%ylim([0 1250])
ax = gca;
ax.XAxis.Exponent = 3;
grid on
grid minor

                             %   Force measurement
figure('Name','Thrust measurements')

errorbar(Data.Run1.N1_1_min_,Data.Run1.F_N_,Data.Error1.F_N_,'Color', [0 0.447 0.741])
hold on
errorbar(Data.Run1.N1_1_min_,Data.Run2.F_N_,Data.Error2.F_N_,'Color', 'r')

title('Thrust measurement','FontSize',14)
legend('Thrust (run 1)','Thrust (run 2)','Location','northwest')
xlabel('Engine speed [1/min]');
ylabel('Thrust [N]');
xlim([0 110000])
%ylim([0 1250])
ax = gca;
ax.XAxis.Exponent = 3;
grid on
grid minor

                            %   Pressure measurement
figure('Name','Pressure measurements')

errorbar(Data.Run1.N1_1_min_,Data.Run1.p3_Pa_,Data.Error1.p3_Pa_,'Color', [0 0.447 0.741])
hold on
errorbar(Data.Run1.N1_1_min_,Data.Run2.p3_Pa_,Data.Error2.p3_Pa_,'Color', 'r')

title('Pressure measurement','FontSize',14)
legend('Pt3 (excess to ambient) (run 1)','Pt3 (excess to ambient) (run 2)','Location','northwest')
xlabel('Engine speed [1/min]');
ylabel('Pressure [Pa]');
xlim([0 110000])
%ylim([0 1250])
ax = gca;
ax.XAxis.Exponent = 3;
ax.YAxis.Exponent = 3;
grid on
grid minor

disp('Done')