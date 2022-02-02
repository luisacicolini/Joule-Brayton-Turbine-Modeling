%close all
clear all
% TURBOJET ENGINE ANALYSIS
%and other stuff 
%ideal gas entropy calculation

format long
load Data_struct.mat
% Data.Run1([6,7,8,9],:) = [];
% Data.Run2([6,7,8,9],:) = [];
% Data.Error1([6,7,8,9],:) = [];
% Data.Error2([6,7,8,9],:) = [];

for i = 1:5

%total temperature
Tt=struct('step0',0,'step2',0,'step3',0,'step4',0,'step9',0);

%total pressure
pt=struct('step0',0,'step2',0,'step3',0,'step4',0,'step9',0);

%static pressure
p=struct('step0',0,'step2',0,'step3',0,'step4',0,'step9',0);

%static temperature
T=struct('step0',0,'step2',0,'step3',0,'step4',0,'step9',0);

%pressure ratio
pratio=struct('inlet',1,'combustion',0.92,'turbine',0); % updated 04/12

%Mach number
Ma=struct('step0',0,'step2',0,'step3',0,'step4',0,'step9',0);

%speed
c=struct('step0',0,'step2',0,'step3',0,'step4',0,'step9',0);

%specific heat capacities
cp=struct('step0',1003.5,'step2',1003.5,'step3',1005,'step4',1142,'step9',1150);
cv=struct('step0',716.5,'step2',716.5,'step3',718,'step4',855,'step9',860);

%isentropic temperatures
Tis=struct('step0',0,'step2',0,'step3',0,'step4',0,'step9',0);

%density values for air
rho=struct('step0',1.225,'step2',1.225,'step3',71.225,'step4',1.225,'step9',1.225);

%entropy
entropy=struct('step0',0,'step2',0,'step3',0,'step4',0,'step9',0);

%calculated with ideal gas assumption, but compressible: ideal gas is more
%dependent on the conditions of the fluid 

%AMBIENT 

T.step0=273.15+5.6; %[K]
p.step0=96900; %[Pa]

%MEASURED INPUT

Tt.step3 = Data.Run1.Tt3_K_(i); %K
Tt.step4 = Data.Run1.Tt4_K_(i); %K
Tt.step9 = Data.Run1.T9_K_(i); %K
vflow_fuel = Data.Run1.Q_Br__m_3_s_(i); %m^3/s
pt.step3=Data.Run1.p3_Pa_(i)+p.step0; %Pa
vflow_air = Data.Run1.Q_L__m_3_s_(i); %m^3/s

%OTHER VARIABLES:

density_fuel=0.79*1000;
%r_comp=0.076;
%height_comp=0.01;
%rg_turb=0.0703;
%rn_turb=0.05;
Hi=43.45e06; %[J/kg]
Rgas=8.314; %[J/mol/K] [constant]
Mm_air=0.029; %[kg/mol]
Mm_fuel=0.162; %[kg/mol]
Rair=Rgas/Mm_air; %[J/kg/K] [relative R]
Rk=Rgas/Mm_fuel; %[J/kg/K] [relative R]
etamech=0.99; %ASSUMPTION
A2=pi*(0.042^2)/4; %[m^2]
A9=pi*(0.050^2)/4; %[m^2]

%ASSIGNING MEASURED VALUES TO CORRESPONDING VARIABLES:

% CALCULATIONS

                            %ADIABATIC INLET 0:2
                            
%MEASURED: Tt.step0,pt.step0,p.step0;
%total values
Tt.step0=T.step0;
pt.step0=p.step0;
Tt.step2=Tt.step0;
pt.step2=pt.step0*pratio.inlet;
entropy.step0=0; %ASSUMPTION entropy at the initial state = 0 J/kgK, for reference
[cp.step0,cv.step0]=ccalc(273.15,'polynomial',1,0); %only for entropy 'zero point' definition

%speed and static values
rho.step0=density(p.step0,T.step0,Mm_air);

c.step2=vflow_air/A2; %m/s
[cp.step2,cv.step2]=ccalc(Tt.step2,'polynomial',1,0);
T.step2=temp(Tt.step2,c.step2,cp.step2); %[K]
p.step2=pt.step2/(1+((c.step2)^2/(2*Rair*T.step2))); %[Pa]

%mass flow

rho.step2=density(p.step2,T.step2,Mm_air); %[kg/m^3]
mflow_air=vflow_air*rho.step2; %[kg/s]
mflow_fuel=vflow_fuel*density_fuel; %[kg/s]
airfuel_ratio=mflow_air/mflow_fuel; %[m^3/s]

%isothermal proces 0t->2t, entropy for an ideal gas
%entropy.step2=entropy.step0+(cp.step2+cp.step0)/2*log(273.15/Tt.step2)-Rair*log(95724.4/pt.step2); %[J/K] 95724.4 is pressure at 0C (x axis definition)
%entropy.step2=entropy.step0+(cp.step2+cp.step0)/2*log(Tt.step2/Tt.step0)-Rair*log(pt.step2/pt.step0);
entropy.step2=entropy.step0+(cp.step2+cp.step0)/2*log(Tt.step2/273.15)-Rair*log(pt.step2/95724.4);

                            %COMPRESSOR 2:3
                            
%MEASURED: Tt.step3,pt.step3;
air=1;
fuel=0;
[cp.step3,cv.step3]=ccalc(Tt.step3,'polynomial',air,fuel); %J/kgK
k=cp.step3/cv.step3; %[-]
Tis.step3=tisentropic(Tt.step2,pt.step3,pt.step2,k); %[K]

%efficiencies:
Teffis=(Tis.step3-Tt.step2)/(Tt.step3-Tt.step2);
%should be isos, but for the real case we have a real Tt
%therefore we can calculate the entropy
%always assuming ideal gas
entropy.step3=entropy.step2+(cp.step2+cp.step3)/2*log(Tt.step3/Tt.step2)-Rair*log(pt.step3/pt.step2);
                                %COMBUSTION 3:4

R=Rgas/(Mm_air*air+Mm_fuel*fuel); %[J/kg/K]
%combustion should be isobaric, but we have the pressure loss
%therefore we can calculate the entropy, always under the ideal gas
%assumption
air=mflow_air/(mflow_air+mflow_fuel); %[-] air fraction
fuel=mflow_fuel/(mflow_air+mflow_fuel); %[-] fuel fraction
pt.step4=pt.step3*pratio.combustion; %[Pa]
[cp.step4,cv.step4]=ccalc(Tt.step4,'polynomial',air,fuel);
cp34=(cp.step3+cp.step4)/2; %average as an polynomial is fine
entropy.step4=entropy.step3+cp34*log(Tt.step4/Tt.step3)-R*log(pt.step4/pt.step3);

% %efficiencies - is it needed?
% T.step3=temp(Tt.step3,%%speed needed,cp.step3);
% T.step4=temp(Tt.step4,%%speed needed,cp.step4);
% effcomb=((mflow_air+mflow_fuel)*cp34*T.step4-mflow_air*cp34*T.step3)/Hi*mflow_fuel
effcomb(i)=((mflow_air+mflow_fuel)*cp34*Tt.step4-mflow_air*cp34*Tt.step3)/(Hi*mflow_fuel);


                                %TURBINE 4:9

%MEASURED: Tt.step4,Tt.step9,p.step9;

% rho.step4 = density(T.step4,p.step4, );
% c.step4=(mflow_air+mflow_fuel)/area(rg_turb,rn_turb,0,'axial')/rho.step4;
k=cp.step4/cv.step3; %[-]

%T.step4=temp(Tt.step4,c.step4,cp.step4); % TO BE UPDATED WITH AREA
% NOT SURE p.step4=p.step3;
%Ma.step4=mach(c.step4,k,T.step4,R);
%pt.step4=ptotal(p.step4,k,Ma.step4);

%step9: impossible to measure T.step9, since we always have stagnation
%points
%}

% ADIABATIC NOZZLE: short time period, too short for heat exchange
% PRESSURE DROP NOZZLE: boudary layer turbolences + pressure losses (there
% is not in free flow) && frictions with the wall (form drag)
% Turbolence in the nozzle: evaluate Reynolds (e.g. dimensional analysis with differentials)
% but we can assume no pressure loss
% possible calculation of critical Reynolds?

%ideal gas: always static temp for density
[cp.step9,cv.step9] = ccalc(Tt.step9,'polynomial',air,fuel);

% GETTING c9, T9, pt9

p.step9=p.step0;

A = ((mflow_air+mflow_fuel)^2*R^2)/(2*cp.step9*A9^2*p.step9^2);
B = 1;
C = -Tt.step9;
e = [A B C];
T.step9=roots(e);

if T.step9(1) > 0
    T.step9 = T.step9(1);
else
    T.step9 = T.step9(2);
end

rho.step9 = density(p.step9,T.step9,(Mm_air*air+Mm_fuel*fuel));
c.step9 = (mflow_air+mflow_fuel)/(rho.step9*A9);
pt.step9 = p.step9+rho.step9/2*c.step9^2;

cp49=(cp.step4+cp.step9)/2;
entropy.step9=entropy.step4+cp49*log(Tt.step9/Tt.step4)-R*log(pt.step9/pt.step4);

Tis.step9=Tt.step4*(pt.step9/pt.step4)^(R/cp.step9);

TeffisT = (Tt.step4-Tt.step9)/(Tt.step4-Tis.step9);

%THRUST COMPARISON

thrust1(i)=(mflow_air+mflow_fuel)*(c.step9-c.step0);

%ENTROPY&TEMPERATURE DEFINITION
c1(i,:)=(struct2array(c));
entropy1(i,:) = (struct2array(entropy));
temperature1(i,:) = (struct2array(Tt));
pressure1(i,:) = (struct2array(pt));
temperature2(i,:) = (struct2array(T));
pressure2(i,:) = (struct2array(p));
air1(i) = air;
fuel1(i) = fuel;

end

datameasured = fopen('data1.txt','w');
fprintf(datameasured, 'RESULTS - MEASURED DATA \n');
fprintf(datameasured,' speed \n');
for r=1:5
    for c=1:5
        fprintf(datameasured,string(c1(r,c)));
        fprintf(datameasured,'-');
    end
    fprintf(datameasured,'\n');
end
fprintf(datameasured,'\n entropy \n');
for r=1:5
    for c=1:5
        fprintf(datameasured,string(entropy1(r,c)));
        fprintf(datameasured,'-');
    end
    fprintf(datameasured,'\n');
end
fprintf(datameasured,'\n total temperature \n');
for r=1:5
    for c=1:5
        fprintf(datameasured,string(temperature1(r,c)));
        fprintf(datameasured,'-');
    end
    fprintf(datameasured,'\n');
end
fprintf(datameasured,'\n static temperature \n');
for r=1:5
    for c=1:5
        fprintf(datameasured,string(temperature2(r,c)));
        fprintf(datameasured,'-');
    end
    fprintf(datameasured,'\n');
end
fprintf(datameasured,'\n total pressure \n');
for r=1:5
    for c=1:5
        fprintf(datameasured,string(pressure1(r,c)));
        fprintf(datameasured,'-');
    end
    fprintf(datameasured,'\n');
end
fprintf(datameasured,'\n static pressure \n');
for r=1:5
    for c=1:5
        fprintf(datameasured,string(pressure2(r,c)));
        fprintf(datameasured,'-');
    end
    fprintf(datameasured,'\n');
end
fprintf(datameasured,'\n air \n');
for r=1:5
        fprintf(datameasured,string(air1(r)));
        fprintf(datameasured,'-');
end
fprintf(datameasured,'\n fuel \n');
for r=1:5
        fprintf(datameasured,string(fuel1(r)));
        fprintf(datameasured,'-');
end
fprintf(datameasured,'\n compressor efficiency \n');
fprintf(datameasured, string(Teffis));
fprintf(datameasured,'\n turbine efficiency \n');
fprintf(datameasured, string(TeffisT));
fclose(datameasured);

% figure('Name','Thrust comparison')
% 
% plot(Data.Run1.N1_1_min_,thrust1,'--o')
% hold on
% errorbar(Data.Run1.N1_1_min_,Data.Run1.F_N_,Data.Error1.F_N_,'-o')
% %plot(Data.Run2.N1_1_min_,thrust2,'--o')
% %errorbar(Data.Run2.N1_1_min_,Data.Run2.F_N_,Data.Error2.F_N_,'-o')
% legend('Calculated thrust - run 1 [N]','Measured thrust - run 1 [N]','Location','northwest')
% 
% figure('Name','Cycle graph');
% hold on
% for k = 1:5
% plot(entropy1(k,2:end),temperature1(k,2:end),'-o')
% plot(polyshape(fliplr(entropy1(k,2:end)),fliplr(temperature1(k,2:end))))
% end
% legend('IDLE','40k rpm','60k rpm','80k rpm','100k rpm','Location','northwest')


%cycleplot(temperature1,entropy1,air1,fuel1,'first run ');
% get indicators straight (efficiencies, etc)