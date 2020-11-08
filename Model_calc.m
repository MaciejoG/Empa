clear all
clc
close all
load('Results_model_v13_tc_varrying_Tdhw')

COP_avg % 1) 3.345575319134551   2) 3.129036076291251
P_hidden % 1) 3.057292599058989e+02   2) 3.801857894216212e+02
P_gw % 1) 2.265239729111218e+02    2) 0

start_day=30;

figure(1) % temperatures in year
plot(time_table(1,:),Troom_table(1,:),'g',time_table(1,:),Tgw_table(1,:),'r',time_table(1,:),Tamb0_table(1,:),'b',time_table(1,:),Tamb1_table(1,:),'k',time_table(1,:),Tamb2_table(1,:),'c',time_table(1,:),Tceiling_table(1,:),'y');
axis([0 time_simulation -10 Thot+1])
title('Temperatures in time');
xlabel('Time t [s]');
ylabel('Temperatures T [C]');
legend('Troom','Tgw','Tsurface','T1.5','T3.5','Tceiling')
grid on
print('Temperatures_in_time','-dpd')

figure(2) % COP in year
plot(time_table(1,:),COP(1,:),'g');
title('COP(t)');
xlabel('Time t [s]');
ylabel('Temperature [C]');
grid on
print('COP','-dpd')

figure(3) % Te in year
plot(time_table(1,:),Te(1,1:time_simulation),'g');
title('Tevap(t)');
xlabel('Time t [s]');
ylabel('Temperature [C]');
grid on
print('Tevap','-dpd')

figure(4) % dTair in year
plot(time_table(1,:),dTair(1,1:time_simulation),'g');
title('dTair(t)');
xlabel('Time t [s]');
ylabel('Temperature [C]');
grid on
print('dTair','-dpd')

figure(5) % Tc in year
plot(time_table(1,:),Tc(1,1:time_simulation),'g');
title('Tcond(t)');
xlabel('Time t [s]');
ylabel('Temperature [C]');
grid on
print('Tcond','-dpd')

figure(6) % TDHW in year
plot(time_table(1,:),T_DHW(1,1:time_simulation),'g');
title('TDHW(t)');
xlabel('Time t [s]');
ylabel('Temperature [C]');
grid on
print('TDHW','-dpd')

figure(7) % temperatures in day 31
plot(time_table(1,(start_day*24*3600):((start_day+1)*24*3600)),Troom_table(1,(start_day*24*3600):((start_day+1)*24*3600)),'g',time_table(1,(start_day*24*3600):((start_day+1)*24*3600)),Tgw_table(1,(start_day*24*3600):((start_day+1)*24*3600)),'r',time_table(1,(start_day*24*3600):((start_day+1)*24*3600)),Tamb0_table(1,(start_day*24*3600):((start_day+1)*24*3600)),'b',time_table(1,(start_day*24*3600):((start_day+1)*24*3600)),Tamb1_table(1,(start_day*24*3600):((start_day+1)*24*3600)),'k',time_table(1,(start_day*24*3600):((start_day+1)*24*3600)),Tamb2_table(1,(start_day*24*3600):((start_day+1)*24*3600)),'c',time_table(1,(start_day*24*3600):((start_day+1)*24*3600)),Tceiling_table(1,(start_day*24*3600):((start_day+1)*24*3600)),'y');
axis([(start_day*24*3600) ((start_day+1)*24*3600) -10 Thot+1])
title('Temperatures in time during day number 31');
xlabel('Time t [s]');
ylabel('Temperatures T [C]');
legend('Troom','Tgw','Tsurface','T1.5','T3.5','Tceiling')
grid on
print('Temperatures_in_time_during_day_number_31','-dpd')

figure(8) % COP day 31
plot(time_table(1,(start_day*24*3600):((start_day+1)*24*3600)),COP(1,(start_day*24*3600):((start_day+1)*24*3600)),'g');
title('COP(t)');
xlabel('Time t [s]');
ylabel('COP [-]');
grid on
print('COP_day_31','-dpd')

figure(9) % Te day 31
plot(time_table(1,(start_day*24*3600):((start_day+1)*24*3600)),Te(1,(start_day*24*3600):((start_day+1)*24*3600)),'g');
title('Tevap(t)');
xlabel('Time t [s]');
ylabel('Temperature [C]');
grid on
print('Tevap_day_31','-dpd')

figure(10) % dTair day 31
plot(time_table(1,(start_day*24*3600):((start_day+1)*24*3600)),dTair(1,(start_day*24*3600):((start_day+1)*24*3600)),'g');
title('dTair(t)');
xlabel('Time t [s]');
ylabel('Temperature [C]');
grid on
print('dTair_day_31','-dpd')

figure(11) % Tc day 31
plot(time_table(1,(start_day*24*3600):((start_day+1)*24*3600)),Tc(1,(start_day*24*3600):((start_day+1)*24*3600)),'g');
title('Tcond(t)');
xlabel('Time t [s]');
ylabel('Temperature [C]');
grid on
print('Tcond_day_31','-dpd')

figure(12) % TDHW day 31
plot(time_table(1,(start_day*24*3600):((start_day+1)*24*3600)),T_DHW(1,(start_day*24*3600):((start_day+1)*24*3600)),'g');
title('TDHW(t)');
xlabel('Time t [s]');
ylabel('Temperature [C]');
grid on
print('TDHW_day_31','-dpd')
