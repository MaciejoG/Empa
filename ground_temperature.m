% A function determining 

clear all
close all
clc

% Importing 'outside air temperature' data from excel file
filename = 'EMPA_Temperatures.xlsx';
sheet = 1;
xlRange = 'B5:B8764';
ylRange = 'L5:L8764';
xdata = xlsread(filename,sheet,xlRange);
ydata = xlsread(filename,sheet,ylRange);

a=3600*1.4/(1920*1850); % a [m2/h] = 3600 * k_gravel/(dens_gravel*cp_gravel), soruce: Waermeatlas - 1.1 Hausenklima und meteorologische parametres..
tyear=365*24;

% Properties used to calculate temperatures in the ground
parts=10; % How many parts is the Tamb(x) curve divided to
for l=1:1:(parts+1) % loop calculating different depth, which will be later used to calculate an avarage temperature in the ground
    x1(l)=0.15+0.2+(2/parts)*(l-1); % Depth(wall) (ceiling(insulation + ceiling(concrete) + gravel(0-2m))
    damper_x1(l)=exp(-x1(l)*sqrt(pi/(a*tyear)));
end
% damper_x1_avg=sum(damper_x1)/(parts+1);
x2=0.2+0.15+2+0.2+0.15+0.5; % depth(ground) (Depth(wall) + ground(concrete) + ground(insulation)+gravel)
damper_x2=exp(-x2*sqrt(pi/(a*tyear)));

figure(1)
plot(x1(:),damper_x1(:))
title('Damping factor as a funtion of depth');
xlabel('Depth [m]');
ylabel('Damping factor [-]');
grid on
print('Damping_factor_as_a_funtion_of_depth','-dpd')

t=1:2*tyear;
TA_amb=10;
theta_m=10;% mean(T_amb);
T_amb=theta_m*ones(1,length(t))+TA_amb*cos(2*pi/tyear*t); % original function at 0 m depth
T3m_x2=theta_m*ones(1,length(t))+TA_amb*damper_x2*cos(2*pi/tyear*t-x2*sqrt(pi/(a*tyear))*ones(1,length(t)));

for t=1:1:2*tyear
    for l=1:1:(parts+1) % calculation of temperatures in the ground at different depths
        T3m_x1(t,l)=theta_m+TA_amb*exp(-x1(l)*sqrt(pi/(a*tyear)))*cos(2*pi/tyear*t-x1(l)*sqrt(pi/(a*tyear)));
    end
%     T3m_x1_avg(t)=sum(T3m_x1)/(parts+1);
    time_table(t)=t;
end

figure(2)
t=1:2*tyear;
plot(t,T_amb,'r')
hold on
plot(time_table(:),T3m_x1(:,1),'b',time_table(:),T3m_x1(:,2),'y',time_table(:),T3m_x1(:,3),'c',time_table(:),T3m_x1(:,4),'m')
hold on
t=1:2*tyear;
plot(t,T3m_x2,'g')
hold off
legend('x1','x2','x3','x4','x5','x6>>x5')
title('Temperatature variations at different depths in time');
xlabel('Time t [s]');
ylabel('Temperatures T [C]');
grid on
print('Temperatature_variations_at_different_depths_in_time','-dpd')


[N,J]=max(T3m_x2(:));
[O,K]=max(T_amb(:));
for i=1:1:(parts+1)
    [M(i),I(i)]=max(T3m_x1(:,i));
    difference_hours_x1(i)=tyear+(I(i)-K);
end
difference_hours_x2=tyear-(K-J);

% DATA TO BE EXPORTED TO MAIN MODEL
for t=1:1:(365*24)
    for i=1:1:(parts+1)
        if t>0&&t<=difference_hours_x1(i)
            ydata_x1(t,i)=damper_x1(i)*ydata(365*24-difference_hours_x1(i)+t);
        else
            ydata_x1(t,i)=damper_x1(i)*ydata(t-difference_hours_x1(i));
        end
    end
    ydata_x1_avg(t)=sum(ydata_x1(t,:))/(parts+1);
end

for t=1:1:(365*24)
    if t>0&&t<=difference_hours_x2
        ydata_x2(t)=damper_x2*ydata(365*24-difference_hours_x2+t);
    else
        ydata_x2(t)=damper_x2*ydata(t-difference_hours_x2);
    end
    time_table2(t)=t;
% ydata_x2(i)=ydata(i-difference_hours_x2);
end

figure(3)
plot(time_table2(:),ydata(:),'r',time_table2(:),ydata_x1_avg(:),'b',time_table2(:),ydata_x2(:),'g')
title('Ground and surface temperatures in time');
xlabel('Time t [s]');
ylabel('Temperatures T [C]');
legend('x0m','x1.5m','x3m')
grid on
print('Ground_and_surface_temperatures','-dpd')

save('ydata_damped_shifted','ydata_x1_avg','ydata_x2')




