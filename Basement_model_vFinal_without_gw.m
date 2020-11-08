clear all
close all
clc
tic; % simulation time measuring

load('ydata_damped_shifted') % loading temperature data (at various depths) calculated in 'ground_temperature.m'

x=[14 17 19.8 22.7 25.1 27.6 30 34 38 42 46 50];
y=[42 44.7 47.2 48.8 50 50.8 51.5 52.1 52.5 52.9 53.1 53.2];
p=polyfit(x,y,4); % curve fitted to the diagram presenting isentropic_efficiency=f(temperature_lift), where temp_lift=t_condensation-t_evaporation; it is used later in the script

% Importing 'outside air temperature' data from excel file
filename = 'EMPA_Temperatures.xlsx';
sheet = 1;
xlRange = 'B5:B8764';
ylRange = 'L5:L8764';
xdata = xlsread(filename,sheet,xlRange);
ydata = xlsread(filename,sheet,ylRange);

m_gw_4_person=142*4; % l/day,  source: graph from Luca, 'Wasserverbrauch im Haushalt pro einwohner und Tag'
shower=0.253; % source: graph from Luca, 'Wasserverbrauch im Haushalt pro einwohner und Tag'
kitchen=0.155; % source: graph from Luca, 'Wasserverbrauch im Haushalt pro einwohner und Tag'
washbasin=0.1130; % source: graph from Luca, 'Wasserverbrauch im Haushalt pro einwohner und Tag'
dishwasher=0.021; % source: graph from Luca, 'Wasserverbrauch im Haushalt pro einwohner und Tag'
washingmashine=0.12; % source: graph from Luca, 'Wasserverbrauch im Haushalt pro einwohner und Tag'
temp_shower=41; % source: 'Recommenden Code of Practise for Safe Water Temperature'
temp_kitchen=48; % source: 'Recommenden Code of Practise for Safe Water Temperature'
temp_washbasin=41; % source: 'Recommenden Code of Practise for Safe Water Temperature'
temp_dishwasher=57; % source: http://products.geappliances.com/appliance/gea-support-search-content?contentId=18924
temp_washingmashine=43.3; % for 'warm water', source: https://www.thespruce.com/laundry-and-water-temperature-1900646
m_c_shower=shower*m_gw_4_person*(60-temp_shower)/50; % mass of cold water (10 C) used in a shower
m_h_shower=shower*m_gw_4_person-m_c_shower; % mass of hot water (temp_shower) used in shower
m_c_kitchen=kitchen*m_gw_4_person*(60-temp_kitchen)/50;
m_h_kitchen=kitchen*m_gw_4_person-m_c_kitchen;
m_c_washbasin=washbasin*m_gw_4_person*(60-temp_washbasin)/50;
m_h_washbasin=washbasin*m_gw_4_person-m_c_washbasin;
m_c_dishwasher=dishwasher*m_gw_4_person*(60-temp_dishwasher)/50;
m_h_dishwasher=dishwasher*m_gw_4_person-m_c_dishwasher;
m_c_washingmashine=washingmashine*m_gw_4_person*(60-temp_washingmashine)/50;
m_h_washingmashine=washingmashine*m_gw_4_person-m_c_washingmashine;
m_gw=(shower+kitchen+washbasin+dishwasher+washingmashine)*m_gw_4_person; % total mass of grey water produced by 4 ppl daily
m_DHW=m_h_shower+m_h_kitchen+m_h_washbasin+m_h_dishwasher+m_h_washingmashine; % totalt mass of DHW used by 4 ppl daily
c_water=4200; % specific heat of water, J/kgK
Thot=(shower*temp_shower+kitchen*temp_kitchen+washbasin*temp_washbasin+dishwasher*temp_dishwasher+washingmashine*temp_washingmashine)/(shower+kitchen+washbasin+dishwasher+washingmashine)*0.55; % temperature of grey water DHW usage, C

% Initial conditions: 
% ROOM
initial_conditions(1)=5.557040889749037; % Initial temperature in the room
% WALLS
initial_conditions(2)=6.189783418450240; % Initial avarage temperature of wall's first layer
initial_conditions(3)=4.403774528630980; % Initial avarage temperature of wall's second layer
initial_conditions(4)=2.523400034528492; % Initial avarage temperature of gravel -> walls
% GROUND
initial_conditions(5)=6.161551919752201; % Initial avarage temperature of ground's first layer
initial_conditions(6)=4.436240311986332; % Initial avarage temperature of ground's second layer
initial_conditions(7)=2.643737767331162; % Initial avarage temperature of gravel -> ground
% CEILING
initial_conditions(8)=7.231039975360592; % Initial avarage temperature of ceiling's first layer
initial_conditions(9)=13.990332663855959; % Initial avarage temperature of ceiling's second layer
initial_conditions(10)=21; % Initial ceiling outside temperature (room -> another room)
% GREY WATER
initial_conditions(11)=14.867435538962756; % Initial grey water temperature
% AMBIENT TEMPERATURE
Tamb0=ydata(end);
% Tamb1, initial
initial_conditions(12)=ydata_x1_avg(1); % avarage temperature in the ground (walls -> gravel)
% Tamb2, initial
initial_conditions(13)=ydata_x2(1); % avarage temperature in the ground (ground -> gravel)

% ROOM'S DIMENSIONS AND PROPERTIES
length=10;
width=10;
height=2;
A_w=2*(length*height)+2*(width*height); % wall area
A_g=length*width; % ground area
A_c=length*width; % ceiling area
V_r=length*width*height; % room volume
dens_air=1.225;
c_air=1006; % specific heat of air, J/kgK
V_air=V_r*dens_air+17*A_g; % '17*A_g' this accounts for termal mass inside the room 17kg/m2 (furniture etc.)
C_r=V_air*c_air;
R_r1=0.13/A_w; % Rsi=0.13 K/W for horizontal heat convection
R_r2=0.17/A_g; % Rsi=0.17 K/W for vertical heat convenction downwards
R_r3=0.17/A_c; % Rsi=0.17 K/W for vertical heat convenction downwards, bc heat is coming from upper room to basement !

% WALLS
d_w(2)=0.2; % concrete layer's thickness 
d_w(3)=0.15; % insulation layer's thickness
d_w(4)=0.5; % gravel's thickness
lambda_w(2)=0.93; % for concrete 0.93 W/(mK) % Source: '2017 Ashrae Handbook, Fundamentals' 33.3
lambda_w(3)=0.036; % for EPS30 0.036 W/(mK) % Source: http://isolite.co.za/more-about-eps/
lambda_w(4)=1.4; % middle value between gravel containing water 2,4 W/(mK) and dry gravel 0,4 W/(mK) % Source: https://www.sciencedirect.com/science/article/pii/S1876610217337931
dens_w(2)=2300; % for concrete 2300 kg/m3 % Source: '2017 Ashrae Handbook, Fundamentals' 33.3
dens_w(3)=30; % for EPS30 30 kg/m3 % Source: http://isolite.co.za/more-about-eps/
dens_w(4)=1920; % for gravel with sand 1920 kg/m3 % source http://www.rfcafe.com/references/general/density-building-materials.htm
c_w(2)=653; % for concrete 653 J/(kgK) % Source: '2017 Ashrae Handbook, Fundamentals' 33.3
c_w(3)=1500; % for EPS30 1500 J/(kgK) % Source: http://isolite.co.za/more-about-eps/
c_w(4)=1850; % for gravel 1850 J/(kgK) % Source: https://www.sciencedirect.com/science/article/pii/S1876610217337931
for i=2:(4)
    R_w(i)=d_w(i)/(lambda_w(i)*A_w); % calculating resistance of walls' layers
    C_w(i)=d_w(i)*A_w*dens_w(i)*c_w(i); % calculating thermals capacitance of walls' layers
end

% GROUND
d_g(5)=0.2; % concrete layer's thickness 
d_g(6)=0.15; % insulation layer's thickness
d_g(7)=0.5; % gravel's thickness
lambda_g(5)=0.93; % for concrete 0.93 W/(mK) % Source: '2017 Ashrae Handbook, Fundamentals' 33.3
lambda_g(6)=0.036; % for EPS30 0.036 W/(mK) % Source: http://isolite.co.za/more-about-eps/
lambda_g(7)=1.4;  % middle value between gravel containing water 2,4 W/(mK) and dry gravel 0,4 W/(mK) % Source: https://www.sciencedirect.com/science/article/pii/S1876610217337931
dens_g(5)=2300; % for concrete 2300 kg/m3 % Source: '2017 Ashrae Handbook, Fundamentals' 33.3
dens_g(6)=30; % for EPS30 30 kg/m3 % Source: http://isolite.co.za/more-about-eps/
dens_g(7)=1920; % for gravel with sand 1920 kg/m3 % source http://www.rfcafe.com/references/general/density-building-materials.htm
c_g(5)=653; % for concrete 653 J/(kgK) % Source: '2017 Ashrae Handbook, Fundamentals' 33.3
c_g(6)=1500; % for EPS30 1500 J/(kgK) % Source: http://isolite.co.za/more-about-eps/
c_g(7)=1850; % for gravel 1850 J/(kgK) % Source: https://www.sciencedirect.com/science/article/pii/S1876610217337931
for i=(5):(7)
    R_g(i)=d_g(i)/(lambda_g(i)*A_g); % calculating resistance of ground's layers 
    C_g(i)=d_g(i)*A_g*dens_g(i)*c_g(i); % calculating thermal capacitance of ground's layers 
end

% CEILING
d_c(8)=0.2; % concrete layer's thickness
d_c(9)=0.15; % insulation layer's thickness
lambda_c(8)=0.93; % for concrete 0.93 W/(mK) % Source: '2017 Ashrae Handbook, Fundamentals' 33.3
lambda_c(9)=0.036; % for EPS30 0.036 W/(mK) % Source: http://isolite.co.za/more-about-eps/
dens_c(8)=2300; % for concrete 2300 kg/m3 % Source: '2017 Ashrae Handbook, Fundamentals' 33.3
dens_c(9)=30; % for EPS30 30 kg/m3 % Source: http://isolite.co.za/more-about-eps/
c_c(8)=653; % for concrete 653 J/(kgK) % Source: '2017 Ashrae Handbook, Fundamentals' 33.3
c_c(9)=1500; % for EPS30 1500 J/(kgK) % Source: http://isolite.co.za/more-about-eps/
for i=(8):(9)
    R_c(i)=d_c(i)/(lambda_c(i)*A_c); % calculating resistance of ceiling's layers 
    C_c(i)=d_c(i)*A_c*dens_c(i)*c_c(i); % calculating thermal capacitance of ceiling's layers 
end
% Room above ceiling
V_ceiling=V_r;
R_c(10)=0.17/A_c; % Room above the besament, because heat is coming from the upper room to basement!
C_ceiling=V_ceiling*dens_air+17*A_c;

% GREY WATER TANK AND FINS
V_gw=m_gw/1000; % tank's volume, m3
r_tank=(V_gw/(3*pi))^(1/3); % tank's radius, m
H_tank=3*r_tank; % tank's height, m
h_air=1.61; % heat convective coefficient for air (calculated using tables) W/(m2K)
k_steel=43; % thermal conductivity for steel in 25C, W/(mK) source - https://www.engineeringtoolbox.com/amp/thermal-conductivity-d_429.html#ampChromeExtension=https%3A%2F%2Fwww.engineeringtoolbox.com%2Fthermal-conductivity-d_429.html
k_pp=0.1; % thermal conductivity for polypropylene (pp) in 25C, W/(mK) source - https://www.engineeringtoolbox.com/amp/thermal-conductivity-d_429.html#ampChromeExtension=https%3A%2F%2Fwww.engineeringtoolbox.com%2Fthermal-conductivity-d_429.html
t_fin=0.002; % fin's thickness, m
L_fin=H_tank; % fin's 'height', m
Lc_fin=L_fin+t_fin/2; % corrected 'height', accounting for uninsolated tip, m
r_fin=0.1; % fin's radius, m
dis=0.05; % distance between fins, m
n_fins=2*pi*r_tank/(dis+t_fin); % number of fins (1 fin per 5 cm)
A_tank=pi*r_tank^2+2*pi*r_tank*H_tank; % tank's surface, m2
A_surf=2*r_fin*Lc_fin+t_fin*Lc_fin; % surface of a single fin, m2
A_root=Lc_fin*t_fin; % surface of a root of a single fin, m2
A_cross=Lc_fin*t_fin; % single fin's cross-section area, m2
P_fin=2*(t_fin+Lc_fin); % single fin's perimeter, m
d_steel=0.003; % thickness of tank's steel layer, m
d_pp=0.003; % thickness of tank's polypropylene layer, m
h_water=245.96; % heat convective coefficient for water (calculated using tables), W/(m2K)
Rconv_water=1/(h_water*A_tank); % convective resistance of water inside the tank, K/W
Rconv_air=1/(h_air*(pi*r_tank^2+(2*pi*r_tank-n_fins*t_fin)*H_tank)); % convective resistance of air surrounding the tank, K/W
Rconv_air1=1/(h_air*(pi*r_tank^2+2*pi*r_tank*H_tank)); % for reference only, convective resistance of air surrounding the tank, IF THERE WERE NO FINS, K/W
Rcond_tank=100000000000000000000000000000000000000000; % conductive restistance of tank's material, K/W
mL=sqrt(h_air*P_fin*r_fin^2/(k_steel*A_cross)); % necessary in order to compute fin's thermal resistance
eff_fin=tanh(mL)/mL; % fin's efficiency, -
E_fin=eff_fin*A_surf/A_cross; % fin effectiveness, -
sigma=5.6703*10^(-8); % The Stefan-Boltzmann Constant, W/(m2K4);
emissivity=0.7; % emissivity for a tank made of steel, % source: https://www.engineeringtoolbox.com/amp/emissivity-coefficients-d_447.html#ampChromeExtension=https%3A%2F%2Fwww.engineeringtoolbox.com%2Femissivity-coefficients-d_447.html
C_gw=m_gw*c_water; % thermal capacitance of grey water
Rrad_tank=1/(4*sigma*emissivity*(pi*r_tank^2+2*pi*r_tank*H_tank)); % for reference only, thermal radiation resistance, IF THERE WERE NO FINS, K/W, Later this expression must be multiplied by '*(((T(1)+T(11))/2)^(-3))'
Rrad_tank_wo_fins=1/(4*sigma*emissivity*(pi*r_tank^2+(2*pi*r_tank-n_fins*t_fin)*H_tank)); % Thermal radiation resistance, TANK'S AREA - FINS' AREA, K/W, Later this expression must be multiplied by '*(((T(1)+T(11))/2)^(-3))'
Rrad_tank_fins=1/(4*sigma*emissivity*A_surf*n_fins); % Thermal radiation resistance, FOR FINS - ALL OF THEM ALTOGETHER, K/W, Later this expression must be multiplied by '*(((T(1)+T(11))/2)^(-3))'
Rcond_fin=(1/(eff_fin*A_surf*h_air))/n_fins; % Thermal resistance of fins (all of them altogether), K/W
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rtotal=Rconv_water+Rcond_tank+1/(1/Rconv_air1+1/(Rrad_tank*((18+273)^(-3))));  % for reference only, total thermal resistance of a tank with no fins at all, K/W
Rtotal_fins=Rconv_water+Rcond_tank+1/(1/Rconv_air+1/Rcond_fin+1/(Rrad_tank_fins*((18+273)^(-3)))+1/(Rrad_tank_wo_fins*((18+273)^(-3))));  % for reference only, total thermal resistance of a tank with fins, K/W
Ratio=Rtotal/Rtotal_fins;  % for reference only, how many times is the resistance of a tank without fins bigger then the resistance of a tank with fins, -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Properties regarding symulation
time_simulation=365*24*3600; % simulation time 
time_charg=8*3600;  % time of heat pump operation
hour=0; % hour counter
day=0; % day counter

% Predefined tables containing results
t_table=zeros(1,time_simulation);
Troom_table=zeros(1,time_simulation);
Tw1_table=zeros(1,time_simulation);
Tw2_table=zeros(1,time_simulation);
Tgravel1_table=zeros(1,time_simulation);
Tg1_table=zeros(1,time_simulation);
Tg2_table=zeros(1,time_simulation);
Tgravel2_table=zeros(1,time_simulation);
Tc1_table=zeros(1,time_simulation);
Tc2_table=zeros(1,time_simulation);
Tceiling_table=zeros(1,time_simulation);
Tgw_table=zeros(1,time_simulation);
Tamb1_table=zeros(1,time_simulation); % temperature in the ground, at wall's depth
Tamb2_table=zeros(1,time_simulation); % temperature in the ground, at grounds's depth
Tamb0_table=zeros(1,time_simulation); % temperature in Duebendorf
Q_hidden_table=zeros(1,time_simulation); % hidden heat (J) transfered from a room above to basement
Q_gw_table=zeros(1,time_simulation); % heat (J) transfered from the grey water tank to basement
dTair=zeros(1,time_simulation); % air's temperature drop in the evaporator
Te=zeros(1,time_simulation); % evaporation temperature, C
Tc=zeros(1,time_simulation); % condensation temperature, C
COP=zeros(1,time_simulation); % COP, -
T_c_in=zeros(1,time_simulation+1); % condenser inlet temperature, C
T_DHW=zeros(1,time_simulation+1); % DHW temperature, C

V_air=510/3600; % avarage air's flow, m3/s, source: http://www.oekotherm.ch/wp-content/uploads/2016/03/MB_Europa-323-DK_de.pdf
m_air=V_air*1.23; % avarage air's mass flow, kg/s
dTair(1)=7; % Troom - Tair_out, we assume the first value, next ones will be calculated, K, source: http://www.oekotherm.ch/wp-content/uploads/2016/03/MB_Europa-323-DK_de.pdf
dtc_k=3; % condenser small temperature difference
dtc=10; % T_c_out - T_c_in, C, we assume constant temperature lift in the condenser as 10 C
dtc_g=dtc+dtc_k; % condenser big temperature difference
dtc_m=dtc/log(dtc_g/dtc_k); % condenser mean temperature difference
dte_k=3; % evaporator small temperature difference
P_h=m_DHW*c_water*(60-10)/time_charg; % Heat flux needed to heat DHW from 10 to 60 C in 'time_charg' time, W, constant
T_DHW(1)=22.5; % the simulation starts at 00:00, that is - during the '22:00 + time_flow every day' period
T_c_in(1)=24.5; % initial condenser inlet temperature, C

    for t=1:1:(time_simulation)
                
        if hour<(6+day*24) || hour>=(22+day*24) % heat pump is running only between 22 pm and 6 am
            Tc(t)=T_c_in(t)+dtc/(1-exp(-dtc/dtc_m)); % condensation temperature is not constant, C
            dte=dTair(t); % Troom - T_air_out = dT_air in evaporator, C
            dte_g=dte_k+dte; % evaporator big temperature difference, C, not constant
            dte_m=dte/log(dte_g/dte_k); % evaporator mean temperature difference, C, not constant
            Te(t)=initial_conditions(1)-dte/(1-exp(-dte/dte_m)); % evaporation temperature, based on Tr, C, not constant
            T_lift=Tc(t)-Te(t); % temp_lift = condensation temp - evaporation temp
            Eff_isentropic=(T_lift^4*p(1)+T_lift^3*p(2)+T_lift^2*p(3)+T_lift^1*p(4)+T_lift^0*p(5))/100; % isentropic efficiency, using a polynomial fitted into data
            COP(t)=Eff_isentropic*(Tc(t)+273)/(Tc(t)-Te(t)); % corrected COP
            P_c=P_h*(COP(t)-1)/COP(t); % Heat flux in evaporator is calculated regarding current COP(t), W
            dTair(t+1)=P_c/(m_air*c_air); % Troom - Tair_out, as before, it depends on COP, C
            P_hp=P_c; % in the set of equation, heat pump's power equals heat flux in evaporator, W
            T_DHW(t+1)=T_DHW(t)+P_h/(m_DHW*c_water); % in a perfectly mixed tank, temperature of DHW is rising linearly
            T_c_in(t+1)=T_DHW(t+1)+2; % we assume that temperature of water entering the condenser is always higher than DHW temperature (at that time) by 2 C
        else % heat pump is turned off
            P_hp=0;
            dTair(t+1)=dTair(1); % Troom - T_air_out, first one assumed for every time heat pump is turned on
            T_DHW(t+1)=10;
            T_c_in(t+1)=10; % at 22:00:01 T_DHW_in=10 C (heat pump is turned on)
        end
        
        if t==((7+day*24)*3600)
           initial_conditions(11)=Thot*0.3+T(end,11)*0.7; % Initial grey water temperature 
        end
        
        if t==((15+day*24)*3600)
           initial_conditions(11)=Thot*0.2+T(end,11)*0.8; % Initial grey water temperature 
        end
        
        if t==((20+day*24)*3600)
           initial_conditions(11)=Thot*0.5+T(end,11)*0.5; % Initial grey water temperature 
        end
        
        if t==(3600*(hour+1)) % hours counter
            hour=hour+1;
            % Temperature on surface, changing with hourly frequencuy
            Tamb0=ydata(hour);
            % Tamb1, changing with hourly frequencuy
            initial_conditions(12)=ydata_x1_avg(hour); % avarage temperature in the ground (walls -> gravel)
            % Tamb2, changing with hourly frequencuy
            initial_conditions(13)=ydata_x2(hour); % avarage temperature in the ground (ground -> gravel)
        end
        
        T_gw_before=initial_conditions(11); % grey water temperature at time 't-1' (necessary for calculation of Qhidden)
        
        tspan=[t-1 t];
        [t,T] = ode45(@(t,T) Basement_model_v12_function(t,T,P_hp,C_r,C_gw,C_w,C_g,C_c,C_ceiling,R_r1,R_r2,R_r3,Rconv_water,Rrad_tank_wo_fins,Rrad_tank_fins,Rcond_tank,Rcond_fin,Rconv_air,R_w,R_g,R_c,n_fins),tspan,initial_conditions);
               
        % ROOM
        initial_conditions(1)=T(end,1); % Initial temperature in the room
        % WALLS
        initial_conditions(2)=T(end,2); % Initial avarage temperature of wall's first layer
        initial_conditions(3)=T(end,3); % Initial avarage temperature of wall's second layer
        initial_conditions(4)=T(end,4); % Initial avarage temperature of gravel -> walls
        % GROUND
        initial_conditions(5)=T(end,5); % Initial avarage temperature of ground's first layer
        initial_conditions(6)=T(end,6); % Initial avarage temperature of ground's second layer
        initial_conditions(7)=T(end,7); % Initial avarage temperature of gravel -> ground
        % CEILING
        initial_conditions(8)=T(end,8); % Initial avarage temperature of ceiling's first layer
        initial_conditions(9)=T(end,9); % Initial avarage temperature of ceiling's second layer
        initial_conditions(10)=21; % Initial ceiling outside temperature (room -> another room)
        % GREY WATER
        initial_conditions(11)=T(end,11); % Initial grey water temperature
               
        time_table(1,t(end))=t(end);
        Troom_table(1,t(end))=T(end,1);
        Tw1_table(1,t(end))=T(end,2);
        Tw2_table(1,t(end))=T(end,3);
        Tgravel1_table(1,t(end))=T(end,4);
        Tg1_table(1,t(end))=T(end,5);
        Tg2_table(1,t(end))=T(end,6);
        Tgravel2_table(1,t(end))=T(end,7);
        Tc1_table(1,t(end))=T(end,8);
        Tc2_table(1,t(end))=T(end,9);
        Tceiling_table(1,t(end))=T(end,10);
        Tgw_table(1,t(end))=T(end,11);
        Tamb1_table(1,t(end))=T(end,12);
        Tamb2_table(1,t(end))=T(end,13);
        Tamb0_table(1,t(end))=Tamb0;
        Q_hidden_table(1,t(end))=C_ceiling*(21-T(end,10));
        Q_gw_table(1,t(end))=C_gw*(T_gw_before-T(end,11));
        
        if hour==(24+day*24) % days counter
        day=day+1;
        end

    end
 
Q_hidden=sum(Q_hidden_table);
Q_gw=sum(Q_gw_table);
P_hidden=Q_hidden/time_simulation; % Heat flux (room above to basement), W
P_gw=Q_gw/time_simulation; % Grey water heat flux, W
COP_avg=sum(COP)/(time_simulation/3); % avarage COP (8h instead of 24h - heat pump operation)

save('Results_model_v13_tc_varrying_Tdhw_without_gw')

start_day=30;

figure(1)
plot(time_table(1,:),Troom_table(1,:),'g',time_table(1,:),Tgw_table(1,:),'r',time_table(1,:),Tamb0_table(1,:),'b',time_table(1,:),Tamb1_table(1,:),'k',time_table(1,:),Tamb2_table(1,:),'c',time_table(1,:),Tceiling_table(1,:),'y');
axis([0 time_simulation -10 Thot+1])
title('Temperatures in time');
xlabel('Time t [s]');
ylabel('Temperatures T [C]');
legend('Troom','Tgw','Tsurface','T1.5','T3.5','Tceiling')
grid on

figure(2)
plot(time_table(1,(start_day*24*3600):((start_day+1)*24*3600)),COP(1,(start_day*24*3600):((start_day+1)*24*3600)),'g');
title('COP(t)');

figure(3)
plot(time_table(1,:),Te(1,1:time_simulation),'g');
title('Tevap(t)');

figure(4)
plot(time_table(1,:),dTair(1,1:time_simulation),'g');
title('dTair(t)');

figure(5)
plot(time_table(1,:),Tc(1,1:time_simulation),'g');
title('Tcond(t)');

figure(6)
plot(time_table(1,:),T_DHW(1,1:time_simulation),'g');
title('TDHWin(t)');

figure(7)
plot(time_table(1,(start_day*24*3600):((start_day+1)*24*3600)),Troom_table(1,(start_day*24*3600):((start_day+1)*24*3600)),'g',time_table(1,(start_day*24*3600):((start_day+1)*24*3600)),Tgw_table(1,(start_day*24*3600):((start_day+1)*24*3600)),'r',time_table(1,(start_day*24*3600):((start_day+1)*24*3600)),Tamb0_table(1,(start_day*24*3600):((start_day+1)*24*3600)),'b',time_table(1,(start_day*24*3600):((start_day+1)*24*3600)),Tamb1_table(1,(start_day*24*3600):((start_day+1)*24*3600)),'k',time_table(1,(start_day*24*3600):((start_day+1)*24*3600)),Tamb2_table(1,(start_day*24*3600):((start_day+1)*24*3600)),'c',time_table(1,(start_day*24*3600):((start_day+1)*24*3600)),Tceiling_table(1,(start_day*24*3600):((start_day+1)*24*3600)),'y');
axis([(start_day*24*3600) ((start_day+1)*24*3600) -10 Thot+1])
title('Temperatures in time during day number 31');
xlabel('Time t [s]');
ylabel('Temperatures T [C]');
legend('Troom','Tgw','Tsurface','T1.5','T3.5','Tceiling')
grid on

figure(8)
plot(time_table(1,(start_day*24*3600):((start_day+1)*24*3600)),COP(1,(start_day*24*3600):((start_day+1)*24*3600)),'g');
title('COP(t)');

figure(9)
plot(time_table(1,(start_day*24*3600):((start_day+1)*24*3600)),Te(1,(start_day*24*3600):((start_day+1)*24*3600)),'g');
title('Tevap(t)');

figure(10)
plot(time_table(1,(start_day*24*3600):((start_day+1)*24*3600)),dTair(1,(start_day*24*3600):((start_day+1)*24*3600)),'g');
title('dTair(t)');

figure(11)
plot(time_table(1,(start_day*24*3600):((start_day+1)*24*3600)),Tc(1,(start_day*24*3600):((start_day+1)*24*3600)),'g');
title('Tcond(t)');

figure(12)
plot(time_table(1,(start_day*24*3600):((start_day+1)*24*3600)),T_DHW(1,(start_day*24*3600):((start_day+1)*24*3600)),'g');
title('TDHWin(t)');

toc;