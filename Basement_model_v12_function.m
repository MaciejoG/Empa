function dTwdt = Basement_model_v12_function(t,T,P_hp,C_r,C_gw,C_w,C_g,C_c,C_ceiling,R_r1,R_r2,R_r3,Rconv_water,Rrad_tank_wo_fins,Rrad_tank_fins,Rcond_tank,Rcond_fin,Rconv_air,R_w,R_g,R_c,n_fins)

dTwdt = [% ROOM
        -P_hp/C_r+(T(11)-T(1))/(C_r*(Rconv_water+Rcond_tank+1/(1/(Rrad_tank_wo_fins*(((T(1)+273+T(11)+273)/2)^(-3)))+1/Rconv_air+1/Rcond_fin+1/(Rrad_tank_fins*(((T(1)+273+T(11)+273)/2)^(-3))))))-(T(1)-T(2))/(C_r*(R_r1+R_w(2)/2))-(T(1)-T(5))/(C_r*(R_r2+R_g(5)/2))-(T(1)-T(8))/(C_r*(R_r3+R_c(8)/2)); % 1
        % WALLS
        1/(C_w(2)*(R_r1+1/2*R_w(2)))*(T(1)-T(2))-1/(C_w(2)*1/2*(R_w(2)+R_w(3)))*(T(2)-T(3)); % 2
        1/(C_w(3)*1/2*(R_w(2)+R_w(3)))*(T(2)-T(3))-1/(C_w(3)*1/2*(R_w(3)+R_w(4)))*(T(3)-T(4)); % 3
        1/(C_w(4)*1/2*(R_w(3)+R_w(4)))*(T(3)-T(4))-1/(C_w(4)*1/2*(R_w(4)))*(T(4)-T(12)); % 4
        % GROUND
        1/(C_g(5)*(R_r2+1/2*R_g(5)))*(T(1)-T(5))-1/(C_g(5)*1/2*(R_g(5)+R_g(6)))*(T(5)-T(6)); % 5
        1/(C_g(6)*1/2*(R_g(5)+R_g(6)))*(T(5)-T(6))-1/(C_g(6)*1/2*(R_g(6)+R_g(7)))*(T(6)-T(7)); % 6
        1/(C_g(7)*1/2*(R_g(6)+R_g(7)))*(T(6)-T(7))-1/(C_g(7)*(1/2*R_g(7)))*(T(7)-T(13)); % 7
        % CEILING
        1/(C_c(8)*(R_r3+1/2*R_c(8)))*(T(1)-T(8))-1/(C_c(8)*1/2*(R_c(8)+R_c(9)))*(T(8)-T(9)); % 8
        1/(C_c(9)*1/2*(R_c(8)+R_c(9)))*(T(8)-T(9))-1/(C_c(9)*(1/2*R_c(9)+R_c(10)))*(T(9)-T(10)); % 9
        % OUTSIDE CEILING
        1/(C_ceiling*(1/2*R_c(9)+R_c(10)))*(T(9)-T(10)); % 10
        % GREY WATER TANK
        -(T(11)-T(1))/(C_gw*(Rconv_water+Rcond_tank+1/(1/(Rrad_tank_wo_fins*((T(1)+273+T(11)+273)/2)^(-3))+1/Rconv_air+1/Rcond_fin+1/(Rrad_tank_fins*((T(1)+273+T(11)+273)/2)^(-3))))); % 11
        % Tamb1
        0; % 12
        % Tamb2
        0]; % 13
    
end

        
        
        
        
        
    
    
