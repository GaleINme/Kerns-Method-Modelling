clear;close

% Goal
fprintf("Before starting please confirm the one-time input parameter has all changed\n")
pause(2)
m_Mate = 107383.4020;   % kg/hr % input prior / one-time input
T_in_Mate = 47.5;         % C % input prior / one-time input
T_out_Mate = 35;        % C % input prior / one-time input
T_in_Cool = 20;         % C % input prior / one-time input
T_out_Cool = 25;        % C % input prior / one-time input

% Step 1 Physical Properties
% Material
Cp_Mate = 2.33;         % kJ/kgC % input prior / one-time inpu
rho_Mate = 648;         % kg/m3 % input prior / one-time inpu
vis_Mate = 0.258*10^-3;  % Ns/m2 % input prior / one-time inpu
k_Mate = 0.113;          % W/mC % input prior / one-time input
fouling_coeff_Mate = 2000; % input prior / one-time input
%Coolant
Cp_Cool = 4.2;          % kJ/kgC % input prior / one-time input
rho_Cool = 997.65;         % kg/m3 % input prior / one-time input
vis_Cool = 0.952*10^-3;   % Ns/m2 % input prior / one-time input
k_Cool = 0.6;          % W/mC % input prior / one-time input
fouling_coeff_Cool = 6000; % input prior / one-time input
confirmation = input("If the Parameter has not yet changed, please enter 'n'","s");
if confirmation == 'n'
    return
end


% Step 2 Duty 
Q = (m_Mate/3600)*Cp_Mate*(T_in_Mate-T_out_Mate);
m_Cool = Q/(Cp_Cool*(T_out_Cool-T_in_Cool));

% Step 3 Assume U
fprintf("\nThe heat duty is %.2f kW and the mass flowrate for coolant is %.2f kg/s\n",Q,m_Cool)
U_guess = input("Input initial guess for U");


% Step 4 Calculate delt_Tlm, F and delt_Tm
% Log Mean Temperature 
delta_T1 = T_in_Mate-T_out_Cool;
delta_T2 = T_out_Mate-T_in_Cool;
delta_Tlm = (delta_T1-delta_T2)/log(delta_T1/delta_T2);
delta_T_Mate = T_in_Mate-T_out_Mate;
delta_T_Cool = T_out_Cool - T_in_Cool;
R = delta_T_Mate/delta_T_Cool;
S = delta_T_Cool/(T_in_Mate - T_in_Cool);
% F_T Single Shell Even Tube Pass (Correction Factor)
F_T_Num = sqrt(R^2+1)*log((1-S)/(1-R*S));
F_T_Den_log = log(((2-S*(R+1-sqrt(R^2+1)))/((2-S*(R+1+sqrt(R^2+1))))));
F_T = F_T_Num/((R-1)*F_T_Den_log);

delta_Tm = F_T*delta_Tlm;
if delta_Tm < 0.75
    fprintf("The delta mean temperature is lower than 0.75, please revise")
    return
end

% Step 5 Heat Transfer Area
H_A_T = Q*1000/(U_guess*delta_Tm);

% Step 6 Choose Material and Pipe Sizing
fprintf("Total heat transfer area = %.2f\n",H_A_T)
k_tube = input("WHats the thermal conductivity of tube?"); % Thermal conductivity of tube wall
OD_Tube = input("WHats the chosen outer diameter of tube? Unit in mm")*10^-3;
thickness_Tube = input("Whats the chosen thickness of tube? Unit in mm")*10^-3;
ID_Tube = OD_Tube - 2*thickness_Tube; %m
length_Tube = input("Whats the chosen length of tube?, m");  

% Step 7 Number of Tube
area_Single_Tube = pi*OD_Tube*length_Tube;
number_Tube = H_A_T/area_Single_Tube;
if mod(ceil(number_Tube),2)==1
    number_Tube = ceil(number_Tube) +1;
else
    number_Tube = ceil(number_Tube);
end

% Step 8 Shell diameter
while true
    fouling_pitch = input("Is the fluid clean, y/n?","s");
    if fouling_pitch == 'y'
        number_of_pass = 2;
        if number_of_pass == 2
            k1 = 0.319;
            n1 = 2.142;
        elseif number_of_pass == 4
            k1 = 0.175;
            n1 = 2.285;
        elseif number_of_pass == 6
            k1 = 0.0743;
            n1 = 2.499;
        elseif number_of_pass == 8
            k1 = 0.0365;
            n1 = 2.675;
        else
            fprinf("No such number, retry!")
        end
        break
    else
        number_of_pass = input("How many passes, '2,4,6,8'?");
        if number_of_pass == 2
            k1 = 0.156;
            n1 = 2.291;
        elseif number_of_pass == 4
            k1 = 0.158;
            n1 = 2.263;
        elseif number_of_pass == 6
            k1 = 0.0402;
            n1 = 2.617;
        elseif number_of_pass == 8
            k1 = 0.0331;
            n1 = 2.643;
        else
            fprinf("No such number, retry!")
        end
        break
    end
end

bundle_Diameter = OD_Tube*1000*(number_Tube/k1)^(1/n1); % mm
fprintf("The bundle diameter is %.2f\n",bundle_Diameter)
clearance = input("Whats the clearance got feom the graph for shell? Unit in mm"); % input needed
shell_Diameter = bundle_Diameter + clearance;
tube_Pitch = 1.25*OD_Tube;

% Step 9 Tube Side Film Coefficient
cross_A_Tube = pi*(ID_Tube^2)/4;
v_Cool = m_Cool/(rho_Cool*(cross_A_Tube*number_Tube/2));
Re_Tube = rho_Cool*v_Cool*ID_Tube/vis_Cool;
fprintf('The reynolds number of tube is %d\n',Re_Tube)
jH_Tube = input("Whats the value of heat transfer factor of tube? The default unit is *10^-3")*10^-3; % input from graph
Pr_Tube = Cp_Cool*10^3*vis_Cool/k_Cool;
Nu_Tube = jH_Tube*Re_Tube*Pr_Tube^0.33;
tube_Film_Coefficient = Nu_Tube*k_Cool/ID_Tube;

% Step 10 Sheel Side Film Coefficient
baffle_coefficient = input("Whats the baffle coefficient for baffle spacing? (usual 0.2)");
baffle_Spacing = baffle_coefficient*shell_Diameter*10^-3;
cross_A_Shell = (tube_Pitch-OD_Tube)*shell_Diameter*10^-3*baffle_Spacing/tube_Pitch;
v_Mate = (m_Mate/3600)/(rho_Mate*cross_A_Shell);
if fouling_pitch == 'y'
    shell_Equi_Diameter = (1.1/OD_Tube)*(tube_Pitch^2-0.917*OD_Tube^2);
else
    shell_Equi_Diameter = (1.27/OD_Tube)*(tube_Pitch^2-0.785*OD_Tube^2);
end
Re_Shell = rho_Mate*v_Mate*shell_Equi_Diameter/vis_Mate;
fprintf('The reynolds number of shell is %d\n',Re_Shell)
jH_Shell = input("Whats the value of heat transfer factor of shell? The default unit is *10^-3")*10^-3; % input from graph
Pr_Shell = Cp_Mate*10^3*vis_Mate/k_Mate;
Nu_Shell = jH_Shell*Re_Shell*Pr_Shell^0.33;
shell_Film_Coefficient = Nu_Shell*k_Mate/shell_Equi_Diameter;

% Step 11 Overall U for Heat Exchanger
U_inverse_1 = OD_Tube/(tube_Film_Coefficient*ID_Tube);
U_inverse_2 = OD_Tube/(fouling_coeff_Cool*ID_Tube);
U_inverse_3 = (OD_Tube*log(OD_Tube/ID_Tube)/(2*k_tube));
U_inverse_4 = (1/fouling_coeff_Mate)+(1/shell_Film_Coefficient);
U_calc = 1/(U_inverse_1+U_inverse_2+U_inverse_3+U_inverse_4);
U_different = (U_calc-U_guess)/U_guess*100;

% Step 12 Calculate Pressure Drop
fprintf('The reynolds number of tube is %d\n',Re_Tube)
tube_fric_factor =input("Whats the value of friction factor of tube? The default unit is *10^-3")*10^-3; % input from fric-Re graph
tube_Pressure_Drop = number_of_pass*(8*tube_fric_factor*(length_Tube/ID_Tube)+2.5)*((rho_Cool*v_Cool^2)/2);
fprintf('The reynolds number of shell is %d\n',Re_Shell)
shell_fric_factor = input("Whats the value of friction factor of shell? Please enter in 0.0xx"); % input from fric-Re graph
shell_Pressure_Drop = 8*shell_fric_factor*((shell_Diameter*10^-3)/shell_Equi_Diameter)*(length_Tube/baffle_Spacing)*((rho_Mate*v_Mate^2)/2);


% Step 13 Cost Estimation
% Print result
if U_different > 30 || U_different < 0
    fprintf("The difference is not within the setting range, it's %.2f %\n",U_different)
else
    if tube_Pressure_Drop > 35e+3
        fprintf("The tube pressure drop is too high, it's %.2f\n",tube_Pressure_Drop)
    else
        if shell_Pressure_Drop > 100e+3
            fprintf("The shell pressure drop is too high, it's %.2f\n",shell_Pressure_Drop)
        else
            fprintf("The result is valid")
        end
    end
end