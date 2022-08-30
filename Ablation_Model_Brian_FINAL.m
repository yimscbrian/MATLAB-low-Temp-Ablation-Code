%% README 

%MATLAB Low Temperature Ablation Code by Brian YIM
% - Adapted from Lorna Barron and Laurance Maskell

% INSTRUCTIONS
% 1. Set run parameters (marked "-MODIFIABLE")
% 2. Make changes to model material if required
% 3. Run code

%% HOUSEKEEPING
clear all
close all

%% INTIALISATION
% Solution Parameters
eps = 1e-6;                     % (-) Solution Tolerance for Gauss-Seidel Solver
maxiter = 20000;                % (-) Maximum number of iterations for Gauss-Seidel Solver
valueholder=1;

% Spatial Discretisation
npi = 100;          % (-) Cells in the radial direction
npj = 182;            % (-) Cells in the tangential direction. The current set-up allows for equal sizes (1deg)


% Time Discretisation
Tend = 60;        % (s) End time                            - MODIFIABLE
dt = 0.1;            % (s) Time step                       - MODIFIABLE
nt = round(Tend/dt);       % (-) Number of time steps

% Freestream Conditions (Wind Tunnel)
P_0 = 50;        % (kPa) Stagnation Pressure of the flow - MODIFIABLE   
T_0 = 300;%500;         % (K) Stagnation temperature of the flow - MODIFIABLE
Ma = 4;%6;         % (-) Mach number of the flow               - MODIFIABLE

% Heatshield Geometry
rin = 0.015;    % (m) Initial radius of model               - MODIFIABLE
rmin = 0.05 * rin;    % (m) Inner radius of heat shield (not a physical feature, but to provide a limit for the solver)
shape = "Standard";   % 'Standard' or 'Elliptical' geometry - MODIFIABLE

% Material Selection
material.mat = "Dry Ice";       % Model material            - MODIFIABLE
% Select 'Dry Ice', 'Naphthalene', or 'Water Ice'
material=getMaterialProperties(material);

% Properties of Air
gamma = 1.4;                % (-) Ratio of specific heats of air
R_air = 287;                    % (J/kgK) Gas constant for air
cp_air = 1004;                  % (J/kgK) Specific Heat Capicity for Air assuming constant pressure
nu_air = 1.568e-05;             % (m2/s) Kinematic viscosity of air
k_air = 0.025;                  % (W/mK) Thermal conductivity of air
k_sg = 1.7415e-04;              % (kg0.5/m) Sutton-Graves constant
Pr = 0.7;                       % Prandtl Number

% Experimental data for pressure correction
useNewtonian = 0;
%caldata=readmatrix('calwithMach.csv'); %Removed as data hardcoded in
%caldata correction factor data 
%Data Layout:
% [8-deg cone Mach,8-deg cone Factor, Hemisphere Mach, Hemisphere Factor, 9-deg cone Mach, 9-deg cone Factor;] 
%First 2 rows left in to match .csv layout
caldata=[NaN,NaN,NaN,NaN,NaN,NaN;...
    NaN,NaN,NaN,NaN,NaN,NaN;...
    0.410934359000000,-0.838509317000000,0.596278433000000,-0.866459627000000,1.49160109900000,-0.124223602000000;...
    0.604556821000000,-0.776397516000000,0.703566883000000,-0.788819876000000,1.88559796900000,-0.0807453420000000;...
    0.901757619000000,-0.534161491000000,0.791049377000000,-0.720496894000000,2.28708042700000,-0.0931677020000000;...
    0.956178385000000,-0.465838509000000,0.892480637000000,-0.596273292000000,2.95161616600000,-0.105590062000000;...
    1.06401097500000,-0.304347826000000,1,-0.459627329000000,3.97010867700000,-0.118012422000000;...
    1.17182369500000,-0.192546584000000,1.09376030000000,-0.285714286000000,4.65226741900000,-0.136645963000000;...
    1.53330567300000,-0.0683229810000000,1.19219502600000,-0.201863354000000,NaN,NaN;...
    2.00629693500000,-0.0341614910000000,1.38267007300000,-0.0683229810000000,NaN,NaN;...
    2.50151777300000,-0.0372670810000000,1.73588503100000,-0.0186335400000000,NaN,NaN;...
    3.50679339300000,-0.0559006210000000,1.89210881000000,-0.00310559000000000,NaN,NaN;...
    4.49463687800000,-0.0621118010000000,2.07665938700000,-0.0124223600000000,NaN,NaN;...
    NaN,NaN,2.35102628700000,-0.0248447200000000,NaN,NaN;...
    NaN,NaN,2.99259471200000,-0.0465838510000000,NaN,NaN;...
    NaN,NaN,3.99757303700000,-0.0465838510000000,NaN,NaN;...
    NaN,NaN,4.98430707200000,-0.0434782610000000,NaN,NaN];
caldata=caldata(3:end,:); %resize data
caldata(isnan(caldata))=0;
calF=griddedInterpolant(caldata(1:length(find(~all(caldata(:,1)==0,2))),1),caldata(1:length(find(~all(caldata(:,1)==0,2))),2));
%% INITIALISE
% Generate initial solution mesh
[r, th, rmax, c, dr, dth] = create_grid(npi, npj, rin, rmin, shape);    % Creates the solution space of the model
[X,Y] = grid_plotting(r, th, rmax, npi, npj, rin, 1);                     % Plots the model and generated grid
theta_cell=(-pi()/2:pi()/(npj-2):pi()/2)';
dth_deg=dth*180/pi();

% Preallocation of various arrays
normqe = zeros(nt,1);               % (-) [nt x 1]          Normalised Conductive Heat
normqrad = zeros(nt,1);             % (-) [nt x 1]          Normalised Radiative Heat
normqsub = zeros(nt,1);             % (-) [nt x 1]          Normalised Sublimative Heat
thickness=zeros(nt,1);
qconv_history = zeros(nt,1);        % (W/m2K) [nt x 1]          Normalised Conductive Heat
qrad_history = zeros(nt,1);   
qcond_history = zeros(nt,1);   
qsub_history = zeros(nt,1);   
qconv_mean_history = zeros(nt,1);        % (W/m2K) [nt x 1]          Normalised Conductive Heat                
qrad_mean_history = zeros(nt,1);        % (W/m2K) [nt x 1]          Normalised Conductive Heat
qcond_mean_history = zeros(nt,1);        % (W/m2K) [nt x 1]          Normalised Conductive Heat
qsub_mean_history = zeros(nt,1);        % (W/m2K) [nt x 1]          Normalised Conductive Heat
delta_r = zeros(nt,npj);             % (m) [nt x npj]        Recession at each timestep
histr = zeros(nt, npi, npj);        % (m) [nt x npi x npj]  History of radial grid
histrdot = zeros (nt,1);            % (m/s) [nt x 1]        History of recession rate
histk = zeros (nt,1);
histrmax = zeros(nt, npj);          % (m) [nt x npj]        History of outer surface dimension
histT = zeros(nt, npi, npj);        % (K) [nt x npi x npj]  History of temperature field
histc = zeros(nt, npj);             % (m-1) [nt x npj]      History of the curvature
mdotstag = zeros(nt, npj);          % (kg/m2s) [nt x npj]        History of surface mass flow
rstag_hist=zeros(nt,1);

T_inf=T_0/(1+(gamma-1)/2*Ma^2);  % (K) Freestream temperature
%T_inf = T_0/(1+(Ma^2*gamma*R_air)/(2*cp_air));
V_inf = Ma*((gamma*R_air*T_inf)^(1/2));         % (m/s) Freestream velocity
p_inf = P_0*1000/(1+(gamma-1)/2*Ma^2)^(gamma/(gamma-1));              % (Pa) Freestream pressure
rho_air = p_inf/(R_air * T_inf);        
rho_inf = rho_air;
rf = Pr^(1/3);                                      % (-) Recovery factor of turbulent flow

% Initialise the temperature field and boundary conditions 
%T_sub=material.T_sub;
T_in=material.T_in;
T = ones(npi, npj) .* T_in;

%% TIME ITERATION
% Run heat transfer from t=0 to t=t_end
for timeIter = 1:nt
 %Obtain surface geometry and paramters at timeIter
 [dS,dydx,theta_ramp,k,k_n]=getSurface(rmax,th,theta_cell,npi,npj);
 %Update Cp and Pressure from surface geometry
 [Cp, Ma_e, p_e, T_e] = modified_newtonian_theory(theta_ramp, Ma, gamma, T_inf, p_inf, T_0, rho_air, R_air, useNewtonian,calF, material,npj);
 %Calculate heat transfer (conductive and sublimative heat flux)
 [ q_cond,q_conv,q_rad,q_sub, T_w,delta_r,mdot,rdot,T_sub] = ...
    heatBalance(rmax, R_air, cp_air, gamma, nu_air, k_air, V_inf, rho_inf, p_inf, Pr, material, npi,npj, p_e, T_e, Ma_e, dS, T, dt, dth, th,k,r,k_n);
               
 %Solve heat equation for heat flux and temperature boundary conditions
      [aN, aE, aS, aW, aP, bP, sP] = Tcoeff(dth, dr, r, c, material, dt, T, rmax, q_cond, npi, npj,T_sub);

      % Gauss-Siedel Solver
      for iter = 1:maxiter
          if(iter == 1)
              % In the first iteration, there is no previously calculated temperature field, and so the current temperature field is provided before it is updated
              Told = T;
          end
          if mod(iter,2) == 1 % check to see if odd or even iteration
              %odd iteration
              [T] = gs_solver(T, aN, aE, aS, aW, aP, bP, npi, npj,T_sub);     % Solving for the temperature field using a Gauss-Seidel Solver

          else
              %even iteration
              [T] = gs_solver_even(T, aN, aE, aS, aW, aP, bP, npi, npj,T_sub); 
          end

          if(iter == 1) 
              norm0 = norm2matrix(Told, T, dr, dth);              % Computing the norm of the two temperature fields
              Told = T;                                           % Overwriting the previous temperature field
          else
              normNext = norm2matrix(Told, T, dr, dth);           % Computing the norm of the two temperature fields
              if(normNext < (eps*norm0))
                  % If the calculated norm is smaller than the applied tolereance, solution has converged
                  break
              else
                  % If the calculated norm is greater than the applied tolerance, calculate next iteration
                  Told = T;
              end
          end
      end


qconv_history(timeIter) = q_conv(ceil(npj/2));                % (W/m2) Storing convetive heat flux
qrad_history(timeIter) = q_rad(ceil(npj/2));
qcond_history(timeIter) = q_cond(ceil(npj/2));
qsub_history(timeIter) = q_sub(ceil(npj/2));

qconv_mean_history(timeIter) = mean(q_conv);                % (W/m2) Storing convetive heat flux
qrad_mean_history(timeIter) = mean(q_rad);
qcond_mean_history(timeIter) = mean(q_cond);
qsub_mean_history(timeIter) = mean(q_sub);

rstag=rmax(round(length(rmax)/2),1);
rstag_hist=rstag;

normqe(timeIter) =  mean(q_cond)/mean(q_conv);          % (-) Normalised conductive heat transfer
normqrad(timeIter) = mean(q_rad)/mean(q_conv);          % (-) Normalised radiative heat transfer
normqsub(timeIter) = mean(q_sub)/mean(q_conv);          % (-) Normalised sublimative heat transfer
thickness(timeIter)= rstag;

% STEP 6. SAVE SURFACE RECESSION AMOUNT TO HISTORY MATRIX
histrdot(timeIter)=min(rdot);
histr(timeIter, :, :) = r;
histrmax(timeIter, :) = rmax;
histT(timeIter, :, :) = T;
histc(timeIter, :) = c;
histk(timeIter)=k_n;
mdotstag(timeIter, 1) = mdot(round(npj/2));%./dt;        % Surface recession at nose
mdotstag(timeIter, 2) = mdot(round(npj/3));%./dt;        % Surface recession at theta = 30 degrees
mdotstag(timeIter, 3) = mdot(round(npj/6));%./dt;        % Surface recession at theta = 60 degrees

% STEP 5. REMESH GRID USING CURVATURE OF THE CICRLCE
[c, dr, r, rmax] = remesh(r, th, rmax, rmin, c, dr, delta_r, npi, npj);

if mod(timeIter,100)==0 %plot curve every 100 timesteps
    figure(8)
    plot(th.*180./pi, q_conv, 'LineWidth', 2)
    M(valueholder) = strcat("t=", num2str(timeIter*dt),"s");
    hold on
    valueholder=valueholder+1;
    grid on
ylabel('Convective Heat Flux $q_{conv}$ [W/m$^2$]', 'Interpreter', 'latex')
xlabel('Angle $\phi$ [$^{\circ}$]', 'Interpreter', 'latex')
legend(M);
set(gca, 'FontSize', 12,'TickLabelInterpreter','latex')
%Plot inital convective heat flux curve (not used)
elseif timeIter==1
    figure(8)
    plot(th.*180./pi, q_conv, 'LineWidth', 2)
    M(valueholder) = strcat("t=", num2str(0),"s");
    hold on
    valueholder=valueholder+1;
    grid on
end

% CHECK FOR MIMIMUM r-VALUE
if(min(rmax)< 0.2 * rin)
    % Checks if the heat shield has completely sublimated away. If so, program exits
    disp('Reached mimimum radius, program aborted');
    break
elseif(c(round(npj/2)) < 0.001*0.25/rin)
    % Checks if the heat shield is no longer circular enough. If so, program exits
    disp('Reached mimimum allowable curvature, program aborted');
    break
end
fprintf('Number of iterations: %i \t Time: %f\n',iter,timeIter*dt);
end


%% PLOTTING OF RESULTS

[Xend,Yend] = grid_plotting(r, th, rmax, npi, npj, rin, 3); % Plotting of model at final timestep
time1 = zeros(1,nt);
time = [1:timeIter] * dt;                                   % Time simulated successfully by model (present in case the model exits early due to a failure condition)
time1(1:length(time))=time;

%Put data into array
data=[time1' qcond_history qconv_history qsub_history qrad_history qcond_mean_history qconv_mean_history qsub_mean_history qrad_mean_history histrdot histrmax(:,round(npj/2)) mdotstag(:,1:3)];

%Create a name string for naming saved figures
nameMatrix =['Ablation' material.mat ' P_0=' P_0 ' T_0=' T_0 ' M=' Ma]; %Video file name


figure(4)   % FIGURE 4 - COMPARING PROPORTIONS OF ENTRY HEATING (CONDUCTION VS RADIATION VS SUBLIMATION)
plot(time,normqe(1:timeIter),  '-r', 'LineWidth', 2)
hold on
grid on
plot(time, normqrad(1:timeIter), '-k', 'LineWidth', 2)
plot(time, normqsub(1:timeIter), '-b', 'LineWidth', 2)
legend('q_{cond}', 'q_{rad}', 'q_{sub}', 'Interpreter', 'tex', 'Location', 'East')
ylabel('$\stackrel{.}{q}$/$\stackrel{.}{q}_{conv}$', 'Interpreter', 'latex')
xlabel('t [s]', 'Interpreter', 'latex')
set(gca, 'FontSize', 12,'TickLabelInterpreter','latex')
exportgraphics(gca,char(strjoin([nameMatrix ' Normalized Heat Flux' '.png'])),"Resolution",300)

figure(5)   % FIGURE 5 - SURFACE MASS FLOW RATES AT VARIOUS ANGLES
plot(time, mdotstag(1:timeIter, 1), 'LineWidth', 2)
hold on
grid on
plot(time, mdotstag(1:timeIter, 2), 'LineWidth', 2)
plot(time, mdotstag(1:timeIter, 3), 'LineWidth', 2)
ylabel('$\stackrel{.}{m}$ [kg/$m^2$s]', 'Interpreter', 'latex')
xlabel('t [s]', 'Interpreter', 'latex')
legend('\phi = 0\circ', '\phi = 30\circ', '\phi = 60\circ', 'Interpreter', 'tex','Location', 'NorthWest')
set(gca, 'FontSize', 12,'TickLabelInterpreter','latex')
exportgraphics(gca,char(strjoin([nameMatrix ' Mass flux' '.png'])),"Resolution",300)

figure(6)   % FIGURE 6 - COMPARING PROPORTIONS OF ENTRY HEATING (CONDUCTION VS RADIATION VS SUBLIMATION)
plot(time, thickness(1:timeIter), '-k', 'LineWidth', 2)
hold on
grid on
%legend('Thickness', 'Interpreter', 'tex', 'Location', 'East','FontSize', 24)
ylabel('Nose Radius $R_N$ [m]', 'Interpreter', 'latex') %'$\stackrel{.}{q}$/$\stackrel{.}{q}_{conv}$'
xlabel('t [s]', 'Interpreter', 'latex')
set(gca, 'FontSize', 12,'TickLabelInterpreter','latex')
exportgraphics(gca,char(strjoin([nameMatrix ' Nose Radius' '.png'])),"Resolution",300)

figure(7)   % FIGURE 7 - SURFACE RECESSION AT NOSE
plot(time, histrdot(1:timeIter, 1)*1000, 'LineWidth', 2)
hold on
grid on
ylabel('Recession Rate $\stackrel{.}{r}$ [mm/s]', 'Interpreter', 'latex')
xlabel('t [s]', 'Interpreter', 'latex')
set(gca, 'FontSize', 12,'TickLabelInterpreter','latex')
exportgraphics(gca,char(strjoin([nameMatrix ' Recession Rate' '.png'])),"Resolution",300)

figure(8)   % FIGURE 8 - CONVECTIVE HEAT FLUX AT TIMESTEP
if mod(timeIter,100)~=0
plot(th.*180./pi, q_conv, 'LineWidth', 2)
M(valueholder) = strcat("t=", num2str(timeIter*dt),"s");
hold on
end
legend off
exportgraphics(gca,char(strjoin([nameMatrix ' Convective Heat flux over time' '.png'])),"Resolution",300)
legend(M)

recessionrate_average=mean(histrdot(timeIter-length(thickness(thickness(1:timeIter)~=rin)):timeIter));
fprintf('Average Recession Rate: %fmm/s \n',recessionrate_average*1000);
fprintf('Maximum Recession Rate: %fmm/s',min(histrdot)*1000);

%Plots Temperature field at all times and saves as video
plotResults(histT,npi,npj,th,histr,timeIter,dt,T,T_in,rin,material,P_0,T_0,Ma, histrdot)


%% FUNCTIONS
function material=getMaterialProperties(material)
% Selects model material
% Properties of the model material
if(strcmp(material.mat, "Water Ice"))
    material.R_abl = 461.52;                             % (J/kgK) Gas constant for water
    material.cp_abl = 4182;                              % (J/kgK) Specific Heat Capicity for water assuming constant pressure
    material.rho_abl = 917;                              % (kg/m3) Density of Ice in solid form
    material.nu_abl = 1e-06;                             % (m2/s) Kinematic viscosity of water
    material.k_abl = 0.598;                              % (W/mK) Thermal conductivity of water (298K)
    material.h_abl = 2.8386e06;                          % (J/kg) Enthalpy of Sublimation of water
    material.T_in=255;%77.36;                            % (K)  Initial Temeprature of model - Typical Freezer Temperature
    material.T_sub=273;
    A=-7.342973E+0;
    B=-7.276391E+03;
    C=6.702455E+01;
    D=4.161914E-06;
    T=50:1:400;
    material.W_abl=18.01528; %(g/mol)
    material.Diffusion_coef=2.36*10^-5;                  %(m^2/s) Diffusion coeffiecient of ablator into air (STP)

elseif(strcmp(material.mat, "Naphthalene"))
    material.R_abl = 64.87;                              % (J/kgK) Gas constant for naphthalene
    material.cp_abl = 1300;                              % (J/kgK) Specific Heat Capicity for napthalene assuming constant pressure
    material.rho_abl = 1140;                             % (kg/m3) Density of naphthalene in solid form
    material.nu_abl = 9.689e-07;                         % (m2/s) Kinematic viscosity of naphthalene
    material.k_abl = 0.1221;                             % (W/mK) Thermal conductivity of naphthalene
    material.h_abl = 5.57e+05;                           % (J/kg) Enthalpy of Sublimation of naphthalene
    material.T_in=298;
    material.T_sub=353.3;
    A=-7.670601E+00;
    B=-8.563701E+03;
    C=6.888815E+01;
    D=2.862167E-06;
    T=250:1:800;
    material.W_abl=128.1705; %(g/mol)
    material.Diffusion_coef=8.6e-6;         %(m^2/s) Diffusion coeffiecient of ablator into air (data assumed 303.2K Temp 1atm Pressure)

elseif(strcmp(material.mat, "Dry Ice"))
    material.R_abl = 188.9;                              % (J/kgK) Gas constant for CO2
    material.cp_abl = 849;                               % (J/kgK) Specific Heat Capicity for CO2 assuming constant pressure
    material.rho_abl = 1001;%1600;                       % (kg/m3) Density of Dry Ice in solid form (1001 is density from Anderson)
    material.nu_abl = 8.3399e-07;                        % (m2/s) Kinematic viscosity of CO2
    material.k_abl = 0.1422;                             % (W/mK) Thermal conductivity of liquid CO2 (250K) [0.0109 W/mK @220K]
    material.h_abl = 5.58e05;                            % (J/kg) Enthalpy of Sublimation of Dry Ice (@180K)
    material.T_in=77.36;                                 % (K) Boiling point of liquid Nitrogen
    material.T_sub=216.5;
    A=-2.403761E+01;
    B=-7.062404E+03;
    C=1.663861E+02;
    D=3.368548E-05;
    T=50:1:400;
    material.W_abl= 44.01; %(g/mol)
    material.Diffusion_coef=1.38*10^-5; %(m^2/s) Diffusion coeffiecient of ablator into air (STP)

elseif(strcmp(material.mat, "Camphor"))
    material.R_abl = 54.61;                              % (J/kgK) Gas constant for CO2
    material.cp_abl = 1871.16;                           % (J/kgK) Specific Heat Capicity for camphor assuming constant pressure
    material.rho_abl = 953;                              % (kg/m3) Density of Camphor in solid form (Anderson)
    material.k_abl = 0.002015;                           % (W/mK) Thermal Conductivity of Camphor (Anderson)
    material.h_abl = 3.52e05;                            % (J/kg) Enthalpy of Sublimation of Camphor (Anderson)
    material.T_in=298;                                   % (K) Room Temperature
    material.T_sub=453.1;
    A=-2.026514E+01;
    B=-1.292476E+04;
    C=1.544389E+02;
    D=9.581395E-06;
    T=250:1:800;
    material.W_abl= 153.23; %(g/mol)
    material.Diffusion_coef=0.0009; %(m^2/s) Diffusion coeffiecient of ablator into air (STP)

end
material.Pvps = exp(A.*log(T) +B./T + C  + D*T.^2)*1000;
material.Pvp_curve=griddedInterpolant(T,material.Pvps);            % (Pa) Vapour Pressure curve
material.Ts=T;
material.alpha_abl = material.k_abl / (material.rho_abl * material.cp_abl);     % (m2/s) Thermal Diffusivity of material in solid form

end

function [X,Y] = grid_plotting(r, th, rmax, npi, npj, rin, fig)
% This function plots the current model shape
%
% INPUTS:
%   r       :   (m)     [npi x npj]     Coordinate field of the model   
%   th      :   (rad)   [npj x 1]       Angular spacing of the model
%   rmax    :   (m)     [npj x 1]       Surface radius
%   npi     :   (-)                     Number of spaces in the radial direction
%   npj     :   (-)                     Number of spaces in the angular or tangential direction 
%   fig     :   (-)                     Figure number parameter
% OUTPUTS:
%   X       :   (m)     [npi x npj]     X-coordinates of the model
%   Y       :   (m)     [npi x npj]     Y-coordinates of the model

    % Cartesian Coordinate conversion
    Xmax = rmax .* cos(th);
    Ymax = rmax .* sin(th);
    X = zeros(npi, npj);
    Y = zeros(npi, npj);    
    for i = 1:npi
        X(i, :) = r(i,:)' .* cos(th);
        Y(i, :) = r(i,:)' .* sin(th);
    end

    figure(fig)
    plot(Xmax, Ymax, 'k',"MarkerSize", 12);
    hold on
    grid on
    title('Solution Grid');
    xlabel('Horizontal Axis (m)');
    ylabel('Vertical Axis (m)');
    axis equal
    ylim([-rin rin])
    xlim([0 rin*2])
    for j = 1:npi
        plot(X(j, :), Y(j, :), 'b--');
    end
    for k = 2:npj-1
        x1 = rmax(k) .* cos(th(k));
        y1 = rmax(k) .* sin(th(k));
        line([0 x1], [0 y1], 'Color', 'red', 'LineStyle', '--');
    end
    line([0 0], [0 rmax(npj)], 'Color', 'red', 'LineStyle', '--');
    line([0 0], [0 -rmax(npj)], 'Color', 'red', 'LineStyle', '--');
end

function [Cp, Ma_e, p_e, T_e] = modified_newtonian_theory(theta_ramp, Ma, gamma, T_inf, p_inf, T_0, rho_air, R_air, useNewtonian,calF, material, npj)

% This function calculates flow parameters using modified newtonian theory
%
% INPUTS:
%   th          :   (rad)   [npj x 1]       Angular spacing of the model
%   Ma          :   (-)                     Mach number of the freestream flow
%   gamma       :   (-)                     Ratio of specific heats for air
%   T_inf       :   (K)                     Temperature of the freestream flow
%   V_inf       :   (m/s)                   Velocity of the freestream flow
%   rho_air     :   (kg/m3)                 Density of the freestream flow
%   R_air       :   (J/kgK)                 Specific gas constant for air
%   cp_air      :   (J/kgk)                 Specific heat capacity of air at constant pressure
%   material    :   ()                      model material
%   useNewtonian:   ()                      Determines whether modified
%   Newtonian Theory with correction factor is used
%   calF        :   ()                      Correction factor Interpolant

% OUTPUTS:
%   Cp      :   (-)     Local pressure coefficient at the surface
%   Ma_e   :   (-)     Mach number at the edge of the boundary layer
%   p_e    :   (Pa)    Pressure at the edge of the boundary layer
%   T_e    :   (K)     Temperature at the edge of the boundary layer


    Cp_max =2/(gamma*Ma^2)*((((gamma+1)*Ma)^2/(4*gamma*Ma^2-2*(gamma-1)))^(gamma/(gamma-1))*((1-gamma+2*gamma*Ma^2)/(gamma+1))-1); %calculate maximum Cp using Rayleigh Pitot formula
    
    if useNewtonian==1 %If ==1
        Cp = Cp_max .* (sin(theta_ramp).^2);
    else
        %HYPERSONIC
        if Ma>=5
            Cp = Cp_max .* (sin(theta_ramp).^2);    % MNT to find pressure coefficient
        elseif (Ma>=1)&&(Ma<3)
            cal=calF(Ma);
            Cp = Cp_max .* ((1-cal)*sin(theta_ramp).^2+cal);
        elseif (Ma>=3)&&(Ma<5)
            %ADD SUPERSONIC AFTER with CP correction factor
            beta_sonic=asin(sqrt(((gamma-3)*Ma^2+(gamma+1)*Ma^4+sqrt(16*gamma*Ma^4+((gamma-3)*Ma^2+(gamma+1)*Ma^4)^2))/(4*gamma*Ma^4))); %beta_sonic (max shock angle to maintain sonic flow post shock) [rads]
            theta_sonic=atan((2*cot(beta_sonic)*(Ma^2*sin(beta_sonic)^2-1))/(Ma^2*(gamma+cos(2*beta_sonic))+2));
            p_ratio_ost_sonic=1./(1+(gamma-1)./2).^(gamma./(gamma-1)).*(((gamma+1)./2.*(Ma.*sin(beta_sonic)).^2).^(gamma./(gamma-1))./(2.*gamma./(gamma+1).*(Ma.*sin(beta_sonic)).^2-(gamma-1)./(gamma+1)).^(1./(gamma-1)));

            %-----------------------
            %Section of code inspired by:
            %Sambit Supriya Dash (2022). Theta-Beta-MachNo Relation (Plot) for Oblique Shock Waves (https://www.mathworks.com/matlabcentral/fileexchange/72590-theta-beta-machno-relation-plot-for-oblique-shock-waves), MATLAB Central File Exchange. Retrieved June 23, 2022.
            %-----------------------
            %thetas & betas correspond to full range of theta-beta-M relation
            betas = 0:(pi/2/180):(pi/2);   %(0-90deg)
            Nr = ((Ma^2)*((sin(betas)).^2))-1;
            Dr = ((gamma+(cos(2*betas)))*Ma^2)+2;
            thetas = atan(2*cot(betas).*Nr./Dr);
            a = max(thetas);                   % max theta for the Mach No.
            b = betas(find(thetas==a));      % find the beta for max. theta
            %Obtains beta for all theta


            betas(thetas<0)='';
            thetas(thetas<0)='';
            theta_weak=thetas;
            beta_weak=betas;
            theta_weak(beta_weak>=beta_sonic)=''; %eliminates strong shocks
            beta_weak(beta_weak>=beta_sonic)='';
            theta_pos=theta_ramp;
            theta_pos(theta_ramp<0)='';
            p_ratio_ost_sonicpoint=1./(1+(gamma-1)./2.*1.^2).^(gamma./(gamma-1)).*(((gamma+1)./2.*(Ma.*sin(beta_sonic)).^2).^(gamma./(gamma-1))./(2.*gamma./(gamma+1).*(Ma.*sin(beta_sonic)).^2-(gamma-1)./(gamma+1)).^(1./(gamma-1)));
            calibration_sonic=(2.*(p_ratio_ost_sonicpoint-1)./(gamma.*Ma^2.*Cp_max)-cos(pi/2-theta_sonic).^2.)/(sin(pi/2-theta_sonic).^2);
            p_ratio_80=1+2.*gamma./(gamma+1).*(Ma.^2.*sin(interp1(theta_weak,beta_weak,1*pi/180)).^2-1);
            calibration_80=(2.*((p_ratio_80)-1)./(gamma.*Ma^2.*Cp_max)-cosd(89).^2.)/(sind(89).^2);
            
            cal=(calibration_80+calibration_sonic)/2; %Average calibration term
            Cp = Cp_max .* (  (1-cal)*sin(theta_ramp).^2+cal);
        elseif Ma<1
            fprintf('FLOW IS SUBSONIC!!!')
            return
        end
    end

    % Flow properties    
    p0inf_pinf=(1+((gamma-1)/2)*Ma^2)^(gamma/(gamma-1));   % (Pa) Freestream stagnation pressure ratio
    pe_pinf = (1 + ((Cp*gamma*Ma.^2)/2));     % (Pa)  Static Pressure ratio at the edge of the boundary layer
    p_e=p_inf*pe_pinf;

    %Stagnation pressure estimated with Normal shock only
    p0e_p0inf = (((gamma+1)*Ma^2)/(2+(gamma-1)*Ma^2))^(gamma/(gamma-1))*((gamma+1)/((2*gamma*Ma^2)-(gamma-1)))^(1/(gamma-1));   % (Pa) Stagnation pressure ratio at the edge of the boundary layer
    
    p0e_pe=p0e_p0inf*p0inf_pinf./pe_pinf;

    Ma_e = (((p0e_pe).^((gamma-1)/gamma)-1).*(2./(gamma-1))).^(1/2);  % (-) Mach number at the edge of the boundary layer
    T_e = T_0./(1+(((gamma-1)/2).*(Ma_e.^2)));                                    % (K) Temperature at the edge of the boundary layer
    
    
end

function [dS,dydx,theta,k,k_n]=getSurface(rmax,THETA,THETA_cell,npi,npj)
X=-rmax.*cos(THETA);                                                                                                                                                                                                                                                                                                                                                 
Y=rmax.*sin(THETA);
rmax_cell=interp1(THETA, rmax, THETA_cell);
X_cell=-rmax_cell.*cos(THETA_cell);
Y_cell=rmax_cell.*sin(THETA_cell);

%calculate dS (Distance between 2 points on surface) npj-1 dS
dS=zeros(npj-1,1);
for j=1:npj-1
    dS(j)=((X(j+1)-X(j))^2+(Y(j+1)-Y(j))^2)^0.5;
end

%npj gradients
dydx=gradient(Y)./gradient(X);

%Alternate way of calculating gradient (same result)
% dydx(1)=(Y(2)-Y(1))/(X(2)-X(1));
% for j=2:npj-1
%     dydx(j,1)=(Y(j+1)-Y(j-1))/(X(j+1)-X(j-1));
% end
% dydx(npj)=(Y(npj-1)-Y(npj))/(X(npj-1)-X(npj));

theta=atan(dydx); %ramp angle of surface (used to calculate modified newtonian theory)

[k]=calculateCurvature(X,Y);
Xnose=[X(1);X(round(npj/6));X(round(npj/3));X(ceil(npj/2));X(round(2*npj/3));X(round(5*npj/6));X(npj)];
Ynose=[Y(1);Y(round(npj/6));Y(round(npj/3));Y(ceil(npj/2));Y(round(2*npj/3));Y(round(5*npj/6));Y(npj)];
k_n=calculateCurvature(Xnose,Ynose);
k_n=k_n(4);
end

function [q]=blendq(rmax,THETA,npj,dS,q,q_FR)
RB=rmax(npj);    
dist=2*RB; % Smooth the stagnation point heating over entire RB (Body Radius) due to small size
    x=rmax.*cos(THETA);                                                                                                                                                                                                                                                                                                                                                 
    y=rmax.*sin(THETA);
    %Find Stagnation point coord
    x_stag=interp1(THETA,rmax,0);
    y_stag=0;

    %Shorten array to start at point after stagnation
    x=x(round(npj/2)+1:npj);
    y=y(round(npj/2)+1:npj);
    dS=dS(round(npj/2):npj-1);
    dS(1)=dS(1)/2; %dS from stagnation point is half
    q=q(round(npj/2)+1:npj);
    q_max=find(q==max(q));
    %q_FR=q_FR(round(npj/2)+1:npj);
    S=zeros(round(npj/2),1);
    S(1)=dS(1);
    rad_dist=((x-x_stag).^2+(y-y_stag).^2).^0.5; 
    diff=abs(dist-rad_dist);
    ref=find(diff==min(diff));
    for j=2:npj/2
        S(j)=S(j-1)+dS(j);
    end
    for i=1:npj/2
        %Exponential blending function from Heat-3D
        xis=4*(S(i,1)-S(1,1))/(S(ref,1)-S(1,1));

        wfs=1-exp(-0.412*xis^2);

        q(i,1)=q_FR+wfs*(q(i,1)-q_FR);
    end

    q=[flipud(q) ; q];
end

function [ q_cond,q_conv,q_rad, q_sub, T_w, delta_r,mdot,rdot,T_sub] = ...
    heatBalance(rmax, R_air, cp_air, gamma, nu_air, k_air, V_inf, rho_inf, p_inf, Pr, material, npi, npj, p_e, T_e, M_e, dS, T, dt, dth,th,k,r,k_n)
% This function calculates the heat Balance on the surface of the model based on fluid and material parameters
%
% INPUTS:
%   k       :   (m-1)       Curvature of the model
%   k_n     :   (m-1)       Curvature of the model at nose (averaged over
%   entire surface)
%   rmax    :   (m)         Surface radius of the model
%   R_abl   :   (J/kgK)     Specific gas constant of the model material
%   rho_abl :   (kg/m3)     Density of the model material
%   nu_abl  :   (m2/s)      Kinematic viscosity of the model material
%   h_abl   :   (J/kg)      Sublimation enthalpy of the model material
%   material:   ()          model material
%   nu_air  :   (m2/s)      Kinematic viscosity of air
%   k_air   :   (W/mK)      Thermal conductivity of air
%   V_inf   :   (m/s)       Freestream velocity of the air
%   T       :   (K)         Temperature of the freestream
%   T_aw    :   (K)         Adiabatic wall temperature
%   Pr      :   (-)         Prandtl number
%   npi     :   (-)         Radial grid amount
%   npj     :   (-)         Tangential grid amount
%   dt      :   (s)         Time spacing
%   dth     :   (rad)       Angular grid spacing
%   r       :   (m)         Raidal matrix
%   th      :   (rad)       Theta matrix
%   p_e     :   (Pa)        Pressure at BL edge
%   T_e     :   (K)         Temperature at BL edge
%   dS      :   (m)         Surface length from nose

%
% OUTPUTS:
%   qconv_history   :   (W/m2)                  Mean convective heat transfer to the surface
%   qrad_history    :   (W/m2)                  Mean radiative heat transfer from the surface
%   normqe          :   (-)                     Mean normalised conductive heat transfer through the surface
%   normqrad        :   (-)                     Mean normalised radiative heat transfer to the surface
%   normqsub        :   (-)                     Mean normalised sublimative heat transfer to the surface
%   deltar          :   (m)         [npj x 1]   Recession rate of the surface
%   qe              :   (W/m2)      [npj x 1]   Conductive heat transfer to the surface
%   mdot            :   (kg/m2s)    [npj x 1]   Mass transfer at the surface
%   rstag           :
%   T_w             :   (K)                     Wall temperature
    % Nondimensional Parameters for Flow and Heat Transfer
    D=2.*max(rmax);
    T_w=round(T(npi-1,:)',3)+0.001;
    ReD = V_inf * D / nu_air;         % (-) Reynolds Number of the flow (ReD)
    NuD = 0.076 * (ReD.^0.7).*(Pr.^0.37);     % (-) Flow Nusselt Number (Could be changed to local Nu number)
    h = NuD * k_air./D;                      % (W/m2K) Heat Transfer Coefficient based on curvature at theta
    %Sc = 7./(T_w.^0.185);           % (-) Schmidt Number
    T_0e = T_e.*(1+(gamma-1)/2.*M_e.^2);
    V_e = M_e.*((gamma*R_air*T_e).^(1/2));         % (m/s) Velocity at edge of BL
    T_aw = T_e.*(1+Pr.^(1/2)*(gamma-1)*M_e.^(2)/2);
    
    
    %Calculate Sublimation temperature (Equilibrium wall temperature)
    %dependent on Stagnation temperature and pressure
    W_air = 28.97; %(g/mol) molecular weight of air
    %Blowing parameters
    B1=cp_air.*(T_0e(ceil(npj/2))-material.Ts)./(material.h_abl-material.cp_abl*(material.Ts-T(npi-2,ceil(npj/2))));
    B2=material.W_abl/W_air.*((material.Pvps./p_e(ceil(npj/2)))./(1-material.Pvps./p_e(ceil(npj/2))));
    %Find intersection of 2 curves to get equilibrium wall temperature and
    %pressure
    [T_sub,B]=curveIntersection([material.Ts;B1],[material.Ts;B2]);
    T_sub=round(T_sub,3);

%% Convective Heat Flux 
    % Considering boundary layer effects, but ignored entropy layer effects
    % for simplicity
    rho_e = p_e./(R_air.*T_e); % Density at edge of BL
    mu_e = 1.716e-5.*(T_e./273.11).^(3/2).*((273.11 + 110.56)./(T_e + 110.56)); %Dynamic Viscosity using Sutherland's law

    mu_w = 1.716e-5.*(T_w./273.11).^(3/2).*((273.11 + 110.56)./(T_w + 110.56));
    rho_w = p_e./(R_air.*T_w); % Find density at wall (p_e=p_w), for compressible (T_w=T_aw)
    
    h_0e = cp_air.*T_0e; 
    h_w = cp_air.*T_w; 
    
    % Fay-Ridell Convective Heat transfer at Stagnation Point
    %q_FR=0.763/(Pr)^0.6*(rho_e(ceil(npj/2))*mu_e(ceil(npj/2)))^0.4.*(rho_w(ceil(npj/2))*mu_w(ceil(npj/2)))^0.1*(h_0e(ceil(npj/2))-h_w(ceil(npj/2)))*(k(ceil(npj/2)).*sqrt(2*(p_e(ceil(npj/2))-p_inf)./(rho_e(ceil(npj/2))))).^0.5;
    
    q_FR=0.763/(Pr)^0.6*(rho_e(ceil(npj/2))*mu_e(ceil(npj/2)))^0.4.*(rho_w(ceil(npj/2))*mu_w(ceil(npj/2)))^0.1*(h_0e(ceil(npj/2))-h_w(ceil(npj/2)))*(k_n.*sqrt(2*(p_e(ceil(npj/2))-p_inf)./(rho_e(ceil(npj/2))))).^0.5;

    %Laminar Flow convective heat flux, based on Zoby et al.
    [q_conv_lam,theta_mom,Re_mom]=calcqlam(T_w,h_w,T_e,gamma,M_e,p_e,R_air,cp_air,rho_e,mu_e,V_e,Pr,dS,rmax,npj,k);

    %   laminar convective heat transfer calculated first (obtain laminar
    %   momentum thickness)

    %Turbulent flow convective heat flux, based on Zoby et al.
    %   currently not in use as calculation of turbulent boundary layer
    %   momentum thickness requires the use of a numerical solver, slow,
    %   also, laminar flow assumed due to small size of model
%     [q_conv_turb,theta_Tmom,Re_Tmom]=calcqturb(T_w,h_w,T_e,gamma,M_e,p_e,R_air,cp_air,rho_e,mu_e,V_e,Pr,theta_mom,dS,rmax,npj);
    
    %Blending Function would vary convective heating from Peak stagnation
    %(Fay-Riddell) to turbulent convective heating (Zoby et al) Using
    %laminar convective heating
    
    [q_conv]=blendq(rmax,th,npj,dS,q_conv_lam,q_FR);

    % Heat Transfer at each Theta
    h_conv = q_conv./ (T_aw - T_w);             % (W/m2) Convective heat transfer coefficient
    
%% Radiative Heat Flux
    % Radiation Heat Transfer Constants
    sigma = 5.67e-08;
    emiss = 0.9;
    %Stefan-Boltzmann
    q_rad = emiss .* sigma .* T_w.^4;    % (W/m2) Radiative heat transfer

%% Sublimative Heat Flux and Conductive Heat Flux
    % Mass Transfer at each theta
    q_cond=zeros(npj,1);
    q_sub=zeros(npj,1);
    % Partial Vapour Pressure of each material
    if(strcmp(material.mat, "Water Ice"))
        
        W_abl=18.01528; %(g/mol)
    elseif(strcmp(material.mat, "Naphthalene"))
        W_abl=128.1705; %(g/mol)
    elseif(strcmp(material.mat, "Dry Ice"))
        W_abl= 44.01; %(g/mol)

    elseif(strcmp(material.mat, "Camphor"))
        W_abl= 153.23; %(g/mol)
    end
    Pvp=material.Pvp_curve(T_sub);

    %Blowing correction factor (0<CH_CH0<1)
    blowing_corr=1; %=1 for blowing correction;
    if blowing_corr==1
    %--Blowing correction factor from Anderson
    CH_CH0=1./(1+(0.84*(W_air/W_abl)^1/3.*B)); %Blowing correction (for laminar BL from Anderson based on Kubota)
    

    %--Blowing correction factor from Charwat and Gross (A Review of Binary
    %Boundary Layer Characteristics) - Not currently used as use of
    %Reynolds number can cause divergence of solution
%     dS1=dS(round(npj/2):npj-1);
%     s(1)=dS1(1)/2;
%     for j=2:round(npj/2)
%         s(j,1)=s(j-1)+dS1(j);
%     end
%     s=[flipud(s); s]; %Surface distance from stagnation poinht
%     ReS=rho_w.*V_e.*s./mu_e; %surface length dependent Reynolds
%     C_H_0=-q_conv./(rho_e.*V_e.*cp_air.*(T_w-T_aw)); %Zero mass transfer stanton number
%     CH_CH0=1-1.82.*((W_air/W_abl)^1/3.*B2.*C_H_0.*sqrt(ReS)); %Blowing correction (for laminar BL from Charwat)

    else
    CH_CH0=1; %No blowing correction 
    end

    %mdot=B2.*1./CH_CH0.*-q_conv./(cp_air.*(T_w-T_aw));    % (kg/m2s) Mass flow of naphthalene away from the surface
    delta_r=zeros(npj,1);
    mdot=zeros(npj,1);
    rdot=zeros(npj,1);
    for j=1:npj
        if T_w(j)>=T_sub
            T_w(j)=T_sub; %Ensure T_w does not exceed T_sub

            % Heat Transfer at each Theta
            q_cond(j)= material.k_abl*(T(npi-2,j)-T(npi-1,j))/(r(npi-2,j)-r(npi-1,j));%Assume 1-D steady state conduction for conductive heat transfer
            q_sub(j) = q_conv(j)-q_rad(j)-q_cond(j);%mdot(j) .* material.h_abl;                       % (W/m2) Sublimative heat transfer

            % Recession Rate
            mdot(j)=q_sub(j)/(material.h_abl)*CH_CH0;
            delta_r(j) = -mdot(j)*dt/(material.rho_abl);      % (m) Change in radius assuming depth is 1m)
            rdot(j)=-mdot(j)/(material.rho_abl);
        elseif T_w(j)<T_sub %No sublimation if wall temperature is less than sublimation temperature
            %Energy Balance consists of convection, conduction into body and
            %radiation
            q_cond(j)=q_conv(j)-q_rad(j);
            delta_r(j) = 0;      % (m) Change in radius assuming depth is 1m)
            rdot(j)=0;
            mdot(j)=0;
            %Calculate Temperature matrix from given conductive heat flux
        end
        %Solve heat equation for heat flux and temperature boundary conditions
    end

end

function [theta_mom]=calcthetalam(muStar,rhoStar,rho_e,V_e,r,dS,npj,k)
%Calculates momentum thickness (Laminar flow)
var=zeros(length(muStar),1);
f=muStar.*rhoStar.*V_e.*r(1).^2;
for i=round(npj/2)+1:npj %Start from stagnation point
    if i==1
    dvar=(f(i)+0)/2*(dS(i-1,1)/2);
    else
    dvar=(f(i)+f(i-1))/2*(dS(i-1,1));
    end
    if isnan(dvar)
        dvar=0;
    end
    var(i,1)=var(i-1,1)+dvar;
end
var=var+flipud(var);

theta_mom=0.664.*sqrt(var)./(rho_e.*V_e.*r);

end

function [theta_Tmom,m,c1,ReT]=calcthetaturb(muStar,rhoStar,rho_e,V_e,r,dS,theta,mu_e,npj)
%Calculates momentum thickness (Turbulent flow)
%Chop arrays into 2 to start at stagnation point
var=zeros(round(npj/2),1);
rho_e=rho_e(round(npj/2)+1:npj);
muStar=muStar(round(npj/2)+1:npj);
rhoStar=rhoStar(round(npj/2)+1:npj);
V_e=V_e(round(npj/2)+1:npj);
dS=dS(round(npj/2):npj-1);
dS(1)=dS(1)/2;
theta=theta(round(npj/2)+1:npj);
mu_e=mu_e(round(npj/2)+1:npj);
r=r(round(npj/2)+1:npj);
ende=10000; % Number of iterations
for c=1:ende
    if c==1
        ReT=rho_e.*V_e.*theta./mu_e;
        theta_Tmom=theta;
    end

    N= 12.67 - 6.5.*log10(ReT(:,c)) + 1.21.*(log10(ReT(:,c))).^2; %use Re_mom from laminar calculation
    m = 2./(N+1);
    c3=1+m;
    c4=1./c3;
    c5=2.2433 + 0.93.*N;
    c1=(1./c5).^(2*N./(N+1)).*(N./((N+1).*(N+2))).^m;
    c2=(1+m).*c1;
    var(1)=0;
    f=muStar.^m .*rhoStar.*V_e.*r.^c3;
    for i=2:round(npj/2)
        dvar=(f(i)+f(i-1))/2*(dS(i,1));
        if isnan(dvar)
            dvar=0;
        end
        var(i,1)=var(i-1,1)+dvar;
    end
    %var=var+flipud(var);
    intT=(c2.*var).^c4./(rho_e.*V_e.*r);
    relax=0.01;
    theta_Tmom(:,c+1) = ((1-relax)*theta_Tmom(:,c) + relax*intT);
    ReT(:,c+1)=rho_e.*V_e.*theta_Tmom(:,c)./mu_e;
    if (c>2)
        resnorm = normalizedChangeOfMultiple(theta_Tmom,ReT);
        if(resnorm <= 1e-6)
            break;
        end
    end
end
theta_Tmom=theta_Tmom(:,end);
theta_Tmom(isnan(theta_Tmom))=0;
ReT=ReT(:,end);
%theta_Tmom=(c2.*var).^c4./(rho_e.*V_e.*r);
theta_Tmom=[flipud(theta_Tmom); theta_Tmom];
ReT=[flipud(ReT); ReT];
c1=[flipud(c1); c1];
m=[flipud(m); m];
    function resnorm = normalizedChangeOfMultiple(varargin)
        resnorm = 0;
        for iterResNorm=1:nargin
            resnorm = max(resnorm,calcNormalizedChange(varargin{iterResNorm}));
        end
    end

    function resnorm = calcNormalizedChange(data)
        resnorms = abs((data(:,end)-data(:,end-1))./data(:,2));
        resnorm = max(resnorms);
    end
end

function [theta_Tmom,m,c1]=calcthetaturb2(muStar,rhoStar,rho_e,V_e,r,dS,theta_mom,npj)
%This function calculates momentum thickness of turbulent boundary layer
%with a simplified assumption that the laminar momentum thickness can be used to calculate the coefficients for obtaining the turbulent momentum thickness 
%Chop arrays into 2 to start at stagnation point
var=zeros(round(npj/2),1);
rho_e=rho_e(round(npj/2)+1:npj);
muStar=muStar(round(npj/2)+1:npj);
rhoStar=rhoStar(round(npj/2)+1:npj);
V_e=V_e(round(npj/2)+1:npj);
dS=dS(round(npj/2):npj-1);
dS(1)=dS(1)/2;
theta_mom=theta_mom(round(npj/2)+1:npj);
r=r(round(npj/2)+1:npj);
%Calculates momentum thickness (Laminar flow)
Re_mom=rho_e.*V_e.*theta_mom./muStar; %Reynolds estimated using laminar momentum thickness
var(1)=0;
N= 12.67 - 6.5.*log10(Re_mom) + 1.21.*(log10(Re_mom)).^2; %use Re_mom from laminar calculation
m = 2./(N+1);
c3=1+m;
c4=1./c3;
c5=2.2433 + 0.93.*N;
c1=(1./c5).^(2*N./(N+1)).*(N./((N+1).*(N+2))).^m;
c1(1)=0;
c2=(1+m).*c1;
f=muStar.^m .*rhoStar.*V_e.*r.^c3;
for i=2:round(npj/2)
    dvar=(f(i)+f(i-1))/2*(dS(i,1));
    if isnan(dvar)
        dvar=0;
    end

    var(i,1)=var(i-1,1)+dvar;
end

theta_Tmom=(c2.*var).^c4./(rho_e.*V_e.*r);
theta_Tmom=[flipud(theta_Tmom); theta_Tmom];
m=[flipud(m); m];
c1=[flipud(c1); c1];
end

function [q_conv_lam,theta_mom,Re_mom]=calcqlam(T_w,h_w,T_e,gamma,M_e,p_e,R_air,cp_air, rho_e,mu_e,V_e,Pr,dS,rmax,npj,k)
%Calculates Laminar Convective heat transfer and laminar momentum thickness
%Ekhart Reference Temperature
TStar=T_e.*(0.45+0.55.*T_w./T_e+0.16.*sqrt(Pr).*(gamma-1)/2.*M_e.^2); %Reference Temperature
muStar=1.716e-5.*(TStar./273.11).^(3/2).*((273.11 + 110.56)./(TStar + 110.56));
rhoStar = p_e./(R_air*TStar);

%Momentum Thickness
theta_mom=calcthetalam(muStar,rhoStar,rho_e,V_e,rmax,dS,npj,k);
Re_mom=rho_e.*V_e.*theta_mom./mu_e;

%Adiabatic wall Temp (Laminar)
h_aw=cp_air*T_e + sqrt(Pr).*V_e.^2./2;
%Heat calculation from Zoby et al. (1979)
q_conv_lam=0.22./Re_mom.*(rhoStar./rho_e).*(muStar./mu_e).*rho_e.*V_e.*(h_aw-h_w).*Pr.^(-0.6);
end

function [q_conv_turb,theta_Tmom,Re_Tmom]=calcqturb(T_w,h_w,T_e,gamma,M_e,p_e,R_air,cp_air,rho_e,mu_e,V_e,Pr,theta_mom,dS,rmax,npj)
%Calculates Turblulent Convective heat transfer and turbulent momentum thickness
%Ekhart Reference Temperature
TStar=T_e.*(0.5.*(1+T_w./T_e)+0.16.*(Pr)^(1/3).*(gamma-1)/2.*M_e.^2); %Reference Temperature
muStar=1.716e-5.*(TStar./273.11).^(3/2).*((273.11 + 110.56)./(TStar + 110.56));
rhoStar = p_e./(R_air*TStar);

%Momentum Thickness
[theta_Tmom,m,c1,Re_Tmom]=calcthetaturb(muStar,rhoStar,rho_e,V_e,rmax,dS,theta_mom,mu_e,npj);
[theta_Tmom2,m2,c12]=calcthetaturb2(muStar,rhoStar,rho_e,V_e,rmax,dS,theta_mom,npj);
Re_Tmom2=rho_e.*V_e.*theta_Tmom2./mu_e;

%Adiabatic wall Temp (Turbulent)
h_aw=cp_air*T_e + Pr.^(1/3).*V_e.^2./2;
%Heat calculation from Zoby et al. (1979)
q_conv_turb=c1.*Re_Tmom.^(-m).*(rhoStar./rho_e).*(muStar./mu_e).^m.*rho_e.*V_e.*(h_aw-h_w).*Pr.^(-0.4);
%From alternative function not solving for momentum thickness iterately
q_conv_turb2=c1.*Re_Tmom2.^(-m).*(rhoStar./rho_e).*(muStar./mu_e).^m.*rho_e.*V_e.*(h_aw-h_w).*Pr.^(-0.4);
end

function [aN, aE, aS, aW, aP, bP, sP] = Tcoeff(dth, dr, r, c, material, dt, T, rmax, q_cond, npi, npj,T_sub)
    %Define Coefficients for N,E,S,W)

    aN = zeros(npi, npj); %varies from 2 to i and 1 to j as ignore points at origin (r=0 as no arclength)
    aE = zeros(npi, npj); 
    aS = zeros(npi, npj); 
    aW = zeros(npi, npj);
    aP = zeros(npi, npj);
    sP = zeros(npi, npj);
    bP = zeros(npi, npj);
    
    %Determine interior Temp coefficients
    for i=2:npi-1
        for j=2:npj-1
            %-------Specify coefficient matrices
        aE(i,j) = dth/dr(j)*(r(i,j)+dr(j)/2);
        aW(i,j) = dth/dr(j)*(r(i,j)-dr(j)/2);
        aN(i,j) = (dr(j+1)/(1+c(j+1)*r(i,j+1))+dr(j)/(1+c(j)*r(i,j)))/(2*dth*(1+c(j)*r(i,j)));
        aS(i,j) = (dr(j-1)/(1+c(j-1)*r(i,j-1))+dr(j)/(1+c(j)*r(i,j)))/(2*dth*(1+c(j)*r(i,j)));
        sP(i,j) = dr(j)*dth*r(i,j)/(material.alpha_abl*dt);
        aP(i,j) = aW(i,j)+aE(i,j)+aN(i,j)+aS(i,j)+sP(i,j);
        bP(i,j) = sP(i,j)*T(i,j);

        end
    end
    % EAST & WEST BC
    for j=2:npj-1
        % Eastern BC - Incoming Heat flux or Constant Temperature

        if T(npi-1,j)<T_sub
            %set constant heat flux condition
            aE(npi-1,j)=0;
            bP(npi-1,j)=sP(npi-1,j)*T(npi-1,j)+rmax(j)*dth*q_cond(j)/material.k_abl;
        else
            %set constant temperature condition
            aE(npi-1,j) = dth/dr(j)*(r(npi-1,j)-dr(j)/2);
        end
        aP(npi-1,j) = aE(npi-1,j)+aN(npi-1,j)+aS(npi-1,j)+aW(npi-1,j)+sP(npi-1,j);
        % Western BC - Constant Temperature adiabatic
        aW(2,j)=0;
        aP(2,j)=aW(2,j)+aN(2,j)+aS(2,j)+aE(2,j)+sP(2,j);
    end
    % NORTH & SOUTH BC
    for i = 2:npi-1
        % Northern BC
        aN(i,npj-1)=0; %
        aP(i,npj-1)=aN(i,npj-1)+aS(i,npj-1)+aE(i,npj-1)+aW(i,npj-1)+sP(i,npj-1);

        % Southern BC
        aS(i,2)=0;
        aP(i,2)=aN(i,2)+aE(i,2)+aW(i,2)+aS(i,2)+sP(i,2);
    end
end

function [phi] = gs_solver(phi, aN, aE, aS, aW, aP, bP, npi, npj,T_sub)
for i = 2:npi-2 %Gauss siedel for interior cells
    for j = 2:npj-1
        if phi(i,j)>=T_sub
            phi(i,j)=T_sub;
        else
            phi(i,j) = (aE(i,j)*phi(i+1,j)+aW(i,j)*phi(i-1,j)+aN(i,j)*phi(i,j+1)+aS(i,j)*phi(i,j-1)+bP(i,j))/aP(i,j);
        end
    end
end
phi(:,1)=phi(:,2);
phi(:,npj)=phi(:,npj-1);

for j=2:npj-1 %Solve considering wall condition
    if phi(npi-1,j)<T_sub
        %Neumann BC
        phi(npi-1,j) = (aW(npi-1,j)*phi(npi-2,j)+aN(npi-1,j)*phi(npi-1,j+1)+aS(npi-1,j)*phi(npi-1,j-1)+bP(npi-1,j))/aP(npi-1,j); % EAST
    else
        %Dirichlet BC
        phi(npi-1,j)=T_sub;
    end
end
%Set south wall temperature T(1) to equal T(2)
%Improves solution for high
%conductivity materials with long runtimes,helps saturate the internal
%temperature close to the internal boundary
phi(:,1) = phi(:,2);
phi(:,npj) = phi(:,npj-1);
end

function [phi] = gs_solver_even(phi, aN, aE, aS, aW, aP, bP, npi, npj,T_sub) % run solver opposite direction on even steps to improve stability.
for i = 2:npi-2 %Gauss siedel for interior cells
    for j = npj-1:-1:2
        if phi(i,j)>=T_sub
            %Dirichlet BC 
            phi(i,j)=T_sub;
        else
            %Neumann BC
            phi(i,j) = (aE(i,j)*phi(i+1,j)+aW(i,j)*phi(i-1,j)+aN(i,j)*phi(i,j+1)+aS(i,j)*phi(i,j-1)+bP(i,j))/aP(i,j);
        end
    end
end


for j = npj-1:-1:2 %Solve considering wall condition
    if phi(npi-1,j)<T_sub
        phi(npi-1,j) = (aW(npi-1,j)*phi(npi-2,j)+aN(npi-1,j)*phi(npi-1,j+1)+aS(npi-1,j)*phi(npi-1,j-1)+bP(npi-1,j))/aP(npi-1,j); % EAST
    else
        phi(npi-1,j)=T_sub; 
    end
end
%Set south wall temperature T(1) to equal T(2) 
%Improves solution for high
%conductivity materials with long runtimes,helps saturate the internal
%temperature close to the internal boundary
phi(:,1)=phi(:,2);
phi(:,npj)=phi(:,npj-1);
end

function norm = norm2matrix(T11, T22, dr, dth)
    ni = size(T11,1);
    nj = size(T11,2);

    norm0 = 0;
    for i = 2:ni-1
        for j = 2:nj-1
            norm0 = norm0 + ((T11(i,j) - T22(i,j))^2);
        end
    end
    norm = sqrt(dr*(floor(ni/2)*dth*norm0));
end

function [c, dr, r, rmax] = remesh(r, th, rmax, rmin, c, ~, deltar, npi, npj)
rmax = rmax + deltar;
for j = 2:npj-1
    x1 = rmax(j-1)*cos(th(j-1));
    y1 = rmax(j-1)*sin(th(j-1));
    x2 = rmax(j)*cos(th(j));
    y2 = rmax(j)*sin(th(j));
    x3 = rmax(j+1)*cos(th(j+1));
    y3 = rmax(j+1)*sin(th(j+1));

    c1 = x1*(y2-y3) - y1*(x2-x3) + x2*y3 - x3*y2;
    c2 = (x1^2+y1^2)*(y3-y2) + (x2^2+y2^2)*(y1-y3) + (x3^2+y3^2)*(y2-y1);
    c3 = (x1^2+y1^2)*(x2-x3) + (x2^2+y2^2)*(x3-x1) + (x3^2+y3^2)*(x1-x2);
    c4 = (x1^2+y1^2)*(x3*y2-x2*y3) + (x2^2+y2^2)*(x1*y3-x3*y1) + (x3^2+y3^2)*(x2*y1-x1*y2);

    c(j) = sqrt((4*c1^2)/(c2^2+c3^2-4*c1*c4)); %Calculate local curvature
end
c(1) = c(2);
c(npj) = c(npj-1);
dr = (rmax - rmin)/(npi-2);
r(1,:) = rmin;
r(2,:) = rmin + 0.5*dr';
for i = 3:npi-1
    r(i,:) = r(i-1,:) + dr';
end
r(npi,:) = r(npi-1,:) + 0.5*dr';
end

function [r, th, rmax, c, dr, dth] = create_grid(npi, npj, rin, rmin, shape)
% This function creates a gird in polar coordinates of a model that will be placed in a wind tunnel
%
% INPUTS:
%   npi     :   (-)     Number of spaces in the radial direction
%   npj     :   (-)     Number of spaces in the angular or tangential direction 
%   rin     :   (m)     Initial outer radius of the model
%   rmin    :   (m)     Minimum radius of the model
%   shape   :   (-)     Shape parameter to diffrentiate between different shape options
%
% OUTPUTS:
%   r       :   (m)     [npi x npj]     Coordinate field of the model   
%   th      :   (rad)   [npj x 1]       Angular spacing of the model
%   rmax    :   (m)     [npj x 1]       Surface radius
%   c       :   (m-1)   [npj x 1]       Surface curvature
%   dr      :   (m)                     Radial grid spacing
%   dth     :   (rad)                   Angular grid spacing

    % Preallocation of arrays
    r = zeros(npi, npj);                    % (m) Spational Location
    th = zeros(npj, 1);                     % (rad) Anglular location

    % Definition of rmax array
    rmax = ones(npj,1) .* rin;              % (m) Maximum Radial Distance

    % Setting up angular space of the solution
    thmax = 0.5 * pi;                       % (rad) Maximum value of theta (pi/2)
    thmin = -0.5 * pi;                      % (rad) Minimum value of theta (-pi/2)
    dth = (thmax - thmin)/(npj-2);          % (rad) Grid spacing of theta

    % Defining Angular spacing
    th(1) = thmin;                          
    th(2) = thmin + 0.5 .* dth;

    for j = 3:npj-1
        th(j) = th(j-1) + dth;              % Following theta = prev theta + delta theta 
    end

    th(npj) = th(npj-1) + 0.5.*dth;         % Final theta has different spacing

    % Altering rmax if shape is an elipse instead of a circle
    if(strcmp(shape,"Ellipse"))
        disp("Ellipse Shape selected")
        rmaxmin = 0.75 * rin;                   % (m) Semi-minor axis of ellipse
        e = sqrt(1 - (rmaxmin.^2)/(rin.^2));    % (-) Eccentricity of ellipse
        for k = 1:npj
            rmax(k) = rmaxmin ./ (sqrt(1 - ((e.*cos(th(k))).^2)));
        end
    end

    % Grid spacing in r and theta
    dr = (rmax-rmin)/(npi-2);               % (m) Grid spacing of r
    r(1,:) = rmin;
    r(2,:) = rmin + 0.5.*dr;
    for i = 3:npi-1
        r(i,:) = r(i-1,:) + dr';
    end
    r(npi,:) = r(npi-1,:) + 0.5 .* dr';

    % Calculate the curvature of the shape
    c = ones(npj, 1) ./ rmax;
end

function [k] = calculateCurvature(X,Y)
%     Code adapted from MATLAB File Exchange Moreno, M. 
%     https://uk.mathworks.com/matlabcentral/fileexchange/107929-discrete-curvature-normals-and-evolute
%---------------------------------------------------------------------------
%     [k] = calculateCurvature(X, Y) Returns the discrete curvature (1D),
%     in the Euclidean space.
%       Inputs:
%       X -> X coordinates of outer surface
%       Y -> Y coordinates of outer surface
% Parse input arguments
    XY = [X,Y];
    [N, p] = size(XY);
    if N < p
        XY = XY';
        [N, p] = size(XY);
    end

    % Obtain pairs of combinations between adjacent triplets
    XY = [XY, zeros(N, 3 - p)];
    A = XY(1 : N - 2, :) - XY(3 : N, :);
    B = XY(2 : N - 1, :) - XY(3 : N, :);
    XY = XY(:, 1 : p);
    % Pre-calculate cross-products
    I = cross(A, B, 2);
    J = cross(B, I, 2);
    K = cross(A, I, 2);
    % Obtain relative circumcentre
    A = A .* A;
    I = I .* I;
    O = (sum(A, 2) .* J - sum(B .* B, 2) .* K) ./ sum(I, 2) / 2;
    O = O(:, 1 : p);
    % Initialise extrapolating coefficients
    I = (XY(2, :) - XY(1, :)) ./ (XY(3, :) - XY(2, :));
    J = (XY(N, :) - XY(N - 1, :)) ./ (XY(N - 1, :) - XY(N - 2, :));
    
    % Calculate end-point extrapolated values
    i = I .* (O(2, :) - O(1, :));
    j = J .* (O(N - 2, :) - O(N - 3, :));
    C = [O(1, :) - i; O; O(N - 2, :) + j];
    if norm(XY(N, :) - XY(1, :)) < 1e-6
        C([1, N], :) = [1; 1] .* (C(1, :) + C(N, :)) / 2;
    end
    % Derive 1-D curvature from radii
    R = sqrt(sum(C .* C, 2));
    i = isnan(R);
    R(i) = Inf;
    k = 1 ./ R;
    
end

function [X,Y] = curveIntersection(L1,L2)
%     Code adapted from MATLAB File Exchange NS 
%     https://www.mathworks.com/matlabcentral/fileexchange/22441-curve-intersections
%---------------------------------------------------------------------------
%   Intersection of curves
%   [X,Y] = curveIntersection(L1,L2) returns the intersection points of two curves L1 
%   and L2. The curves L1,L2 can be either closed or open and are described
%   by two-row-matrices, where each row contains its x- and y- coordinates.
%   The intersection of groups of curves (e.g. contour lines, multiply 
%   connected regions etc) can also be computed by separating them with a
%   column of NaNs as for example
%
%         L  = [x11 x12 x13 ... NaN x21 x22 x23 ...;
%               y11 y12 y13 ... NaN y21 y22 y23 ...]
%
%   P has the same structure as L1 and L2, and its rows correspond to the
%   x- and y- coordinates of the intersection points of L1 and L2. If no
%   intersections are found, the returned P is empty.
%
%   Author : NS
%   Version: 3.0, 21 Sept. 2010
%   Two words about the algorithm: Most of the code is self-explanatory.
%   The only trick lies in the calculation of C1 and C2. To be brief, this
%   is essentially the two-dimensional analog of the condition that needs
%   to be satisfied by a function F(x) that has a zero in the interval
%   [a,b], namely
%           F(a)*F(b) <= 0
%   C1 and C2 exactly do this for each segment of curves 1 and 2
%   respectively. If this condition is satisfied simultaneously for two
%   segments then we know that they will cross at some point. 
%   Each factor of the 'C' arrays is essentially a matrix containing 
%   the numerators of the signed distances between points of one curve
%   and line segments of the other.
    %...Argument checks and assignment of L2

    hF = @le;
   
       
    %...Preliminary stuff
    x1  = L1(1,:)';  x2 = L2(1,:);
    y1  = L1(2,:)';  y2 = L2(2,:);
    dx1 = diff(x1); dy1 = diff(y1);
    dx2 = diff(x2); dy2 = diff(y2);
    
    %...Determine 'signed distances'   
    S1 = dx1.*y1(1:end-1) - dy1.*x1(1:end-1);
    S2 = dx2.*y2(1:end-1) - dy2.*x2(1:end-1);
    
    C1 = feval(hF,D(bsxfun(@times,dx1,y2)-bsxfun(@times,dy1,x2),S1),0);
    C2 = feval(hF,D((bsxfun(@times,y1,dx2)-bsxfun(@times,x1,dy2))',S2'),0)';
    %...Obtain the segments where an intersection is expected
    [i,j] = find(C1 & C2); 
    if isempty(i),P = zeros(2,0);return; end;
    
    %...Transpose and prepare for output
    i=i'; dx2=dx2'; dy2=dy2'; S2 = S2';
    L = dy2(j).*dx1(i) - dy1(i).*dx2(j);
    i = i(L~=0); j=j(L~=0); L=L(L~=0);  %...Avoid divisions by 0
    
    %...Solve system of eqs to get the common points
    P = unique([dx2(j).*S1(i) - dx1(i).*S2(j), ...
                dy2(j).*S1(i) - dy1(i).*S2(j)]./[L L],'rows')';
    X=P(1);
    Y=P(2); 
    function u = D(x,y)
        u = bsxfun(@minus,x(:,1:end-1),y).*bsxfun(@minus,x(:,2:end),y);
    end
end

function plotResults(histT,npi,npj,th,histr,timeIter,dt,T,T_in,rin,material,P_0,T_0,Ma, histrdot)
%
%**** Purpose: Plotting time history video of internal temperature field
%and shape
%
name =char(strjoin(['Ablation' material.mat ' P_0=' P_0 ' T_0=' T_0 ' M=' Ma '.avi'])); %Video file name
vidObj = VideoWriter(name);
vidObj.FrameRate = 30;
open(vidObj);
rin_mm=ceil(rin*1000);
histT(:,npi,:) = histT(:,npi-1,:);
Thgrid = zeros(npi,npj);
for i = 1:npi
    Thgrid(i,:) = th';
end
figure(10)
nt=100; %Number of timesteps to show;
for i = 1:round(timeIter/nt):timeIter
    Rgrid = squeeze(histr(i,:,:));
    X = Rgrid.*cos(Thgrid).*1000;
    Y = Rgrid.*sin(Thgrid).*1000;
    clf
    surface(X,Y,squeeze(histT(i,:,:)))
    hold on
    set(gca,'FontSize',12,'TickLabelInterpreter','latex')
    xlim([0 rin_mm])
    ylim([-rin_mm rin_mm])
    view(2)
    daspect([1 1 1])
    xlabel('$x$ [mm]','fontsize',12,'interpreter', 'latex')
    ylabel('$y$ [mm]','fontsize',12,'interpreter', 'latex')
    shading interp
    lighting phong
    colorbar
    caxis([T_in max(max(T))])
    set(colorbar,'TickLabelInterpreter','latex')
    text(rin_mm*1.05,rin_mm*1.1,'T [K]','FontSize',12,'interpreter', 'latex')
    text(1.4*rin_mm,0.9*rin_mm,sprintf('t = %.2fs',i*dt),'interpreter', 'latex')
    text(1.4*rin_mm,-0.9*rin_mm,['$$\stackrel{.}{r}_{N}$$' convertStringsToChars(sprintf(' = %.5fmm/s',histrdot(i)*1000))],'interpreter', 'latex')
    pause(0.00001)
    mov = getframe(gcf); %add current figure as frame to mov
    writeVideo(vidObj,mov); %write frame to vidObj
end
% Close the file.
close(vidObj);
end