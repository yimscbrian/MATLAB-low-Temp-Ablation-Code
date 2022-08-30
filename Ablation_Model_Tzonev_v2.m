% AUTHOR: PHILIP TZONEV (based on the code of Laurence Maskell)
% OBJECTIVE: Calculation of temperature distribution and recession rates 
% for a low-temperature sublimator heat shield model.

% HOUSEKEEPING
clear
clc
close all
set(0,'defaulttextinterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
tic

%% USER INPUT
[material, shape, temp, press, speed, npi, ~, Tend, dt] = userInput();

%% INTIALISATION
% Solution Parameters
eps = 1e-6;                     % (-) Solution Tolerance
maxiter = 20000;                % (-) Maximum number of iterations for Gauss-Seidel Solver

% Spatial Discretisation
npi = str2double(npi);          % (-) Cells in the radial direction
npj = round(pi*npi);            % (-) Cells in the tangential direction. The current set-up allows for equal sizes 

% Time Discretisation
Tend = str2double(Tend);        % (s) End time
dt = str2double(dt);            % (s) Time step
nt = round(Tend/dt);            % (-) Number of time steps

% Freestream Conditions (Wind Tunnel)
P_s = str2double(press);        % (kPa) Stagnation Pressure of the flow    
T_s = str2double(temp);         % (K) Stagnation temperature of the flow
Ma = str2double(speed);         % (-) Mach number of the flow

% Heatshield Geometry
rin = 0.05;                    % (m) Initial OR of heat shield
rmin = 0.1 * rin;               % (m) Inner radius of heat shield (not a physical feature, but to provide a limit for the solver)
T_in = 293;                     % (K) Initial temperature of the heatshield
input = 0;                      % (-) Parameter that control initial temperature distribution

% Properties of Air
gamma_air = 1.4;                % (-) Ratio of specific heats of air
R_air = 287;                    % (J/kgK) Gas constant for air
cp_air = 1004;                  % (J/kgK) Specific Heat Capicity for Air assuming constant pressure
rho_air = 0.29;                 % (kg/m3) Density of Air
nu_air = 1.568e-05;             % (m2/s) Kinematic viscosity of air
k_air = 0.025;                  % (W/mK) Thermal conductivity of air
k_sg = 1.7415e-04;              % (kg0.5/m) Sutton-Graves constant
Pr = 0.7;                       % Prandtl Number

% Properties of the model material
if(strcmp(material, "Water Ice"))
    R_nap = 461.52;                             % (J/kgK) Gas constant for naphthalene
    cp_nap = 4182;                              % (J/kgK) Specific Heat Capicity for napthalene assuming constant pressure
    rho_nap = 917;                              % (kg/m3) Density of naphthalene in solid form
    nu_nap = 1e-06;                             % (m2/s) Kinematic viscosity of naphthalene
    k_nap = 0.598;                              % (W/mK) Thermal conductivity of naphthalene
    h_nap = 2.8386e06;                          % (J/kg) Enthalpy of Sublimation of naphthalene
    alpha_nap = k_nap / (rho_nap * cp_nap);     % (m2/s) Thermal Diffusivity of naphthalene in solid form
elseif(strcmp(material, "Naphthalene"))
    R_nap = 64.87;                              % (J/kgK) Gas constant for naphthalene
    cp_nap = 1300;                              % (J/kgK) Specific Heat Capicity for napthalene assuming constant pressure
    rho_nap = 1140;                             % (kg/m3) Density of naphthalene in solid form
    nu_nap = 9.689e-07;                         % (m2/s) Kinematic viscosity of naphthalene
    k_nap = 0.1221;                             % (W/mK) Thermal conductivity of naphthalene
    h_nap = 5.57e05;                            % (J/kg) Enthalpy of Sublimation of naphthalene
    alpha_nap = k_nap / (rho_nap * cp_nap);     % (m2/s) Thermal Diffusivity of naphthalene in solid form
else
    R_nap = 64.87;                              % (J/kgK) Gas constant for naphthalene
    cp_nap = 1300;                              % (J/kgK) Specific Heat Capicity for napthalene assuming constant pressure
    rho_nap = 1000;                             % (kg/m3) Density of naphthalene in solid form
    nu_nap = 9.689e-07;                         % (m2/s) Kinematic viscosity of naphthalene
    k_nap = 0.34;                               % (W/mK) Thermal conductivity of naphthalene
    h_nap = 5.57e05;                            % (J/kg) Enthalpy of Sublimation of naphthalene
    alpha_nap = k_nap / (rho_nap * cp_nap);     % (m2/s) Thermal Diffusivity of naphthalene in solid form
end

%% Create and plot the solution space
[r, th, rmax, c, dr, dth] = create_grid(npi, npj, rin, rmin, shape);    % Creates the solution space of the model
[X,Y] = grid_plotting(r, th, rmax, npi, npj, 1);                        % Plots the model and generated grid

% Preallocation of various arrays
normqe = zeros(nt,1);               % (-) [nt x 1]          Normalised Conductive Heat
normqrad = zeros(nt,1);             % (-) [nt x 1]          Normalised Radiative Heat
normqsub = zeros(nt,1);             % (-) [nt x 1]          Normalised Sublimative Heat
qconv_history = zeros(nt,1);        % (W/m2K) [nt x 1]          Normalised Conductive Heat
deltar = zeros(nt,npj);             % (m) [nt x npj]        Recession at each timestep
histr = zeros(nt, npi, npj);        % (m) [nt x npi x npj]  History of radial grid
histrmax = zeros(nt, npj);          % (m) [nt x npj]        History of outer surface dimension
histT = zeros(nt, npi, npj);        % (K) [nt x npi x npj]  History of temperature field
histc = zeros(nt, npj);             % (m-1) [nt x npj]      History of the curvature
mdotstag = zeros(nt, npj);          % (kg/m2s) [nt x npj]        History of surface mass flow

%% STEP 1. MODIFIED NEWTONIAN THEORY AND STAGNATION TEMPERATURE
% Calculating Freestream quantities
T_inf = T_s/(1+(Ma^2*gamma_air*R_air)/(2*cp_air));  % (K) Freestream temperature
V_inf = Ma*((gamma_air*R_air*T_inf)^(1/2));         % (m/s) Freestream velocity
rf = Pr^(1/3);                                      % (-) Recovery factor of turbulent flow

[Cp, Ma_ee, p_ee, T_ee] = modified_newtonian_theory(th, Ma, gamma_air, T_inf, V_inf, rho_air, R_air, cp_air);

T_aw = T_ee.*(1+rf*(gamma_air-1)*Ma_ee.^(2)/2); %pg 508 Viscous Fluid Flow Frank White

% Initialise the temperature field and boundary conditions
T = initial_temp_distribution(T_in, input, npi, npj);
%% TIME LOOP (STEPS 2 - 5)
for timeIter = 1:nt
    % STEP 2. SURFACE HEATING AND MASS FLOW
    [qconv_temp, normqe_temp, normqrad_temp, normqsub_temp, deltar_temp, qe, mdot] = ...
        heatIn(c, rmax, R_nap, rho_nap, nu_nap, h_nap, nu_air, k_air, V_inf, T, T_aw, Pr, npi, dt, dth);
    % Storing results in history arrays
    normqe(timeIter) = normqe_temp;
    normqrad(timeIter) = normqrad_temp;
    normqsub(timeIter) = normqsub_temp;
    qconv_history(timeIter) = qconv_temp;
    deltar(timeIter,:) = deltar_temp';

    % STEP 3. TEMPERATURE COEFFICIENT MATRICIES
    [aN, aE, aS, aW, aP, bP, sP] = Tcoefficiants(dth, dr, r, c, alpha_nap, dt, T, rmax, qe, k_nap, npi, npj);

    % STEP 4. THE GAUSS-SEIDEL SOLVER
    for iter = 1:maxiter
        if(iter == 1)
            % In the first iteration, there is no previously calculated temperature field, and so the current temperature field is provided before it is updated
            Told = T; 
        end

        T = gs_solver(T, aN, aE, aS, aW, aP, bP, npi, npj);     % Solving for the temperature field using a Gauss-Seidel Solver

        if(iter == 1)
            norm0 = norm2matrix(Told, T, dr, dth);              % Computing the norm of the two temperature fields
            Told = T;                                           % Overwriting the previous temperature field
        else
            normNext = norm2matrix(Told, T, dr, dth);           % Computing the norm of the two temperature fields
            if(normNext < (eps*norm0))
                % If the calculated norm is smaller than the applied tolereance, solution has converged
                break
            else
                % If the calculated norm is greater than the applied tolerance, another iteration is computed
                Told = T;   
            end
        end
    end

    % STEP 6. SAVE SURFACE RECESSION AMOUNT TO HISTORY MATRIX
    histr(timeIter, :, :) = r;
    histrmax(timeIter, :) = rmax;
    histT(timeIter, :, :) = T;
    histc(timeIter, :) = c;
    mdotstag(timeIter, 1) = mdot(npj/2)./dt;        % Surface recession at nose
    mdotstag(timeIter, 2) = mdot(npj/3)./dt;        % Surface recession at theta = 30 degrees
    mdotstag(timeIter, 3) = mdot(npj/6)./dt;        % Surface recession at theta = 60 degrees

    % STEP 5. REMESH GRID USING CURVATURE OF THE CICRLCE
    [c, dr, r, rmax] = remesh(r, th, rmax, rmin, c, dr, deltar_temp, npi, npj);

    % CHECK FOR MIMIMUM r-VALUE
    if(max(rmax)< 0.2 * rin)
        % Checks if the heat shield has completely sublimated away. If so, program exits
        disp('Reached mimimum radius, program aborted');
        break
    elseif(c(round(npj/2)) < 10)
        % Checks if the heat shield is no` longer circular enough. If so, program exits
        disp('Reached mimimum allowable curvature, program aborted');
        break
    end
    fprintf('Number of iterations: %i \t Time: %f\n',iter,timeIter*dt);
end

toc

%% PLOTTING OF RESULTS

[Xend,Yend] = grid_plotting(r, th, rmax, npi, npj, 2);      % Plotting of model at final timestep
time = [1:timeIter] * dt;                                   % Time simulated successfully by model (present in case the model exits early due to a failure condition)

figure(4)   % FIGURE 4 - COMPARING PROPORTIONS OF ENTRY HEATING (CONDUCTION VS RADIATION VS SUBLIMATION)
plot((time),normqe(1:timeIter),  '-r', 'LineWidth', 4)
hold on
grid on
plot(time, normqrad(1:timeIter), '-k', 'LineWidth', 4)
plot(time, normqsub(1:timeIter), '-b', 'LineWidth', 4)
legend('q_{cond}', 'q_{rad}', 'q_{sub}', 'Interpreter', 'tex', 'Location', 'East','FontSize', 24)
ylabel('$\stackrel{.}{q}$/$\stackrel{.}{q}_{conv}$', 'Interpreter', 'latex', 'FontSize', 24)
xlabel('t [s]', 'FontSize', 24)
set(gca, 'FontSize', 14)

figure(5)   % FIGURE 5 - SURFACE MASS FLOW RATES AT VARIOUS ANGLES
plot(time, mdotstag(1:timeIter, 1), 'LineWidth', 2)
hold on
grid on
plot(time, mdotstag(1:timeIter, 2), 'LineWidth', 2)
plot(time, mdotstag(1:timeIter, 3), 'LineWidth', 2)
ylabel('$\stackrel{.}{m}$ [kg/$m^2$s]', 'Interpreter', 'latex')
xlabel('t [s]')
legend('\theta = 0\circ', '\theta = 30\circ', '\theta = 60\circ', 'Interpreter', 'tex','Location', 'NorthWest')
set(gca, 'FontSize', 12)

%% FUNCTIONS
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

function [Cp, Ma_ee, p_ee, T_ee] = modified_newtonian_theory(th, Ma, gamma_air, T_inf, V_inf, rho_air, R_air, cp_air)
% This function calculates flow parameters using modified newtonian theory
%
% INPUTS:
%   th          :   (rad)   [npj x 1]       Angular spacing of the model
%   Ma          :   (-)                     Mach number of the freestream flow
%   gamma_air   :   (-)                     Ratio of specific heats for air
%   T_inf       :   (K)                     Temperature of the freestream flow
%   V_inf       :   (m/s)                   Velocity of the freestream flow
%   rho_air     :   (kg/m3)                 Density of the freestream flow
%   R_air       :   (J/kgK)                 Specific gas constant for air
%   cp_air      :   (J/kgk)                 Specific heat capacity of air at constant pressure
%
% OUTPUTS:
%   Cp      :   (-)     Local pressure coefficient at the surface
%   Ma_ee   :   (-)     Mach number at the edge of the boundary layer
%   p_ee    :   (Pa)    Pressure at the edge of the boundary layer
%   T_ee    :   (K)     Temperature at the edge of the boundary layer

    % CHECK METHODOLOGY FOR OBTAINING BOUNDARY LAYER EDGE VALUES
    angle = abs(th) - (pi/2);               % (rad) Local angle of each section with the freestream
    absolute_angle = abs(angle);            % (rad) MNT does not depend on sign and so we remove the negatives

    Cp_max =2/(gamma_air*Ma^2)*((((gamma_air+1)*Ma)^2/(4*gamma_air*Ma^2-2*(gamma_air-1)))^(gamma_air/(gamma_air-1))*((1-gamma_air+2*gamma_air*Ma^2)/(gamma_air+1))-1); %calculate maximum Cp using Rayleigh Pitot formula
    Cp = Cp_max .* (sin(absolute_angle).^2);    % MNT to find pressure coefficient

    % Flow properties
    T_0 = T_inf + ((V_inf^2)/(2*cp_air));   % (K) Stagnation temperature of the flow
    p_inf = rho_air * R_air * T_inf;        % (Pa) Freestream pressure
    p_0_inf = p_inf*(1+((gamma_air-1)/2)*Ma^2)^(gamma_air/(gamma_air-1));   % (Pa) Freestream stagnation pressure
    p_0_e = p_0_inf*(((gamma_air+1)*Ma^2)/(2+(gamma_air-1)*Ma^2))^(gamma_air/(gamma_air+1))*((gamma_air+1)/((2*gamma_air*Ma^2)-(gamma_air-1)))^(1/(gamma_air-1));   % (Pa) Stagnation pressure at the edge of the boundary layer
    p_ee = p_inf*(1 + ((Cp*gamma_air*Ma^2)/2));     % (Pa)  Static Pressure at the edge of the boundary layer

    Ma_ee = (((p_0_e./p_ee).^((gamma_air-1)/gamma_air)-1).*(2./(gamma_air-1))).^(1/2);  % (-) Mach number at the edge of the boundary layer
    T_ee = T_0./(1+(((gamma_air-1)/2).*(Ma_ee.^2)));                                    % (K) Temperature at the edge of the boundary layer
end

function T = initial_temp_distribution(T_in, ~, npi, npj)
% This function creates an initial temperature distribution
%
% INPUTS:
%   T_in    :   (K)     Angular spacing of the model
%   input   :   (-)     Input parameter to control initial temperature distribution
%   npi     :   (-)     Number of spaces in the radial direction
%   npj     :   (-)     Number of spaces in the angular or tangential direction 
%
% OUTPUTS:
%   T   :   (K)     [npi x npj] Temperature field

    % Assumes a uniform temperature distribution of temperature T_in
    T = ones(npi, npj) .* T_in;
end

function [qconv_history, normqe, normqrad, normqsub, deltar, qe, mdot] = ...
    heatIn(c, rmax, R_nap, rho_nap, nu_nap, h_nap, nu_air, k_air, V_inf, T, T_aw, Pr, npi, dt, dth)
% This function calculates the heat entering the surface of the model based on fluid and material parameters
%
% INPUTS:
%   c       :   (m-1)       Curvature of the model
%   rmax    :   (m)         Surface radius of the model
%   R_nap   :   (J/kgK)     Specific gas constant of the model material
%   rho_nap :   (kg/m3)     Density of the model material
%   nu_nap  :   (m2/s)      Kinematic viscosity of the model material
%   h_nap   :   (J/kg)      Sublimation enthalpy of the model material
%   nu_air  :   (m2/s)      Kinematic viscosity of air
%   k_air   :   (W/mK)      Thermal conductivity of air
%   V_inf   :   (m/s)       Freestream velocity of the air
%   T       :   (K)         % Temperature of the freestream
%   T_aw    :   (K)         % Adiabatic wall temperature
%   Pr      :   (-)         % Prandtl number
%   npi     :   (-)         % Radial grid amount
%   dt      :   (s)         % Time spacing
%   dth     :   (rad)       % Angular grid spacing
%
% OUTPUTS:
%   qconv_history   :   (W/m2)                  Mean convective heat transfer to the surface
%   normqe          :   (-)                     Mean normalised conductive heat transfer through the surface
%   normqrad        :   (-)                     Mean normalised radiative heat transfer to the surface
%   normqsub        :   (-)                     Mean normalised sublimative heat transfer to the surface
%   deltar          :   (m)         [npj x 1]   Recession rate of the surface
%   qe              :   (W/m2)      [npj x 1]   Conductive heat transfer to the surface
%   mdot            :   (kg/m2s)    [npj x 1]   Mass transfer at the surface

    % Radiation Heat Transfer Constants
    sigma = 5.67e-08;
    emiss = 0.9;

    D = 2./c;           % (m) Effective Diameter at each theta location

    % Nondimensional Parameters for Flow and Heat Transfer
    Re = V_inf * 2 * rmax / nu_air;         % (-) Reynolds Number of the flow
    Nu = 0.076 * (Re.^0.7).*(Pr.^0.37);     % (-) Flow Nusselt Number (Could be changed to local Nu number)
    h = Nu * k_air./D;                      % (W/m2K) Heat Transfer Coefficient based on curvature at theta
    Sc = 7./(T(npi-1,:)'.^0.185);           % (-) Schmidt Number
    Sh = Nu.*(Sc./Pr).^0.37;                % (-) Sherwood Number

    % Mass Transfer at each theta
    pw = exp(31.23252 - 8587.36./T(npi-1,:)');              % (Pa) Vapour Pressure of Naphthalene given by Sogin
    mdot = nu_nap.*Sh.*pw./(R_nap.*Sc.*D.*T(npi-1,:)');     % (kg/m2s) Mass flow of naphthalene away from the surface

    % Heat Transfer at each Theta
    qconv = h.* (T_aw - T(npi,:)');             % (W/m2) Convective heat transfer to surface
    qrad = emiss .* sigma .* T(npi-1,:)'.^4;    % (W/m2) Radiative heat transfer (IS T(npi-1,:) CORRECT?)
    qsub = mdot .* h_nap;                       % (W/m2) Sublimative heat transfer

    % Calculate Heat conducted INTO naphthalene
    qe = qconv - qrad - qsub;
    qconv_history = mean(qconv);                % (W/m2) Storing convetive heat flux
    normqe = mean(qe)/mean(qconv);              % (-) Normalised entering heat transfer
    normqrad = mean(qrad)/mean(qconv);          % (-) Normalised radiative heat transfer
    normqsub = mean(qsub)/mean(qconv);          % (-) Normalised sublimative heat transfer

    % Recession Rate
    deltar = -mdot.*dt./(rho_nap.*rmax.*dth);      % (m) Change in radius assuming depth is 1m)
end

function [aN, aE, aS, aW, aP, bP, sP] = Tcoefficiants(dth, dr, r, c, alpha_nap, dt, T, rmax, qe, k_nap, npi, npj)
    aN = zeros(npi-1, npj-1);
    aE = zeros(npi-1, npj-1);
    aS = zeros(npi-1, npj-1);
    aW = zeros(npi-1, npj-1);
    aP = zeros(npi-1, npj-1);
    sP = zeros(npi-1, npj-1);
    bP = zeros(npi-1, npj-1);
    
    for i=2:npi-1
        for j=2:npj-1
            %-------Specify coefficient matrices
            aE(i,j) = dth/dr(j)*(r(i,j)+dr(j)/2);
            aW(i,j) = dth/dr(j)*(r(i,j)-dr(j)/2);
            aN(i,j) = (dr(j+1)/(1+c(j+1)*r(i,j+1))+dr(j)/(1+c(j)*r(i,j)))/(2*dth*(1+c(j)*r(i,j)));
            aS(i,j) = (dr(j-1)/(1+c(j-1)*r(i,j-1))+dr(j)/(1+c(j)*r(i,j)))/(2*dth*(1+c(j)*r(i,j)));
            sP(i,j) = dr(j)*dth*r(i,j)/(alpha_nap*dt);
            aP(i,j) = aW(i,j)+aE(i,j)+aN(i,j)+aS(i,j)+sP(i,j);
            bP(i,j) = sP(i,j)*T(i,j);
        end
    end
    %---Eastern and western boundary conditions
    for j = 2:npj-1
        %-----Prescribe heat flux
        aE(npi-1,j)=0;
        aP(npi-1,j)=aW(npi-1,j)+aN(npi-1,j)+aS(npi-1,j)+sP(npi-1,j);
        bP(npi-1,j)=bP(npi-1,j)+rmax(j)*dth*qe(j)/k_nap;
        %-----No temperature gradient over the western face
        aW(2,j)=0;
        aP(2,j)=aE(2,j)+aN(2,j)+aS(2,j)+sP(2,j);
    end
    %---No temperture gradient over the northern and southern face
    for i = 2:npi-1
        %---- Southern boundary
        aN(i,npj-1)=0;
        aP(i,npj-1)=aS(i,npj-1)+aE(i,npj-1)+aW(i,npj-1)+sP(i,npj-1);
        %---- Northern boundary
        aS(i,2)=0;
        aP(i,2)=aN(i,2)+aE(i,2)+aW(i,2)+sP(i,2);
    end
end

function phi = gs_solver(phi, aN, aE, aS, aW, aP, bP, iend, jend)
    for i = 2:iend-1
        for j = 2:jend-1
            phi(i,j) = (aE(i,j)*phi(i+1,j)+aW(i,j)*phi(i-1,j)+aN(i,j)*phi(i,j+1)+aS(i,j)*phi(i,j-1)+bP(i,j))/aP(i,j);
        end
    end
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
    
        c(j) = sqrt((4*c1^2)/(c2^2+c3^2-4*c1*c4)); %Caalculate local curvature
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

function [X,Y] = grid_plotting(r, th, rmax, npi, npj, fig)
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
    axis equal
    title('Solution Grid');
    xlabel('Horizontal Axis (m)');
    ylabel('Vertical Axis (m)');
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

function [material, shape, temp, press, speed, npi, npj, Tend, dt] = userInput()
% This function gathers user input using a dialog window with various menu options
%
% DOCUMENTATION: https://blogs.mathworks.com/pick/2009/12/25/input-dialog-box-on-steroids/
    hfig = figure('Name', 'User Input', 'CloseRequestFcn', @close_req_fun, 'menu', 'none');

    % MODEL MATERIAL
    model_material = {'Water Ice', 'Dry Ice', 'Naphthalene', 'Camphor', 'Naphthalene (OLD)'};
    material = 'Naphthalene';

    % MODEL SHAPE
    model_geometry = {'Hemisphere', 'Ellipse'};
    shape = 'Hemisphere';

    % WIND TUNNEL PARAMETERS
    temp = '380';
    press = '100';
    speed = '5';

    % SOLUTION SPACE CREATION
    npi = 23;
    npj = 72;
    Tend = 10;
    dt = 0.001;

    % DIALOG BOX CREATION
    set(hfig,'menu','none')
    % MODEL PARAMETERS
    uicontrol('Style', 'Text', 'String', "Please select the model material:", ...
        'Parent', hfig, 'Units', 'normalized', ...
        'Position', [0.2, 0.9, 0.6, 0.1]);
    material_menu = uicontrol('Style', 'popupmenu', 'String', model_material, ...
        'Parent', hfig, 'Units', 'Normalized', ...
        'Position', [0.2, 0.87, 0.6, 0.1 ]);

    uicontrol('Style', 'Text', 'String', "Please select the model geometry:", ...
        'Parent', hfig, 'Units', 'normalized', ...
        'Position', [0.2, 0.8, 0.6, 0.1]);
    geometry_menu = uicontrol('Style', 'popupmenu', 'String', model_geometry, ...
        'Parent', hfig, 'Units', 'normalized', ...
        'Position', [0.2, 0.77, 0.6, 0.1]);

    % WIND TUNNEL PARAMETERS
    uicontrol('Style', 'Text', 'String', "Wind Tunnel Stagnation Temperature (K):", ...
        'Parent', hfig, 'Units', 'normalized', ...
        'Position', [0.2, 0.7, 0.6, 0.1]);
    temperature = uicontrol('Style', 'Edit', 'String', temp, ...
        'Parent',hfig,'Units','Normalized', ...
        'Position', [0.2, 0.71, 0.6, 0.05]);

    uicontrol('Style', 'Text', 'String', "Wind Tunnel Stagnation Pressure (kPa):", ...
        'Parent', hfig, 'Units', 'normalized', ...
        'Position', [0.2, 0.6, 0.6, 0.1]);
    pressure = uicontrol('Style', 'Edit', 'String', press, ...
        'Parent',hfig,'Units','Normalized', ...
        'Position', [0.2, 0.61, 0.6, 0.05]);

    uicontrol('Style', 'Text', 'String', "Wind Tunnel Mach Speed (-):", ...
        'Parent', hfig, 'Units', 'normalized', ...
        'Position', [0.2, 0.5, 0.6, 0.1]);
    speed_menu = uicontrol('Style', 'Edit', 'String', speed, ...
        'Parent',hfig,'Units','Normalized', ...
        'Position', [0.2, 0.51, 0.6, 0.05]);

    % SOLUTION SPACE CREATION
    uicontrol('Style', 'Text', 'String', "Radial Grid Length (-):", ...
        'Parent', hfig, 'Units', 'normalized', ...
        'Position', [0.2, 0.4, 0.6, 0.1]);
    npi_menu = uicontrol('Style', 'Edit', 'String', npi, ...
        'Parent',hfig,'Units','Normalized', ...
        'Position', [0.2, 0.41, 0.6, 0.05]);

    uicontrol('Style', 'Text', 'String', "Tangential Grid Length (-):", ...
        'Parent', hfig, 'Units', 'normalized', ...
        'Position', [0.2, 0.3, 0.6, 0.1]);
    npj_menu = uicontrol('Style', 'Edit', 'String', npj, ...
        'Parent',hfig,'Units','Normalized', ...
        'Position', [0.2, 0.31, 0.6, 0.05]);

    uicontrol('Style', 'Text', 'String', "Run Time (s):", ...
        'Parent', hfig, 'Units', 'normalized', ...
        'Position', [0.2, 0.2, 0.6, 0.1]);
    Tend_menu = uicontrol('Style', 'Edit', 'String', Tend, ...
        'Parent',hfig,'Units','Normalized', ...
        'Position', [0.2, 0.21, 0.6, 0.05]);

    uicontrol('Style', 'Text', 'String', "Time Step (s):", ...
        'Parent', hfig, 'Units', 'normalized', ...
        'Position', [0.2, 0.1, 0.6, 0.1]);
    dt_menu = uicontrol('Style', 'Edit', 'String', dt, ...
        'Parent',hfig,'Units','Normalized', ...
        'Position', [0.2, 0.11, 0.6, 0.05]);


    % BUTTON CONTROLS
    uicontrol('Style', 'pushbutton', 'String', 'Accept', ...
        'Parent',hfig,'Units','Normalized', ...
        'Position', [0.1, 0.05, 0.3, 0.05],...
        'Callback','close(gcbf)');
    cancel = uicontrol('Style', 'pushbutton', 'String', 'Cancel', ...
        'Parent',hfig,'Units','Normalized', ...
        'Position', [0.6, 0.05, 0.3, 0.05],...
        'Tag','0','Callback',@cancelfun);
    %wait for figure being closed (with OK button or window close)
    uiwait(hfig)
    %figure is now closing
    if strcmp(cancel.Tag,'0')%not canceled, get actual inputs
        material = model_material{material_menu.Value};
        shape = model_geometry{geometry_menu.Value};
        temp = temperature.String;
        press = pressure.String;
        speed = speed_menu.String;
        npi = npi_menu.String;
        npj = npj_menu.String;
        Tend = Tend_menu.String;
        dt = dt_menu.String;
    end
    %actually close the figure
    delete(hfig)
    
    function cancelfun(h,~)
        set(h,'Tag','1')
        uiresume
    end
    
    function close_req_fun(~,~)
        uiresume
    end
end