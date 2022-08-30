%**** Author: Laurence Maskell
%**** Objective: Calculation of temperature distribution and
% recession rates for a naphthalene heat shield
% model.
%
%**** References:
% [1] H.K. Versteeg & W. Malalasekera:
% An Introduction to Computational Fluid Dynamics: The Finite Volume Method,
% 2nd edition, Pearson / Prentice Hall, 2007.
% [2] S.V. Patankar: Numerical Heat Transfer and Fluid Flow, McGraw-Hill, 1980.
clear
clc
close all
set(0,'defaulttextinterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
tic
% Preallocate movie structure.
global vidObj
vidObj = VideoWriter('ablation.avi');
vidObj.FrameRate = 40;
open(vidObj);

%-Domain size in theta and r
global r
global th
global rmax histrmax

%-Grid spacing
global npi npj % number of grid points in r and theta
global RmaxEnd % Radius at start and minimum radius at end
global deltar c

%-Time discretization
global Tend dt
global nt

%-Solution matrices
global T Told % temperature fields
global aW aE aN aS aP % coefficient matrcies (relaxation parameters?)
global sP bP

%-Boundary conditions
global Tw Te Tn Ts % Dirichlet boundary temperatures

%-Thermodynamics
global alpha
global k_sg R_air rho_air nu_air V_inf k_air %k_sg Sutton-Graves parameter
global Pr R_nap h_nap rho_nap cp_nap k_nap nu_nap
global T_aw T_s

%-Time Histories
global histr histT histc histrdot histmdot
global timeIter mdot
global qconv qe qrad qsub rarcl qww
global normqe normqrad normqsub

%-Inputs
global Ma
global cp_air
global gamma
global T_inf
global Ma_ee V_ee T_ee

%-Solution parameters
eps=1e-6; % Tolerance
maxiter=20000;

%---Spatial discretization
npi = 23; % cells in r-direction
npj = round(pi*npi); % cells in theta-direction (scaling for aspect ratio close to unity)

%---Time discretization
Tend = 10; % end time
dt = 0.001; % time step
nt = round(Tend/dt); % number of time steps

%---Freestream windtunnel conditions
T_s = 330; % stagnation temperature [K] 380
Ma = 5; % mach number

%---Properties of air
gamma = 1.4; % ratio of specific heats of air []
R_air = 287; % gas constant for air [J/kgK]
cp_air = 1004; % specific heat (const p) (solid) [J/kgK]
rho_air = 0.29; % density [kg/m3] 0.29
nu_air = 1.568e-05; % kinematic viscosity [m2/s]
k_air = 0.025; % conductivity [W/mK]
k_sg = 1.7415e-04; % Sutton-Graves constant [kg0.5/m]
Pr = 0.7; % Prandtl no.

%--- Properties of Naphthalene
rho_nap = 1000; % density (solid) [kg/m3]
cp_nap = 1300; % specific heat (const p) (solid) [J/kgK]
k_nap = 0.34; % conductivity (solid) [W/mK]
alpha = k_nap / (rho_nap * cp_nap); % thermal diffusivity (solid) [m2/s]
h_nap = 5.57e+05; % sublimation enthalpy (J/kg)
nu_nap = 9.689e-07; % kinematic viscosity (vapour) [m2/s]
R_nap = 64.87; % gas constant [J/kgK]

%-Create domain
grid()

%--- Potential flow solutions temperature
%-Freestream
T_inf = T_s/(1+(Ma^2*gamma*R_air)/(2*cp_air)); % freestream temperature [K]
V_inf = Ma*((gamma*R_air*T_inf)^(1/2)); % freestream velocity [m/s]

cp()


rf=Pr^(1/3); %Recovery factor (turbulent flow recovery factor equation, Dr Bruce Ch 1 notes, pg 28)
T_aw = T_ee.*(1+rf*(gamma-1)*Ma_ee.^(2)/2); %pg 508 Viscous Fluid Flow Frank White

%---Initialize the temperture field and boundary conditions
Tw=293; %Changed from 293 > 300
Ts=293;
Te=T_aw;
Tn=293;

%---Set initial condition
T=293*ones(npi,npj);



%-Allocating remaining arrays:
aE=zeros(npi,npj); aW=zeros(npi,npj); aN=zeros(npi,npj);
aS=zeros(npi,npj); aP=zeros(npi,npj); sP=zeros(npi,npj);
bP=zeros(npi,npj); Told=zeros(npi,npj);
histr=zeros(nt,npi,npj); histT=zeros(nt,npi,npj);
histrdot = zeros(nt,3); histmdot = zeros(nt,3); histc = zeros(nt,1);
normqe=zeros(nt,1); normqrad=zeros(nt,1); normqsub=zeros(nt,1);


%-Implement boundary conditions
bound()


%-Time loop
for timeIter=1:nt

    %---Calculate the boundary condition qe
    heatin()

    %---Update coefficient matrices - (uses old time level for bP)
    Tcoeff()

    %---Gauss Seidel Solver for numerical system
    for iter=1:maxiter

        %------Store initial temperature distribution
        if (iter==1)
            Told = T;
        end

        T = solve(T,bP,1,1,npi,npj);

        %----- Cacluate 2-norm of current and old temperature field
        if (iter==1)
            norm0=norm2(Told,T);
            Told = T;
        else
            normEnd=norm2(Told,T);
            %---------Check for convergence and update old temperture field
            if (normEnd<eps*norm0)
                break
            else
                Told=T;
            end
        end
    end

    %---Save results to matrices
    histr(timeIter,:,:) = r;
    histrmax(timeIter,:) = rmax;
    histT(timeIter,:,:) = T;
    for i = 3:-1:1
        index = i*(npj/6);
        histrdot(timeIter,4-i) = deltar(index)/dt;
        histmdot(timeIter,4-i) = mdot(index)/dt;
    end
    histc(timeIter) = c(round(npj/2));


    %---Remesh in r-direction
    remesh()

    %---Check for reached minimum r-value
    if (max(rmax)<RmaxEnd)
        disp('Reached minimum radius, program aborted')
        break
    elseif (c(round(npj/2))< 10)
        disp('Reached minimum allowable curvature, program aborted')
        break
    end

    fprintf('Number of iterations: %i \t Time: %f\n',iter,timeIter*dt)
end


plotResults()

figure(4)
plot(normqe,'-r', 'LineWidth', 4)
hold on
plot(normqrad, '-k', 'LineWidth', 4)
hold on
plot(normqsub, '-b', 'LineWidth', 4)
hold off
legend('q_{cond}', 'q_{rad}', 'q_{sub}', 'Interpreter', 'tex', 'Location', 'East','FontSize', 24)

ylabel('$\stackrel{.}{q}$/$\stackrel{.}{q}_{conv}$', 'Interpreter', 'latex', 'FontSize', 24)
xlabel('t [s]', 'FontSize', 24)
set(gca, 'FontSize', 14)

figure(6)
plot(histmdot, 'LineWidth', 2)
ylabel('$\stackrel{.}{m}$ [kg/$m^2$s]', 'Interpreter', 'latex')
xlabel('t [s]')
legend('\theta = 0\circ', '\theta = 30\circ', '\theta = 60\circ', 'Interpreter', 'tex','Location', 'NorthWest')
set(gca, 'FontSize', 12)

toc

function grid()
%
%**** Purpose: Defining the initial grid
%
global r th rmax thmax c
global rmin thmin npi npj dr dth Rmax0 RmaxEnd
global qratio qstag qconv Re Nu Sc Sh h Tfluid pw mdot deltar qe


%---Discretization Allocation
r = zeros(npi,npj); th = zeros(npj,1); % Spatial location
rmax = zeros(npj,1); % grid spacing in r

%---Thermo allocations
qratio=zeros(npj,1); qstag=zeros(npj,1); qconv=zeros(npj,1);
Re=zeros(npj,1); Nu=zeros(npj,1); Sc=zeros(npj,1); Sh=zeros(npj,1); h=zeros(npj,1);
Tfluid=zeros(npj,1); pw=zeros(npj,1); mdot=zeros(npj,1);
qe=zeros(npj,1); deltar=zeros(npj,1);

%---Length of calculation domain in r and theta direction
rmax(:) = 0.025;
Rmax0 = max(rmax);
RmaxEnd = 0.5*max(rmax);
rmin = 0.1*rmax(:);
thmax = 0.5*pi;
thmin = -0.5*pi;
c = 1./rmax;

%---Grid spacing in r and theta
dr = (rmax-rmin)./(npi-2);
dth = (thmax-thmin)./(npj-2);

%---Length variable for the scalar points in the r and theta direction
r(1,:) = rmin;
r(2,:) = rmin+0.5*dr;
th(1) = thmin;
th(2) = thmin+0.5*dth;

for i=3:npi-1
    r(i,:)=r(i-1,:)+dr';
end
r(npi,:)=r(npi-1,:)+0.5*dr';
for j=3:npj-1
    th(j)=th(j-1)+dth;
end
th(npj)=th(npj-1)+0.5*dth;
end

function cp()
%
% Find local angles that the flow is being turned by
% so Cp can be calculated at each point at the surface.
global rmax deltar npi npj th gamma Ma T_inf V_inf cp_air rho_air R_air X Y Thgrid histr timeIter p_e T_s
%global T_ee p_inf
global Ma_ee T_ee V_ee R_qc p_0_1 p_w T_0_e
global qww T_aw Pr T
rmax = rmax + deltar;




x = rmax.*cos(th); %Need to change this value to max X - from plot results section
y = rmax.*sin(th); %Need to change this value to max Y -

%x_cp = X(npi, :)';
%y_cp = Y(npi, :)';

%R_qc = (x_cp.^2 + y_cp.^2).^(1/2);

% coord = [x, y];
ang = abs(th) - (pi/2); %thetaval
absang = abs(ang);
%{
h = diff(th(1:2)); %thetaval
dx = -gradient(x_cp, h);
dy = -gradient(y_cp, h);
slope = (-atan2d(dy, dx) -(0)); %flipud flips the matrix upside down, so the angle system is above stag point is +ve
slope(slope>(90)) = slope(slope>(90)) - 180;
slope = abs(slope);

figure(3)
plot(x_cp, y_cp, '-+')
hold on
%quiver(x_cp, y_cp, dx, dy, '-r');
hold off
axis('equal')
angstr = sprintfc('\\angle %.0f \\circ', slope);
text(x_cp, y_cp, angstr, 'HorizontalAlignment', 'left', 'Interpreter', 'tex')
xlabel('x [m]', 'FontSize',16)
ylabel('y [m]', 'FontSize', 16)
%}

Cp_max = 2/(gamma*Ma^2)*((((gamma+1)*Ma)^2/(4*gamma*Ma^2-2*(gamma-1)))^(gamma/(gamma-1))*((1-gamma+2*gamma*Ma^2)/(gamma+1))-1); %calculate maximum Cp using Rayleigh Pitot formula
%Cp = Cp_max.*sind(slope).^2;
Cp = Cp_max.*sin(absang).^2;

T_0 = T_inf + ((V_inf^2)/(2*cp_air));
p_inf = rho_air*R_air*T_inf; %Perfect gas law
p_0_inf = p_inf*(1+((gamma-1)/2)*Ma^2)^(gamma/(gamma-1));
p_0_1 = p_0_inf*(((gamma+1)*Ma^2)/(2+(gamma-1)*Ma^2))^(gamma/(gamma+1))*((gamma+1)/((2*gamma*Ma^2)-(gamma-1)))^(1/(gamma-1));

p_e = p_inf*(1 + ((Cp*gamma*Ma^2)/2));

p_w = Cp.*(0.5.*rho_air.*V_inf.^2) + p_inf;

Ma_ee = (((p_0_1./p_e).^((gamma-1)/gamma)-1).*(2./(gamma-1))).^(1/2);
T_ee = T_0./(1+(((gamma-1)/2).*(Ma_ee.^2)));

T_0_e = (1 + ((gamma-1)/2).*Ma_ee.^2).*T_ee;

V_ee = Ma_ee.*(gamma*R_air.*T_ee).^(1/2);
%{
rho_e = p_e./(R_air.*T_0_e); % Density at edge of BL (removed T_ee as we require stagnation dege density)

mu_e = 1.7894e-5.*(T_0_e./273.11).^(3/2).*((273.11 + 110.56)./(T_0_e + 110.56));

mu_w = 1.7894e-5.*(T(npi-1,:)'./273.11).^(3/2).*((273.11 + 110.56)./(T(npi-1,:)' + 110.56));

rho_w = p_w./(R_air.*303); % Find density at wall

duedx = (1./R_qc).*sqrt((2.*(p_e - p_inf))./rho_e); %Need to find Stagnation

h_e = cp_air.*T_0_e;

h_w = cp_air.*303; %Sublimation temp of Naphthalene

qww = (0.763*Pr.^(-0.6).*((rho_w.*mu_w).^(0.1)).*((rho_e.*mu_e).^(0.4)).*(h_e*h_w).*((duedx).^(1/2)));

figure(5)
rarcl = (X(npi, :)')./(R_qc);
plot(rarcl, Ma_ee, ':bs');
xlabel('x/R_{N}', 'Interpreter', 'tex', 'FontSize', 16)
ylabel('Ma_{e}', 'Interpreter', 'tex', 'FontSize', 16)
set(gca, 'FontSize', 14)
hold on
expy = [2.93 2.67 2.228 1.76 1.5 1.223 1.0 0.80 0.62 0.386 0.23 0.026];
expx = [0 0.112 0.24 0.37 0.50 0.61 0.68 0.80 0.88 0.94 0.98 1];
plot(expx, expy, '-r*')
legend('Numerical Model Results', 'Winkler et al Theoretical Results', 'FontSize', 16);
%}
end



function remesh()
%
%**** Purpose: Remesh the grid in r-direction
%
global r th rmax c rmin npi npj dr deltar

%---Change local radius and curvatures
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
    % https://uk.mathworks.com/matlabcentral/answers/284245-matlab-code-for-computingcurvature-equation
end
c(1)=c(2); % First two
c(npj)=c(npj-1); % & Last two are the same value

%---Grid spacing in r-direction
dr = (rmax-rmin)/(npi-2); %Specify new spacing


%---Discretization in r
r(1,:) = rmin; % Define all rows of row 1 to be rmin, which = 0.0025.
r(2,:) = rmin+0.5*dr; % Define all rows of row 2 to be rmin + 0.5dr
for i=3:npi-1
    r(i,:)=r(i-1,:)+dr';
end
r(npi,:)=r(npi-1,:)+0.5*dr';
end



function bound()
%
%**** Purpose: Prescribe Dirichlet BCs (not in effect if von Neumann is used in Tcoeff())
%
global npi npj T Tw Te Tn Ts

%---West and East BC
for j = 2:npj-1
    T(1,j) = Tw;
    T(npi,j) = Te(j);
end

%---North and South BC
for i = 2:npi-1
    T(i,1)= Ts;
    T(i,npj) = Tn;
end
end

function heatin()
%
%**** Purpose: caluculate the input heat to the surface at each theta location
%
global qconv qrad qsub qe
global normqe normqrad normqsub
global rmax c dth dt Re Nu Sc Sh h
global pw mdot deltar V_inf nu_air k_air Pr R_nap cp_air
global nu_nap h_nap rho_nap T npi
global timeIter
global Rgrid histr p_e p_inf T_ee R_air npi R_qc T_aw
global p_w T_0_e qww V_ee T_inf

sigma = 5.67e-08;
emiss = 0.9;

% for i = timeIter:timeIter
% Rgrid = squeeze(histr(i,:,:));
% end

% Effective diameter at each theta location
D = 2./c;


%---Find convective heat transfer coefficiient
Re = V_inf*2*rmax/nu_air; % Reynolds no.
Nu = 0.076*(Re.^0.7)*(Pr^0.37); % Average Nusselt no. (CHANGE THIS TO LOCAL VALUE)
h = Nu*k_air./D; % Average coefficient based on curvature at theta

%h = qww./(T(npi-1,:)'-T_inf);
%Nu = (h.*D)./k_air;



%---Find convective heat transfer at each theta location
qconv = h.*(T(npi,:)' - T(npi-1,:)');

%---Calculate non-dimensional numbers
Sc = 7./(T(npi-1,:)'.^0.185); % Schmidt no
Sh = Nu.*(Sc./Pr).^0.37; % Sherwood no

%---Get mass transfer rate at each theta
pw = exp(31.23252 - 8587.36./T(npi-1,:)'); % Naphthalene vapour pressur (IMPROVE?)
mdot = nu_nap*Sh.*pw./(R_nap*Sc.*D.*T(npi-1,:)'); % Mass transfer rate (per unit area)

qrad = emiss*sigma*T(npi-1,:)'.^4;
qsub = mdot*h_nap;

%---Finally calculate heat conducted into the Naphthalene
qe = qconv - qrad - qsub; %Replaced qconv with qww
normqe(timeIter) = mean(qe)/mean(qconv); %Replaced qconv with qww
normqrad(timeIter) = mean(qrad)/mean(qconv); %Replaced qconv with qww
normqsub(timeIter) = mean(qsub)/mean(qconv); %Replaced qconv with qww


%---Also compute change in radius of each theta location (assume depth 1m)
deltar = -mdot*dt./(rho_nap*rmax*dth);
end



function Tcoeff()
%
%**** Purpose: To calculate the coefficients for the T equation.
% Page 72 - 73 Numerical Heat Transfer and fluid flow
% Unsteady heat equation in a coordinate system for plane flows
global r c npi npj dr dth dt T aW aE aN aS aP sP bP qe alpha k_nap rmax
for i=2:npi-1
    for j=2:npj-1
        %-------Specify coefficient matrices
        aE(i,j) = dth/dr(j)*(r(i,j)+dr(j)/2);
        aW(i,j) = dth/dr(j)*(r(i,j)-dr(j)/2);
        aN(i,j) = (dr(j+1)/(1+c(j+1)*r(i,j+1))+dr(j)/(1+c(j)*r(i,j)))/(2*dth*(1+c(j)*r(i,j)));
        aS(i,j) = (dr(j-1)/(1+c(j-1)*r(i,j+1))+dr(j)/(1+c(j)*r(i,j)))/(2*dth*(1+c(j)*r(i,j)));
        sP(i,j) = dr(j)*dth*r(i,j)/(alpha*dt);
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
function [phi] = solve(phi,b,istart,jstart,iend,jend)
%
%**** Purpose: To solve the algebraic equation with Gauss-Seidel solver
%
global aE aW aN aS aP
%---Solving from left to right
for i=istart+1:iend-1
    %------Solving from down to up
    for j=jstart+1:jend-1
        phi(i,j)=(aE(i,j)*phi(i+1,j)+aW(i,j)*phi(i+1,j)+aN(i,j)*phi(i,j+1)+aS(i,j)*phi(i,j-1)+b(i,j))/aP(i,j);
    end
end
end
function [norm] = norm2(T11,T22)
%
%**** Purpose: To calculate the 2-norm of the difference of two matrices
%
global dr dth
nj = size(T11,2);
ni = size(T11,1);
norm=0;
%---Calculate 2-norm
for i=2:ni-1
    for j=2:nj-1
        norm=norm+(T11(i,j)-T22(i,j))^2;
    end
end
norm=sqrt(dr(floor(ni/2))*dth*norm);
end
function plotResults()
%
%**** Purpose: Plotting time history
%
global npi npj th timeIter histr histT vidObj T_s dt X Y Thgrid Rgrid
histT(:,npi,:) = histT(:,npi-1,:);
Thgrid = zeros(npi,npj);
for i = 1:npi
    Thgrid(i,:) = th';
end
figure(1)
for i = timeIter:timeIter
    Rgrid = squeeze(histr(i,:,:));
    X = Rgrid.*cos(Thgrid);
    Y = Rgrid.*sin(Thgrid);
    clf
    surface(X,Y,squeeze(histT(i,:,:)))
    hold on
    set(gca,'FontSize',12)
    xlim([0 0.026])
    ylim([-0.026 0.026])
    view(2)
    daspect([1 1 1])
    xlabel('x [m]','fontsize',12)
    ylabel('y [m]','fontsize',12)
    shading interp
    lighting phong
    colorbar
    caxis([290 400])
    text(0.02,0.021,'T [K]','FontSize',12)
    text(0.02,0.024,sprintf('t = %.2fs',i*dt))
    pause(0.00001)
    mov = getframe(gcf); %add current figure as frame to mov
    writeVideo(vidObj,mov); %write frame to vidObj
end
% Close the file.
close(vidObj);
end