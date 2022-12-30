% Heat and mass transfer in a laminar flow
% Based on 'Analysis Of Heat And Mass Transfer' by Eckert, 1986
% 
% Aerosol growth by condensation
% Based on 'Atmospheric Chemistry and Physics: From Air Pollution to Climate
% Change' by Seinfeld & Pandis, 2016
%
% Created by Egor Demidov (ed242@njit.edu)

%% Model parameters (SI units)

clear all;
clc;
close all;

global rg;
global mw;
global di;

% Length of the growth tube
l = 1; % m
% Radius of the growth tube
w = 2.3e-3; % m
% Volumetric flow rate
q = 0.3 * 1e-3 / 60; % m^3/s
% Mean flow velocity
v_mean = q / (pi * w^2); % m/s
% Thermal diffusivity
alpha = 2.17e-5; % m^2/s
% Mass diffusivity
di = 2.42e-5; % m^2/s
% Antoine equation with appropriate coefficients (Pa, K)
antoine = @(t) 1e5.*10.^(6.20963-2354.731./(t+7.559));
% Universal gas constant
rg = 8.314; % J/(mol*K)
% Molar mass of air
mw_air = 28.97e-3; % kg/mol
% Density of coating material
rho = 997; % kg/mol
% Molar mass of coating material
mw = 18.015e-3; % kg/mol
% Surface tension of coating material
gamma = 72e-3; % N/m
% Name of file with temperature profile
filename = 'temp_profile.csv';
% Initial particle diameter
dp0 = 240e-9;  % m
% Integration linearization density
d_int = 1e3; % 1/m
% Number of thin shells for averaging
n_shells = 150;
% Plot experimental data
plot_experimental = true;
% Name of file with experimental data
exp_filename = 'experimental_data.csv';
% Export growth data
export_data = false;
% Export growth file name
out_filename = 'growth.csv';
% Initial condition for heat trasfer
init_heat = 25 + 273.15; % K
% Initial condition for mass transfer
init_mass = antoine(init_heat)/(rg*init_heat); % mol/m^3
% Residual tolerance
tolerance = 1e-6;

% Mesh size factor (lower value - higher computation time and accuracy)
mesh_s = 0.00075;
% WARNING: An appropriate mesh size needs to be picked for each problem
% individually. A mesh too coarse may lead to an incorrect solution, while
% a mesh too fine can cause the solver to throw an error. The optimal size
% for this particular problem was found to be 0.00075

%% Load temperature profile (Signal processing toolbox required)
% Load data from file
data = readtable(filename,'NumHeaderLines',1);
% Axial positions
z_data = transpose(table2array(data(:,1)))./100; % m
% Temperatures
t_data = transpose(table2array(data(:,2)))+273.15; % K
% Interpolate data with a cubic spline
t_spl = spline(z_data,[0 t_data 0]);

%% Nonconstant coefficient/boundary condition handles
% Diffusion coefficient
c = @(location, state) location.x;
% Source term for heat transfer
f_h = @(location, state) -state.uy(1,:).*location.x.*2.*v_mean.*(1-(location.x/w).^2)./alpha;
% Source term for mass tranfer
f_m = @(location, state) -state.uy(1,:).*location.x.*2.*v_mean.*(1-(location.x/w).^2)./di;
% Boundary condition for heat transfer
bc_h = @(location,state) ppval(t_spl,location.y);
% Boundary condition for mass transfer
bc_m = @(location,state) antoine(ppval(t_spl,location.y))./(rg.*ppval(t_spl,location.y));

%% Geometry initialization
% Construct a rectangular domain (x=0 is the centerline, axial symmetry
% assumed)
R1 = [3,4,0,w,w,0,l,l,0,0]';
gm = [R1];
sf = 'R1';
ns = char('R1');
ns = ns';
g = decsg(gm,sf,ns);

%% Heat transfer solution (PDE toolbox required)
% Initialize the PDE interface
heat_model = createpde();
geometryFromEdges(heat_model,g);

% Specify the initial condition
setInitialConditions(heat_model,init_heat);
% Set solver residual tolerance
mass_model.SolverOptions.ResidualTolerance = tolerance;

% Specify the coefficients
specifyCoefficients(heat_model,'m',0,'d',0,'c',c,'a',0,'f',f_h);

% Specify the boundary conditions
applyBoundaryCondition(heat_model,'Dirichlet','Edge',[3,2],'u',bc_h); % Inlet and wall temperature
applyBoundaryCondition(heat_model,'Neumann','Edge',[1,4],'q',0,'g',0); % Zero flux at other boundaries

% Build a finite element mesh
fprintf('Generating mesh for heat transfer\n')
generateMesh(heat_model,'Hmax',mesh_s);

% Solve the stationary system
fprintf('Solving for temperature\n')
heat_results = solvepde(heat_model);

% Plot the results
figure;
pdeplot(heat_model,'XYData',heat_results.NodalSolution-273.15,'ColorMap','hot');
pbaspect([1 3 1]);
title('Temperature, °C');

%% Mass transfer solution (PDE toolbox required)
% Initialize the PDE interface
mass_model = createpde();
geometryFromEdges(mass_model,g);

% Specify the initial condition
setInitialConditions(mass_model,init_mass);
% Set solver residual tolerance
mass_model.SolverOptions.ResidualTolerance = tolerance;

% Specify the coefficients
specifyCoefficients(mass_model,'m',0,'d',0,'c',c,'a',0,'f',f_m);

% Specify the boundary conditions
applyBoundaryCondition(mass_model,'Dirichlet','Edge',[3,2],'u',bc_m); % Inlet and wall concentration
applyBoundaryCondition(mass_model,'Neumann','Edge',[1,4],'q',0,'g',0); % Zero flux at other boundaries

% Build a finite element mesh
fprintf('Generating mesh for mass transfer\n')
generateMesh(mass_model,'Hmax',mesh_s);

% Solve the stationary system
fprintf('Solving for concentration\n')
mass_results = solvepde(mass_model);

% Plot the results
figure;
pdeplot(mass_model,'XYData',mass_results.NodalSolution,'ColorMap','jet');
pbaspect([1 3 1]);
title('Concentration, mol/m^3');

%% Supersaturation surface plot
figure;
pdeplot(mass_model,'XYData',mass_results.NodalSolution.*heat_results.NodalSolution.*rg./antoine(heat_results.NodalSolution)-1,'ColorMap','jet');
%title('Supersaturation');
pbaspect([1 3 1]);
xlabel('r, m');
ylabel('z, m');

%% Cut line plots for the paper
ys = linspace(0,l,l*d_int);
xs = zeros(1,l*d_int);
ptsm = interpolateSolution(mass_results,xs,ys);
ptsh = interpolateSolution(heat_results,xs,ys);

figure;
plot(pi.*alpha.*ys./q, ptsm.*rg.*ptsh./antoine(ptsh)-1,'-k','LineWidth',1.5);
hold on;
xs = 0.5*w.*ones(1,l*d_int);
ptsm = interpolateSolution(mass_results,xs,ys);
ptsh = interpolateSolution(heat_results,xs,ys);
plot(pi.*alpha.*ys./q, ptsm.*rg.*ptsh./antoine(ptsh)-1,'--k','LineWidth',1.5);
xs = 0.8*w.*ones(1,l*d_int);
ptsm = interpolateSolution(mass_results,xs,ys);
ptsh = interpolateSolution(heat_results,xs,ys);
plot(pi.*alpha.*ys./q, ptsm.*rg.*ptsh./antoine(ptsh)-1,':k','LineWidth',1.5);
legend('Centerline','r/R=0.5','r/R=0.8');
xlim([0 6]);
ylabel('Supersaturation');
xlabel('πα_tz/Q')

%% Growth calculation
% Evaluation points
ys = linspace(0,l,l*d_int);
solutions = struct();

% Solve the growth ODE at each thin shell
fprintf('Calculating growth at each thin shell\n');
for pos=(w/(n_shells-1)):((w-2*w/(n_shells-1))/(n_shells-1)):(w-w/(n_shells-1))
    xs = pos.*ones(1,l*d_int);
    
    % Evaluate concentration and temperature at each point
    ptsm = interpolateSolution(mass_results,xs,ys);
    ptsh = interpolateSolution(heat_results,xs,ys);
    
    % Solution buffers
    zs = [];
    rs = [];
    ss = [];
    
    % Calculate flow velocity
    v_z = 2*v_mean*(1-(pos/w)^2);

    % Solve the ODE
    rp = dp0 / 2;
    for i = 1:(length(ptsm)-1)
        zs(end+1)=ys(i);
        rs(end+1)=rp;
        ss(end+1)=ptsm(i)*rg*ptsh(i)/antoine(ptsh(i))-1;
        del_r = (ptsm(i)-antoine(ptsh(i))*exp(2*gamma*mw/(rho*rp*rg*ptsh(i)))/(rg*ptsh(i)))*trans_corr(rp, mfp(ptsh(i)))*di*mw/(rp*rho);
        del_t = (ys(i+1)-ys(i))/v_z;
        if del_r < 0 && rp <= dp0/2
            del_r = 0;
        end
        if rp + del_r * del_t < dp0/2
            del_r = (dp0/2 - rp) / del_t;
        end
        rp = rp + del_r*del_t;
    end
    
    % Append the buffers to a global structure
    solutions(end+1).z = zs;
    solutions(end).r = rs;
    solutions(end).ss = ss;
end

%% Weighted average calculation
% Calculate weighted average droplet radius and supersaturation at each
% axial position
fprintf('Calculating weighted average\n');
r_mean = [];
ss_mean = [];

for i = 1:length(solutions(2).z)
    r_buff = 0;
    ss_buff = 0;
    for n = 2:(n_shells+1)
        n_alg = n-1;
        r_buff = r_buff + solutions(n).r(i)*(2*n_alg-1);
        ss_buff = ss_buff + solutions(n).ss(i)*(2*n_alg-1);
    end
    r_mean(end+1) = r_buff/n_shells^2;
    ss_mean(end+1) = ss_buff/n_shells^2;
end

%% Growth plot and export
% Load experimental data if needed
if plot_experimental
    data = readtable(exp_filename,'NumHeaderLines',1);
    exp_data_z = transpose(table2array(data(:,1))); % cm
    exp_data_d = transpose(table2array(data(:,2))); % nm
end

figure;
yyaxis left;

% Plot experimental data if needed
if plot_experimental
    scatter(exp_data_z,exp_data_d,'sr','filled');
end

hold on;
plot(solutions(2).z*100,r_mean*2e9,'-');
ylabel('Diameter, nm');
xlabel('Axial position (z), cm');
yyaxis right;
plot(solutions(2).z*100,ss_mean)
ylabel('Supersaturation');

% Add legends
if plot_experimental
    legend('Experimental','Model','Supersaturation');
else
    legend('Model','Supersaturation');
end

hold off;

% Export solution if needed
if export_data
    out_table = table(transpose(solutions(2).z*100),transpose(r_mean*2e9),transpose(ss_mean));
    out_table.Properties.VariableNames = {'z, cm','D, nm','ss'};
    writetable(out_table,out_filename);
end

%% Growth function definitions
% Transition correction factor (Fuchs & Sutugin)
function out = trans_corr(rp,lambda)
    mc = 1; % Mass accomodation coefficient
    kn = lambda/rp; % Knudsen number
    out = (1+kn)/(1+0.3773*kn+1.33*kn*(1+kn)/mc);
end

% Mean free path, m
function out = mfp(t)
    global rg;
    global mw;
    global di;
    cbar = (8.0*rg*t/(pi*mw))^0.5; % Average molecular velocity
    out = 3*di/cbar;
end
