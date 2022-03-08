clear
% Electrode voltage
V = 10e4;
% Beam Conditions
charge = -1.602e-19;% Particle charge in C
mass = 9.109e-31;   % Particle mass in kg
numRays = 100; % Number of particles to launch
beamEnergy = 2e5;% eV of accelerating voltage
c = 2.99e8;     % Speed of light in SI

% Aperture and image planes
apPlace = 15;  % Position of the aperture
apSize = .5;    % Diameter of the aperture
imPlane = 20;    % Position of the image plane (terminal plane for particles)

% Calculate relativistic beam velocity
% Electron at 100 eV travels at approx 0.02c
ev2J = 1.6022e-19; % electron volt to joule conv factor
beJ = beamEnergy*ev2J; % beam energy in Joules
v0beam = 2.99e8*(1-(beJ/(mass*c^2)+1)^-2)^(1/2); % Beam velocity in m/s

% Define scale of geometry (bounds of model)
xmin = 0;
xmax = imPlane+1;
ymin = 0;
ymax = 20;

% Calculate electrode placement relative to outer bounds
uAp = 0.55*ymax;    % Upper edge of apertures
lAp = 0.45*ymax;    % Lower edge of apertures
L3 = 1/2*xmax;      % First edge of second element
L4 = L3+.1;         % Second edge of second element
L1 = L3-1;      % First edge of first element
L2 = L1+.1;    % Second edge of first element
L5 = L3+1;          % First edge of third element
L6 = L5+.1;         % Second edge of third element
emagmodel = createpde();

% Generate polygons which define element geometry
R2 = [3,4,L1,L2,L2,L1,uAp,uAp,ymax,ymax]';
R3 = [3,4,L1,L2,L2,L1, ymin, ymin,lAp,lAp]';
R5 = [3,4,L3,L4,L4,L3,uAp,uAp,ymax,ymax]';
R6 = [3,4,L3,L4,L4,L3, ymin, ymin,lAp,lAp]';
R7 = [3,4,L5,L6,L6,L5,uAp,uAp,ymax,ymax]';
R8 = [3,4,L5,L6,L6,L5, ymin, ymin,lAp,lAp]';
R4 = [3,4,xmin,xmax,xmax,xmin,ymax,ymax,ymin,ymin]';
gm = [R2,R3,R4,R5,R6,R7,R8];

% Create logical definition of regions and create geometry
ns = char('R2','R3','R4','R5','R6','R7','R8');
sf = 'R4-(R2+R3+R5+R6+R7+R8)';
ns = ns';
g = decsg(gm,sf,ns);
geom = geometryFromEdges(emagmodel,g);

% Apply static potentials to model edges (dirichlet)
applyBoundaryCondition(emagmodel,'dirichlet','Edge', [1:emagmodel.Geometry.NumEdges],'u',0);
applyBoundaryCondition(emagmodel,'dirichlet','Edge',[7,9,8,10,12,25],'u',V);
% Tell pdesolve the form of the equation to be solved.
% Solving elliptic Poisson eqn
specifyCoefficients(emagmodel,'m',0,'d',0,'c',1,'a',0,'f',0);
generateMesh(emagmodel,'Hmax',.1); % Generate Mesh over model geometry
% Results are to be used in the ray tracing
results = solvepde(emagmodel);
u = results.NodalSolution;

% Plot iso lines for voltage contours
figure('Name','Calculated Ray Trajectory Sample')
pdeplot(emagmodel,'XYData', u,'XYStyle','off','Contour','on','title',"Planar Einzel Lens")
hold on;
xlim([xmin, xmax]);
ylim([0.4*ymax,0.6*ymax]);
dOut = cell(numRays,1); % Init cell array to store the ODE solutions of arbitrary length

% Add parameters to the ODE solver
options = odeset('Events',@(t,y) myEventsFcn(t,y,xmax,ymax,apPlace,apSize,imPlane),'RelTol',1e-11);
% Generate initial conditions and trace particles through lens
% Run in parallel using parfor (300% faster than single thread runtime for large N)
% The loop can be changed betweeen parfor (parallel) and for (sequential)
parfor i = 1:numRays
    y0 = (i+19*numRays/2-1)*ymax/(20*numRays); % Create y spacing between launched particles
    vx0 = normrnd(v0beam,1e4); % Gaussian chromaticism in vx
    vy0 = normrnd(0,1e4); % Gaussian drift velocity in y direction
    [t,y] = ode45(@(t,y) trace(t,y,results,mass,charge),[0,1e-1],[xmax*0.05;y0;vx0;vy0],options); % Trace ray by solving ODE (RK45)
    dOut{i}(1,:) = t;
    dOut{i}(2,:) = y(:,1);
    dOut{i}(3,:) = y(:,2);
    dOut{i}(4,:) = y(:,3);
    dOut{i}(5,:) = y(:,4);
end

% Plot n sample rays from the calculated rays.
n = numRays;
for k = 1:n
    k = floor(k*numRays/n);
    plot(dOut{k}(2,:),dOut{k}(3,:))
end
% Plot n curves of energy vs. position on beam axis
figure('Name', 'Particle Energy (eV)')
hold on
xlabel('Position on Beam Axis (m)')
for k = 1:n
    L = floor(k*numRays/n);
    energy = 0.5*mass*(dOut{L}(4,:).^2+dOut{L}(5,:).^2)/ev2J;
    plot(dOut{L}(2,:),energy)
end

% Get beam density and plot at source and at image plane
figure('Name','Source & Image Current Density, N = ' + string(numRays))
hold on
xlabel('y position (m)')

% Grab init and final y values from traces
y0s = zeros(numRays,1);
yfs = zeros(numRays,1);
for j = 1:numRays
    y0s(j) = dOut{j}(3,1);
    if dOut{j}(2,end) > apPlace+0.1 % assign yfs if terminal is > aperture position
        yfs(j) = dOut{j}(3,end);    % Otherwise yfs = 0
    end
end
yfs = yfs(yfs>0);   % Reduce to only rays which pass the aperture position
% Show density of particles on source and image plane
[f0, xi] = ksdensity(y0s,'Bandwidth', 0.02);
[ff, xf] = ksdensity(yfs,'Bandwidth', 0.02);
plot(xi-10,f0,xf-10,ff)
legend({'Source Density', 'Image Density'},'Location','northwest')
%--------------------------------------------------------------------------
% Functions for particle tracing ODE:
% Define particle eqns of motion
function dydt = trace(~,y, res,m,q)
qp = [y(1);y(2)]; % Current position vector
[Ex, Ey] = evaluateGradient(res,qp); % Interpolate electric field at qp
dydt = [y(3); y(4); q/m*Ex; q/m*Ey]; % Linear system of eqns for Lorentz force
end
% Create event which stops ODE solver
% If particles reach the image plane or hit the aperture, terminate ODE
function [value,isterminal,direction] = myEventsFcn(~,y,xmax,ymax,apPlace,apSize,imPlane)
value = 1;
if y(1)>imPlane % check if particle reaches image plane
    value = 0;
end
% Check if particle hits the aperture
if y(1) > apPlace && y(1) < apPlace+0.1 && abs(y(2)-ymax/2) > apSize/2
    value = 0;
end
isterminal = 1;
direction = 0;
end