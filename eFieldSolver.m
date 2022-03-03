clear
% Electrode voltage
V = 13;
% Beam Conditions
charge = -1.602e-19;
mass = 9.109e-31;
numRays = 200; % Number of particles to launch 
beamEnergy = 10; % eV of accelerating voltage
c = 2.99e8; %SOL in SI
% Calculate relativistic beam velocity
ev2J = 1.6022e-19;
beJ = beamEnergy*ev2J;
v0beam = 2.99e8*(1-(beJ/(mass*c^2)+1)^-2)^(1/2);
% Define scale of geometry
xmin = 0;
xmax = 8;
ymin = 0;
ymax = 20;
% Calculate electrode placement
uAp = 0.55*ymax;
lAp = 0.45*ymax;
L1 = 3/8*xmax;
L2 = 3.1/8*xmax;
L3 = 1/2*xmax;
L4 = 1.025/2*xmax;
L5 = 5/8*xmax;
L6 = 5.1/8*xmax;
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
% Results are globally defined to be used in the ray tracing
results = solvepde(emagmodel);
u = results.NodalSolution;
% Plot iso lines for voltage contours
figure(1)
pdeplot(emagmodel,'XYData', u,'XYStyle','off','Contour','on','title',"Planar Einzel Lens")
options = odeset('Events',@(t,y) myEventsFcn(t,y,xmax),'RelTol',1e-11);
hold on;
xlim([xmin, xmax]);
ylim([0.4*ymax,0.6*ymax]);
dOut = cell(numRays,1); % Init cell array to store the ODE solutions of arbitrary length

% Generate initial conditions and trace particles through lens
% Run in parallel using parfor (300% faster than single thread runtime for large N)
parfor i = 1:numRays
    y0 = (i+19*numRays/2-1)*ymax/(20*numRays); % Create y spacing between particles
    vx0 = normrnd(v0beam,1e4); % Gaussian chromaticism in vx
    vy0 = normrnd(0,1e4); % Gaussian chromaticism in vy
    [t,y] = ode45(@(t,y) trace(t,y,results,mass,charge),[0,1e-1],[xmax*0.05;y0;vx0;vy0],options); % Trace ray by solving ODE (RK45)
    dOut{i}(1,:) = t;
    dOut{i}(2,:) = y(:,1);
    dOut{i}(3,:) = y(:,2);
end
% Plot 20 sample rays from the calculated rays.
for k = 0:19
    k = floor(k*numRays/20+1);
    plot(dOut{k}(2,:),dOut{k}(3,:))
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
    yfs(j) = dOut{j}(3,end);
end
% Show density of particles on source and image plane
ksdensity(y0s,'Bandwidth', 0.01);
ksdensity(yfs,'Bandwidth', 0.01);
legend({'Source Density', 'Image Density'},'Location','northwest')
% Define particle eqns of motion
function dydt = trace(~,y, res,m,q)
qp = [y(1);y(2)]; % Current position vector
[Ex, Ey] = evaluateGradient(res,qp); % Interpolate electric field at qp
dydt = [y(3); y(4); q/m*Ex; q/m*Ey]; % Linear system of eqns for Lorentz force
end
% Create event which stops ODE solver
function [value,isterminal,direction] = myEventsFcn(~,y,xmax)
value = 1;
if y(1)>(0.95*xmax) % check if particle is nearing edge of geometry
    value = 0;
end
isterminal = 1;
direction = 0;
end