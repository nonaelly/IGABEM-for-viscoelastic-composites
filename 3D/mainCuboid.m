% =========================================================================
% A 3D cuboid under tension. 
% The results are compared with analytical solutions.
%
% Author: Wang Zhetong
% Date: 2025-06-06
% =========================================================================
%% clear
clear;  close all;  clc;    
disp(datetime('now'))
%% addpath
addpath('./cfiles_New','./IGA_mesh','./Matrices_assemble','./Picture_plot','./solve','.')
%% Define Material Parameters.
E = 1280;            % Instantaneous elastic modulus (MPa)
nu = 0.4;            % Instantaneous Poisson's ratio
tau = 2.5;           % relaxation times
mu = E/(2*(1+nu));   % shear modulus
K = E/(3*(1-2*nu));  % volumetric modulus
visG = 0.25;         % nomalized shear modulus
visK = 0;            % nomalized volume modulus
%% Define Geometry Parameters.
fprintf(['----------------------------------\n' 'Generate NURBS\n']);

d = 3; % Dimension
% Cuboid for example.
% Center points (0,0,0); length (a,b,c) .
infoGeo = {[15,20,25,15,20,25]};
geomData = {{'cub', infoGeo{1}}};

% NURBS parameters
numPlot = 1;
refinement = [2, 3, 4];
numBoundary = length(infoGeo);  % Number of boundaries
lamda = 1;
nurbsStr = getStructure3D(geomData, E, nu, refinement, lamda, numPlot);
[nurbsStr,collocNew] = generateGlobalCollocPoints(nurbsStr, numPlot);

% Inner points
numIP = [1, 1, 1];
innerPts = getInnerPoints3D(geomData{1}, numIP);
plot3(innerPts(:, 1), innerPts(:, 2), innerPts(:, 3), '.k', ...
    'MarkerSize', 20, 'DisplayName','Inner Points')

% Applied points for RIM
appliedPts = innerPts;
for i = 1:numBoundary
    [appliedPts] = [collocNew(numBoundary - i + 1).collocPts; appliedPts];
end
%% Assemble elastic system matrices.
fprintf(['----------------------------------\n' 'Assemble elastic matrix\n']);

numGauR = 12; % for singular integral
numGauS = 12; % for regular integral

% IGA Transform Part
[tranU,tranT] = transformMatrix3D(nurbsStr, collocNew);

% Elastic Part
% Boundary
tElaBou = tic;
[HElaBou, GElaBou] = elasticMatrixBoundary3D(nurbsStr, collocNew, numGauR, numGauS, tranU, appliedPts);
toc(tElaBou)

tElaIn = tic;
% Inside
[HElaIn, GElaIn] = elasticMatrixInside3D(nurbsStr, numGauR, innerPts, collocNew);
toc(tElaIn)
%% Solve elastic problem.
fprintf(['----------------------------------\n' 'Solve elastic problem\n']);

% Boundary condition.
BC = ["dy","px","f","dx","f","dz"];
pressure = 1;   
dirichlet = 1e-3;
[preU, preT, unU, unT, vecT, vecU] = setBC3D(BC,nurbsStr,pressure,dirichlet,collocNew);

% Assembly H and G; Transform the displacement and traction.
H = HElaBou/tranU; 
G = GElaBou/tranT; 
Q = zeros(size(vecU,1),1);
[Disp, Trac] = solveBoundaryCondition(H, G, Q, preU, preT, unU, unT, vecT, vecU);
%% Deformed picture.
fprintf(['----------------------------------\n' 'Deformed picture\n']);

factor = 150;
numPlot = numPlot + 1;

% Extend u and t to all control points.
[ctrlDisp, ~] = expandToControlPoints(nurbsStr, collocNew, tranU\Disp, tranT\Trac);

plotDeformed3D(nurbsStr, ctrlDisp, factor, numPlot)
%% Stress contour.
% Extend u and t to all control points.
[ctrlDisp, ctrlTrac] = expandToControlPoints(nurbsStr, collocNew, tranU\Disp, tranT\Trac);

% Compute stress on the boundary.
[stress, strain]= getStressBoundary(nurbsStr,ctrlDisp,ctrlTrac);

% Plot stress contour.
numPlot = numPlot + 1;
plotStress3D(nurbsStr, [ctrlDisp, ctrlDisp], 'xx', numPlot, 1:6)
numPlot = numPlot + 1;
plotStress3D(nurbsStr, [ctrlDisp, ctrlDisp], 'yy', numPlot, 1:6)
numPlot = numPlot + 1;
plotStress3D(nurbsStr, [ctrlDisp, ctrlDisp], 'zz', numPlot, 1:6)
%% Assemble viscoelastic system matrices.
fprintf(['----------------------------------\n' 'Assemble viscoelastic matrix\n']);

% Support radius.
dA = 100;

% Transform matrix in radial integral method.
tRim = tic;
[rimPhi] = rimMatrixPhi(appliedPts, appliedPts, dA, d);
toc(tRim)
rimInvPhi = inv(rimPhi);
rimInvPhi = rimInvPhi(:, 1 : end-d*(1+d));

% Stinffness matrix for viscoelastic problem.
tVis = tic;
[HVis1, HVis2] = viscoelasticMatrix3D(nurbsStr, numGauR, numGauS, appliedPts, dA, collocNew);
toc(tVis)
%% Solve viscoelastic problem.
% In paper 'A modified RI-IGABEM with only weakly singular integral for
% viscoelastic analysis', the displacement on boundary and inside are
% calculated together.

fprintf(['----------------------------------\n' 'Solve viscoelastic process\n']);

% Time parameters
deltaT = 0.01; 
numT = 30/deltaT;
eTaoI = exp(-deltaT./tau);	% e^(-Δt/τi)
t = 0 : deltaT : numT*deltaT;

[indU, ~, indNumU, indNumT, indNumB] = initializeIndices(nurbsStr, collocNew);
dispTime = zeros(length(t), size(appliedPts,1) * d);
tracTime = zeros(length(t), indNumT(end) * d);

% Boundary condition.
BC = ["dy","px","f","dx","f","dz"];
pressure = 1;  
dirichlet = 1e-3;
[preU, preT, unU, unT, vecT, vecU] = setBC3D(BC,nurbsStr,pressure,dirichlet,collocNew);

% Assembly H and G (different for elastic and viscoelastic problem).
HVis = ( 2* mu*(HVis1 - 1/3*HVis2) + K*HVis2 )*rimInvPhi;
GVis = [GElaBou; GElaIn]/tranT;

% Change the displacement unknown (inner points).
unU = [unU; (1 : size(innerPts,1)*d)' + indNumU(end) * d];

% Matrix for viscoelastic problem.
HG = (HVis1 - 1/3*HVis2)*rimInvPhi;
HK = HVis2*rimInvPhi;
QGi = zeros(d * size(appliedPts,1), length(visG));
QKi = zeros(d * size(appliedPts,1), length(visK));
[Q, U1, U2] = deal(zeros(d*size(appliedPts,1), 1));

% Viscoelastic progress
for i = 1 : length(t)
    [Disp, Trac] = solveBoundaryCondition(HVis, GVis, Q, preU, preT, unU, unT, vecT, vecU);
    
    % The displacement on the colloction points.
    U2 = U1;
    U1 = Disp;
	[QGi, QKi, Q] = QViscoelastic(QGi, QKi, U1, U2, HG, HK, deltaT, eTaoI, i, mu, K, tau, visG, visK);
    QOut = Q(1 : size(vecU,1));
    QIn = Q(size(vecU,1)+1 : end);

    dispTime(i, :) = Disp';
	tracTime(i, :) = Trac';
end
%% Plot viscoelastic progress.
fprintf(['----------------------------------\n' 'Plot viscoelastic process\n']);

% (a, 0) u1.
tA = 0 : 2 : 30;
GInf = mu*(1-visG);
G1 = mu*visG;
dispA = CuboidAnalysis3D([2*infoGeo{1}(4),0,0], 2*infoGeo{1}(4:6), pressure, GInf, G1, nu, tA, tau);

load color_RGBK.mat
numPlot = numPlot + 1;
figure(numPlot);
indPts = (refinement(1) + 2)*d - d + 1;
plot(t, dispTime(:,indPts),'Color',color_RGBK{1},'LineWidth',1.5,'LineStyle','-.');
hold on
plot(tA,dispA(:,1),'ko','MarkerSize',5,'LineWidth',1)

% photo imformation
lgd= legend('BEM','Analytical Solution');
xtitle = 'Time (s)';    
ytitle = 'Displacement (mm)';
setPlot(xtitle, ytitle, lgd)
%% Rmpath.
rmpath('./cfiles_New','./IGA_mesh','./Matrices_assemble','./Picture_plot','./solve','.')