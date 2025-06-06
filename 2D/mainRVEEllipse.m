% =========================================================================
% A 2D plane strain RVE with elliptical inclusion.
% The results are compared with those using COMSOL
%
% Author: Wang Zhetong
% =========================================================================
%% clear
clear;  close all;  clc;    
disp(datetime('now'))
%% addpath
addpath('./cfiles_New','./IGA_mesh','./Matrices_assemble','./Picture_plot','./solve','.')
%% Define Material Parameters.
E=[1280,12800];          % Instantaneous elastic modulus (matrix, inclusion)
nu=1/3;                  % Instantaneous Poisson's ratio
tau=10;                  % Relaxation time
mu=E/(2*(1+nu));         % Shear modulus for each material
K=E/(3*(1-2*nu));        % Volumetric modulus for each material
visG=0.25;               % Normalized shear modulus
visK=0;                  % Normalized volumetric modulus
%% Define Geometry Parameters.
fprintf(['----------------------------------\n' 'Generate NURBS\n'])

d=2;  % Problem dimension

% Geometry setup: RVE with elliptical inclusion
% Rectangle: [x0, y0, width, height]; Ellipse: [x0, y0, axis_x, axis_y]
infoGeo={[0.5,0.5,0.5,0.5], [0.5,0.5,0.3,0.2]};
geomData=generateGeometryRVE({'rec','ell'}, infoGeo);

% NURBS parameters
numPlot=1;
refinement=cell(1,2);
ref=[9,6];  % Refinement levels for outer and inner boundary
for i=1:length(infoGeo)
    refinement{i}=ref(i)*ones(geomData(i).numSeg,1);
end
numBoundary=length(infoGeo);         % Number of boundaries
p=2*ones(numBoundary,1);             % NURBS order (quadratic)
isClosed=ones(numBoundary,1);        % Closed curves
nurbsStr=getStructure2D(geomData,p,E,nu,refinement,isClosed,numPlot);

% Generate inner points
% numInner = [numX, numY, numR, numTheta, scale, tolerance]
numInner=[10,10,3,4*4,0.2,0.1];
innerPts=getInnerPointsRVE("RVE",numInner,infoGeo);
plot(innerPts(:,1),innerPts(:,2),'.k','MarkerSize',20,'DisplayName','Inner Points')

% Collect applied points for RIM (collocation + inner)
appliedPts=innerPts;
for i=1:numBoundary
    [appliedPts]=[nurbsStr(numBoundary-i+1).collocPts;appliedPts];
end
%% Assemble elastic system matrices.
fprintf(['----------------------------------\n' 'Assemble elastic matrix for RVE\n']);

numGauR=20;                   % Number of Gauss points for singular integral
numGauS=20;                   % Number of Gauss points for regular integral

% IGA transform matrix (displacement/traction)
[tranU,tranT]=transformMatrix2D(nurbsStr);

% Assemble boundary integral matrices
tElaBou=tic;
[HElaBou,GElaBou,HElaBouOri,GElaBouOri]=elasticMatrixBoundaryRVE2D(nurbsStr,numGauR,numGauS,appliedPts);
toc(tElaBou)

% Assemble domain integral matrices
tElaIn=tic;
[HElaIn,GElaIn]=elasticMatrixInsideRVE2D(nurbsStr,numGauR,innerPts);
toc(tElaIn)
%% Solve elastic problem using Periodic Boundary Condition.
fprintf(['----------------------------------\n' 'Solve elastic problem\n']);

% Define macro strain
epsilon0=1e-3;
PBC=[0,1;1,0]*epsilon0;                % ε₁₂=ε₂₁=1e-3

% Apply periodic boundary condition
[preU,preT,unU,unT,vecT,vecU,U32,U14,T32,T14]=setPBC(PBC,nurbsStr);

% Assemble full system with transformation and periodicity
[H,G,vecP]=AssembleMatrix(HElaBou,GElaBou,tranU,tranT,U32,U14,T32,T14,PBC);
HOri=HElaBouOri/tranU;
GOri=GElaBouOri/tranT;

% Solve boundary system
QOut=zeros(size(vecU,1),1)+vecP;
[Disp,Trac]=solveBoundaryCondition(H,G,QOut,preU,preT,unU,unT,vecT,vecU);

% Recover physical displacement and traction
[Disp,Trac]=recoverUT(Disp,Trac,U32,U14,T32,T14,HOri,GOri,QOut,vecP,PBC,nurbsStr);
%% Deformed picture.
fprintf(['----------------------------------\n' 'Deformed picture\n']);

factor=100;                        % Scaling factor for visualization
numPlot=numPlot+1;
plotDeformed(nurbsStr,tranU\Disp,factor)
%% Assemble viscoelastic system matrices.
fprintf(['----------------------------------\n' 'Assemble viscoelastic matrix\n']);

dA=2;                           % Support radius

% RIM transform matrix
tRim=tic;
rimPhi=rimMatrixPhi(appliedPts,appliedPts,dA,d);
toc(tRim)
rimInvPhi=inv(rimPhi);
rimInvPhi=rimInvPhi(:,1:end-d*(1+d));

% Stiffness matrix for viscoelastic problem
tVis=tic;
[HVis1,HVis2]=viscoelasticMatrix2D(nurbsStr,numGauR,numGauS,appliedPts,dA);
toc(tVis)
%% Solve viscoelastic problem.
% In the paper "IGABEM for the Homogenization of Linear Viscoelastic Composites",
% displacement on boundary is solved first.

fprintf(['----------------------------------\n' 'Solve viscoelastic process\n']);

deltaT=0.1;                     % Time step size
numT=50/deltaT;                 % Number of steps
eTaoI=exp(-deltaT./tau);        % e^(-Δt/τi)
t=0:deltaT:numT*deltaT;         % Time vector

epsilon0=1e-3;                  % Macroscopic strain magnitude
avgStrain={[1 0;0 0];[0 0;0 1];[0 1/2;1/2 0]};  % ε11, ε22, ε12
avgStress=zeros(length(t),3,3);                % Stress storage

for k=1:3
    PBC=avgStrain{k}*epsilon0;
    [preU,preT,unU,unT,vecT,vecU,U32,U14,T32,T14]=setPBC(PBC,nurbsStr);

    % Assemble boundary system with periodicity and transformation
    [H,G,vecP]=AssembleMatrix(HElaBou,GElaBou,tranU,tranT,U32,U14,T32,T14,PBC);
    HOri=HElaBouOri/tranU;
    GOri=GElaBouOri/tranT;
    HIn=HElaIn/tranU;
    GIn=GElaIn/tranT;

    % Matrix decomposition
    HG=(HVis1-1/3*HVis2)*rimInvPhi;
    HK=HVis2*rimInvPhi;
    QGi=zeros(d*size(appliedPts,1),length(visG));
    QKi=zeros(d*size(appliedPts,1),length(visK));
    [U1,U2]=deal(zeros(d*size(appliedPts,1),1));
    QOut=zeros(size(vecU,1),1)+vecP;
    QIn=zeros(size(innerPts,1)*d,1);

    % Time integration
    for i=1:length(t)
        [Disp,Trac]=solveBoundaryCondition(H,G,QOut,preU,preT,unU,unT,vecT,vecU);
        [Disp,Trac]=recoverUT(Disp,Trac,U32,U14,T32,T14,HOri,GOri,QOut,vecP,PBC,nurbsStr);

        % Displacement inside domain
        DispIn=-HIn*Disp+GIn*Trac+QIn;
        U2=U1;
        U1=[Disp;DispIn];

        [QGi,QKi,Q]=QViscoelastic(QGi,QKi,U1,U2,HG,HK,deltaT,eTaoI,i,mu(1),K(1),tau,visG,visK);
        QOut=Q(1:size(vecU,1))+vecP;
        QIn=Q(size(vecU,1)+1:end);

        % Homogenized stress
        avgStress(i,:,k)=getAverageStress(nurbsStr,tranT\Trac);
    end
end
%% Get the exponential fitting in Abaqus.
n=1;
s1111=[avgStress(2:n:end, 1, 1)/avgStress(1, 1, 1) t(2:n:end)'];
s2222=[avgStress(2:n:end, 2, 2)/avgStress(1, 2, 2) t(2:n:end)'];
s1212=[avgStress(2:n:end, 3, 3)/avgStress(1, 3, 3) t(2:n:end)'];
%% Plot Average Stress in viscoelastic progress.
fprintf(['----------------------------------\n' 'Plot viscoelastic process\n']);

filename1={'./txtData/CASE1_E11.txt';'./txtData/CASE1_E22.txt';'./txtData/CASE1_E12.txt'};

load color_RGBK.mat
lineStys={'-.','--',':'};
shapeStys={'s','o','^'};

numPlot=numPlot+1;
figure(numPlot);

idxFEM=[2;3;4];         % Column indices in FEM file
idxBEM=[1;2;3];         % Stress component index
n=2;                    % FEM sampling interval

for i=1:3
    % Plot IGABEM results
    plot(t',avgStress(:,idxBEM(i),i)/epsilon0,'Color',color_RGBK{i},'LineWidth',1.5,...
        'LineStyle',lineStys{i})
    hold on

    % Plot FEM results
    sFEM=readmatrix(filename1{i});
    plot(sFEM(1:n:end-1,1),sFEM(1:n:end-1,idxFEM(i))/epsilon0,...
        'MarkerEdgeColor',color_RGBK{i},'Marker',shapeStys{i},'LineStyle','none',...
        'MarkerSize',5,'LineWidth',1)
    hold on
end

% Plot info
lgd=legend('$${\rm{\bar{C}_{1111}\,(IGABEM)}}$$','$${\rm{\bar{C}_{1111}\,(COMSOL)}}$$',...
    '$${\rm{\bar{C}_{2222}\,(IGABEM)}}$$','$${\rm{\bar{C}_{2222}\,(COMSOL)}}$$',...
    '$${\rm{\bar{C}_{1212}\,(IGABEM)}}$$','$${\rm{\bar{C}_{1212}\,(COMSOL)}}$$');
xtitle='$${\mathrm{Time(days)}}$$';
ytitle='$${\mathrm{\bar{C}_{ijkl}(MPa)}}$$';
setPlot(xtitle,ytitle,lgd)
%% Rmpath
rmpath('./cfiles_New','./IGA_mesh','./Matrices_assemble','./Picture_plot','./solve','.')




