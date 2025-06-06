% =========================================================================
% A 2D plane strain cylinder under pressure. The results are compared with
% those in "RI-IGABEM for 2D Viscoelastic Problems and Its Application to
% Solid Propellant Grains"
%
% Author: Wang Zhetong
% =========================================================================
%% clear
clear;  close all;  clc;
disp(datetime('now'))
%% addpath
addpath('./cfiles_New','./IGA_mesh','./Matrices_assemble','./Picture_plot','./solve','.')
%% Define Material Parameters.
mu=1e4;                % shear modulus
K=2e4;                 % volumetric modulus
E=9*K*mu/(3*K+mu);     % Instantaneous elastic modulus (MPa)
nu=(3*K-2*mu)/(2*(3*K+mu)); % Instantaneous Poisson's ratio
tau=1/2.3979;          % relaxation time
visG=0.99;             % normalized shear modulus
visK=0;                % normalized bulk modulus
%% Define Geometry Parameters.
fprintf(['----------------------------------\n' 'Generate NURBS\n']);

d=2;                          % problem dimension

% Define geometry: thick-walled cylinder (center (0,0), inner R=5, outer R=10)
infoGeo = {[0,0,5,10]};
geomData = generateGeometry({'cyl'}, infoGeo);

% NURBS parameters
numPlot=1;
refinement={9*ones(geomData.numSeg,1)};      % uniform refinement per segment
numBoundary=length(infoGeo);                 % number of boundary curves
p=2*ones(numBoundary,1);                     % NURBS order
isClosed=ones(numBoundary,1);                % closed curve indicator
nurbsStr=getStructure2D(geomData,p,E,nu,refinement,isClosed,numPlot);

% Different inner point configurations to test convergence
numIPs = {[3 3], [5 5], [8 8], [10 20]};
lineStys = {'-.', '--', '-', ':'};
load color_RGBK.mat

for idx=1:4
    % Generate inner points
    numIP = numIPs{idx};
    innerPts = getInnerPoints("cyl", numIP, infoGeo{1});
    
    if idx==1
        figure(1)
        plot(innerPts(:,1),innerPts(:,2),'.k','MarkerSize',20,'DisplayName','Inner Points')
    end

    % Generate applied points for RIM (boundary collocation + inner points)
    appliedPts=innerPts;
    for i=1:numBoundary
        [appliedPts]=[nurbsStr(numBoundary-i+1).collocPts;appliedPts];
    end
    %% Assemble elastic system matrices.
    fprintf(['----------------------------------\n' 'Assemble elastic matrix\n']);

    numGauR=12;      % Gauss points for singular integral
    numGauS=12;      % Gauss points for regular integral

    [tranU,tranT]=transformMatrix2D(nurbsStr);  % local-global transform

    tElaBou=tic;
    [HElaBou,GElaBou]=elasticMatrixBoundary2D(nurbsStr,numGauR,numGauS,appliedPts);
    toc(tElaBou)

    tElaIn=tic;
    [HElaIn,GElaIn]=elasticMatrixInside2D(nurbsStr,numGauR,innerPts);
    toc(tElaIn)
    %% Assemble viscoelastic system matrices.
    fprintf(['----------------------------------\n' 'Assemble viscoelastic matrix\n']);

    dA=10;  % support radius for RIM

    tRim=tic;
    rimPhi=rimMatrixPhi(appliedPts,appliedPts,dA,d);  % radial basis matrix
    toc(tRim)

    rimInvPhi=inv(rimPhi);                   % inverse of support matrix
    rimInvPhi=rimInvPhi(:,1:end-d*(1+d));    % drop nullspace cols

    tVis=tic;
    [HVis1,HVis2]=viscoelasticMatrix2D(nurbsStr,numGauR,numGauS,appliedPts,dA); % stiffness matrices
    toc(tVis)
    %% Solve viscoelastic problem.
    fprintf(['----------------------------------\n' 'Solve viscoelastic process\n']);

    dt=0.03;
    tEnd=300;
    numT=tEnd/dt;
    eTaoI=exp(-dt./tau);     % e^(-Δt/τi)
    t=0:dt:numT*dt;
    dispTime=zeros(length(t),(size(innerPts,1)+nurbsStr.numCollocU)*d);
    tracTime=zeros(length(t),nurbsStr.numCollocT*d);

    % Boundary conditions
    BC = ["dy","f","dx","p"];
    pressure = 10;
    dirichlet = 0;
    [preU,preT,unU,unT,vecT,vecU]=setBC(BC,nurbsStr,pressure,dirichlet);

    % Final system matrices
    HVis=(2*mu*(HVis1-1/3*HVis2)+K*HVis2)*rimInvPhi;
    GVis=[GElaBou;GElaIn]/tranT;
    [unU]=[unU;(1:size(innerPts,1)*d)'+nurbsStr.numCollocU*d];

    HG=(HVis1-1/3*HVis2)*rimInvPhi;
    HK=HVis2*rimInvPhi;
    QGi=zeros(d*size(appliedPts,1),length(visG));
    QKi=zeros(d*size(appliedPts,1),length(visK));
    [Q,U1,U2]=deal(zeros(d*size(appliedPts,1),1));

    % Time-stepping loop
    for i=1:length(t)
        [Disp,Trac]=solveBoundaryCondition(HVis,GVis,Q,preU,preT,unU,unT,vecT,vecU);

        U2=U1;
        U1=Disp;
        [QGi,QKi,Q]=QViscoelastic(QGi,QKi,U1,U2,HG,HK,dt,eTaoI,i,mu,K,tau,visG,visK);
        dispTime(i,:)=Disp';
        tracTime(i,:)=Trac';
    end

    %% Plot displacement at specific inner point
    fprintf(['----------------------------------\n' 'Plot viscoelastic process\n']);

    if idx==1
        numPlot=numPlot+1;
    end
    figure(numPlot);
    idxPt=6;  % choose middle inner point
    plot(t,dispTime(:,idxPt*2-1),'Color',color_RGBK{idx},'LineWidth',1.5 ...
        ,'LineStyle',lineStys{idx});
    hold on
end
%% Plot analytical solution
tA=0:10:tEnd;
GInf=mu*(1-visG);             % long-term shear modulus
G1=mu*visG;                   % transient part of shear modulus

% Analytical solution at point (7.5, 0)
dispA = CylinderAnalysis([7.5,0], infoGeo{1}(3), infoGeo{1}(4), -pressure, GInf, G1, nu, tA, tau);

% Plot analytical result
plot(tA,dispA(:,1),'ko','MarkerSize',5,'LineWidth',1)
hold on

% Add legend and axis labels
lgd= legend('9 nodes','25 nodes','64 nodes','200 nodes', 'Analytical Results');
xtitle='Time (s)';   ytitle='Displacement (mm)';
setPlot(xtitle,ytitle,lgd)
%% Reference
% Plot reference data from previous validated results
load ref_in_pts.mat

numPlot=numPlot+1;
figure(numPlot);
for idx=1:4
    plot(ref_data{idx}(:,1),ref_data{idx}(:,2),'Color',color_RGBK{idx}, ...
        'LineWidth',1.5,'LineStyle',lineStys{idx})
    hold on
end

% Plot analytical result again
plot(tA,dispA(:,1),'ko','MarkerSize',5,'LineWidth',1)
hold on

% Add legend and axis labels
lgd= legend('9 nodes (ref)','25 nodes (ref)','64 nodes (ref)', ...
    '200 nodes (ref)', 'Analytical Results');
xtitle='Time (s)';   ytitle='Displacement (mm)';
setPlot(xtitle,ytitle,lgd)
%% Rmpath
% Clean up environment
rmpath('./cfiles_New','./IGA_mesh','./Matrices_assemble','./Picture_plot','./solve','.')


