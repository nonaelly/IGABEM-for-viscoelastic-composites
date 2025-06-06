% =========================================================================
% A 2D plane strain plate under tension. The results are compared with
% those in "RI-IGABEM for 2D Viscoelastic Problems and Its Application to
% Solid Propellant Grains"
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
mu=1e4;                % shear modulus
K=2e4;                 % volumetric modulus
E=9*K*mu/(3*K+mu);     % Instantaneous elastic modulus (MPa)
nu=(3*K-2*mu)/(2*(3*K+mu)); % Instantaneous Poisson's ratio
tau=1/2.3979;          % relaxation times
visG=0.99;             % normalized shear modulus
visK=0;                % normalized volume modulus
%% Define Geometry Parameters.
fprintf(['----------------------------------\n' 'Generate NURBS\n']);

d=2;                          % Dimension

% Plate (0,0)-(10,5)
infoGeo={[0,0,10,5]};
geomData=generateGeometry({'rec'}, infoGeo);

% NURBS parameters
numPlot=1;
refinement={3*ones(geomData.numSeg,1)};
numBoundary=length(infoGeo);           % Number of boundaries
p=2*ones(numBoundary,1);               % NURBS order
isClosed=ones(numBoundary,1);          % Indicates whether the curve is closed or open
nurbsStr=getStructure2D(geomData,p,E,nu,refinement,isClosed,numPlot);

% different time step size
dts=[0.3,0.8,1,3];
lineStys = {'-.', '--', '-', ':'};
load color_RGBK.mat
for idx=1:4
    % Inner points
    innerPts=[6,2.5];
    if idx==1
        figure(1)
        plot(innerPts(:,1),innerPts(:,2),'.k','MarkerSize',20, ...
            'DisplayName','Inner Points')
    end

    % Applied points for RIM
    appliedPts=innerPts;
    for i=1:numBoundary
        [appliedPts]=[nurbsStr(numBoundary-i+1).collocPts;appliedPts];
    end
    %% Assemble elastic system matrices.
    fprintf(['----------------------------------\n' 'Assemble elastic matrix\n']);

    numGauR=12;                  % for singular integral
    numGauS=12;                  % for regular integral

    % IGA Transform Part
    [tranU,tranT]=transformMatrix2D(nurbsStr);

    % Elastic Part
    % Boundary
    tElaBou=tic;
    [HElaBou,GElaBou]=elasticMatrixBoundary2D(nurbsStr,numGauR,numGauS,appliedPts);
    toc(tElaBou)

    tElaIn=tic;
    % Inside
    [HElaIn,GElaIn]=elasticMatrixInside2D(nurbsStr,numGauR,innerPts);
    toc(tElaIn)
    %% Assemble viscoelastic system matrices.
    fprintf(['----------------------------------\n' 'Assemble viscoelastic matrix\n']);

    dA=10;                       % Support radius

    % Transform matrix in radial integral method.
    tRim=tic;
    rimPhi=rimMatrixPhi(appliedPts,appliedPts,dA,d);
    toc(tRim)
    rimInvPhi=inv(rimPhi);
    rimInvPhi=rimInvPhi(:,1:end-d*(1+d));

    % Stiffness matrix for viscoelastic problem.
    tVis=tic;
    [HVis1,HVis2]=viscoelasticMatrix2D(nurbsStr,numGauR,numGauS,appliedPts,dA);
    toc(tVis)
    %% Solve viscoelastic problem.
    fprintf(['----------------------------------\n' 'Solve viscoelastic process\n']);

    dt=dts(idx);
    tEnd=300;
    numT=tEnd/dt;
    eTaoI=exp(-dt./tau);     % e^(-Δt/τi)
    t=0:dt:numT*dt;
    dispTime=zeros(length(t),(size(innerPts,1)+nurbsStr.numCollocU)*d);
    tracTime=zeros(length(t),nurbsStr.numCollocT*d);

    % Boundary condition
    BC=["dy","p","f","dx"];
    pressure=-100;
    dirichlet=0;
    [preU,preT,unU,unT,vecT,vecU]=setBC(BC,nurbsStr,pressure,dirichlet);

    % Assembly H and G
    HVis=(2*mu*(HVis1-1/3*HVis2)+K*HVis2)*rimInvPhi;
    GVis=[GElaBou;GElaIn]/tranT;

    % Add displacement unknown (inner points)
    [unU]=[unU;(1:size(innerPts,1)*d)'+nurbsStr.numCollocU*d];

    % Matrix for viscoelastic problem
    HG=(HVis1-1/3*HVis2)*rimInvPhi;
    HK=HVis2*rimInvPhi;
    QGi=zeros(d*size(appliedPts,1),length(visG));
    QKi=zeros(d*size(appliedPts,1),length(visK));
    [Q,U1,U2]=deal(zeros(d*size(appliedPts,1),1));

    % Viscoelastic time loop
    for i=1:length(t)
        [Disp,Trac]=solveBoundaryCondition(HVis,GVis,Q,preU,preT,unU,unT,vecT,vecU);

        U2=U1;
        U1=Disp;
        [QGi,QKi,Q]=QViscoelastic(QGi,QKi,U1,U2,HG,HK,dt,eTaoI,i,mu,K,tau,visG,visK);
        dispTime(i,:)=Disp';
        tracTime(i,:)=Trac';
    end
    %% Plot viscoelastic progress
    fprintf(['----------------------------------\n' 'Plot viscoelastic process\n']);

    if idx==1
        numPlot=numPlot+1;
    end
    figure(numPlot);
    idxPt=size(appliedPts,1);
    plot(t,dispTime(:,idxPt*2-1),'Color',color_RGBK{idx},'LineWidth',1.5, ...
        'LineStyle',lineStys{idx});
    hold on
end
%% Plot analytical solution
tA=0:10:tEnd;
GInf=mu*(1-visG);
G1=mu*visG;
dispA=RectangleAnalysis([6,2.5],infoGeo{1}(3),infoGeo{1}(4),pressure,GInf,G1,nu,tA,tau);

plot(tA,dispA(:,1),'ko','MarkerSize',5,'LineWidth',1)
hold on

% Plot info
lgd=legend('$${\rm{\Delta t=0.3}}$$','$${\rm{\Delta t=0.8}}$$', ...
    '$${\rm{\Delta t=1.0}}$$','$${\rm{\Delta t=3.0}}$$','Analytical Results');
xtitle='Time (s)';   ytitle='Displacement (mm)';
setPlot(xtitle,ytitle,lgd)

%% Reference
load ref_time_step.mat

numPlot=numPlot+1;
figure(numPlot);
for idx=1:4
    plot(ref_data{idx}(:,1),ref_data{idx}(:,2),'Color',color_RGBK{idx}, ...
        'LineWidth',1.5,'LineStyle',lineStys{idx})
    hold on
end
plot(tA,dispA(:,1),'ko','MarkerSize',5,'LineWidth',1)
hold on

% Plot info
lgd=legend('$${\rm{\Delta t=0.3 (ref)}}$$','$${\rm{\Delta t=0.8 (ref)}}$$', ...
    '$${\rm{\Delta t=1.0 (ref)}}$$','$${\rm{\Delta t=3.0 (ref)}}$$','Analytical Results');
xtitle='Time (s)';   ytitle='Displacement (mm)';
setPlot(xtitle,ytitle,lgd)

%% Rmpath
rmpath('./cfiles_New','./IGA_mesh','./Matrices_assemble','./Picture_plot','./solve','.')