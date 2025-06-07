% =========================================================================
% RVE with a spherical inclusion.
% The results are compared with those using COMSOL and those using 
% asymptotic homogenization method (AHM).
% 
% AHM: On the effective behavior of viscoelastic composites in three 
% dimensions
% 
% Author: [Wang Zhetong]
% =========================================================================
%% clear
clear;close all;clc;
disp(datetime('now'))
%% addpath
addpath('./cfiles_New','./IGA_mesh','./Matrices_assemble','./Picture_plot','./solve','.')
%% Define Material Parameters.
E=[4.082e3;168.4e3];       % Instantaneous elastic modulus (MPa)
% the actual nu for inclusion is 0.443
nu=0.311*ones(length(E),1);   % Instantaneous Poisson's ratio
tau=159.81/4.082;             % relaxation times (min)
visG=1;                       % normalized shear modulus
visK=1;                       % normalized volume modulus
Vs=[10 20 30];                % Volume fractions
PBs={[0 0 1;0 0 0;1 0 0]/2,[1 0 0;0 0 0;0 0 0]}; % Average strains
nameE={'e13','e11'};         % Names for different load cases
%% 6 cases
for idxV=1:3
    %% Define Geometry Parameters.
    fprintf(['----------------------------------\n' 'Generate NURBS\n']);

    dim=3;                        % Dimension
    % A cube matrix with a spherical inclusion.
    matL=[0.5,0.5,0.5];           % Matrix size
    matO=[0.5,0.5,0.5];           % Matrix center
    V=Vs(idxV)/100;               % Inclusion volume fraction
    incO=[0.5,0.5,0.5];           % Inclusion center
    incR=(3*V/(4*pi))^(1/3);      % Inclusion radius
    infoGeo={{'cub',[matL,matO]};{'sph6',[incR,incO]}};

    % NURBS parameters
    numPlot=1;
    refinement=[5,5,5;2*ones(size(incO,1),3)];
    % Number of boundaries
    numBoundary=size(infoGeo,1);
    % Discontinuous element parameter.
    lamda=[0.97;0.965*ones(size(incO,1),1)];
    nurbsStr=getStructure3D(infoGeo,E,nu,refinement,lamda,numPlot);
    [nurbsStr,collocNew]=generateGlobalCollocPoints(nurbsStr,numPlot);

    % Inner points
    innerPoints=[];

    % Applied points for RIM
    appliedPts=innerPoints;
    for i=1:numBoundary
        [appliedPts]=[collocNew(numBoundary-i+1).collocPts;appliedPts];
    end
    %% Assemble elastic system matrices.
    fprintf(['----------------------------------\n' 'Assemble elastic matrix\n']);

    % IGA Transform Part
    [tranU,tranT]=transformMatrix3D(nurbsStr,collocNew);

    numGauR=12; % for singular integral
    numGauS=12; % for regular integral

    % Elastic Part
    % Boundary
    tElaBou=tic;
    [HElaBou,GElaBou,HOrig,GOrig]=elasticMatrixBoundaryRVE3D(nurbsStr,...
        collocNew,numGauR,numGauS,tranU,appliedPts);
    toc(tElaBou)

    % Inside
    tElaIn=tic;
    [HElaIn,GElaIn]=elasticMatrixInsideRVE3D(nurbsStr,numGauR,innerPoints,collocNew);
    toc(tElaIn)

    for idxBC=1:2
        %% Solve viscoelastic problem.
        fprintf(['----------------------------------\n' 'Solve viscoelastic process\n']);

        % Time parameters
        dt=0.5;
        nt=150/dt;
        eTaoI=exp(-dt./tau); % e^(-Δt/τi)
        t=0:dt:nt*dt;

        % Boundary condition.
        PB=PBs{idxBC};
        epsilon0=1e-3;
        epsilonPB=epsilon0*PB;

        % Set periodic boundary conditions
        [preU,preT,unU,unT,vecT,vecU,face325U,face146U,face325T,face146T]=...
            setPeriodicBoundary3D(nurbsStr,collocNew,epsilonPB);

        % Assemble system matrices
        [H,G,HO,GO,HIn,GIn]=assembleSystemMatrices(HElaBou,GElaBou,...
            HOrig,GOrig,HElaIn,GElaIn,tranU,tranT);

        % Calculate the vector for the periodic boundary
        [H,G,vecP]=rebulidPeriodicMatrix(H,G,face325U,face146U,face325T, ...
            face146T,epsilonPB);

        % Viscoelastic matrix
        % According to Ref 'A boundary element methodology for viscoelastic
        % analysis: Part I with cells.', viscoelastic matrix can be
        % descirbed by elastic one for the special viscoelastic
        % parameter.(nu keeps constant)
        HVis=HElaBou/tranU;
        Qi=zeros(dim*size(appliedPts,1),length(visG));
        [U1,U2]=deal(zeros(dim*size(appliedPts,1),1));
        Q=zeros(size(collocNew,1),1)+vecP;
        averStress=zeros(length(t),6);

        % Viscoelastic progress
        for i=1:length(t)
            % Solve boundary condition
            [Disp,Trac]=solveBoundaryCondition(H,G,Q,preU,preT,unU,unT,vecT,vecU);
            
            % Apply periodic boundary condition to results
            [Disp,Trac]=computePeriodicValues(Disp,Trac,face325U,face146U,face325T,...
                face146T,epsilonPB,nurbsStr,collocNew,HO,GO,Q,vecP);
            
            % The displacement on the colloction points.
            U2=U1;
            U1=Disp;
            Q=zeros(size(collocNew,1),1)+vecP;
            Qi=eTaoI'.*Qi+(visG.*(1-(tau/dt-1).*(1-eTaoI)))'.*(HVis*U1)...
                -(i>1)*(visG.*(1-tau/dt.*(1-eTaoI)))'.*(HVis*U2);
            Q=sum(Qi,2)+Q;

            % Compute average stress.
            averStress(i,:)=getAverageStress3D(nurbsStr,tranT\Trac);
        end

        filename=['./matData/V',num2str(Vs(idxV)),'_',nameE{idxBC},'.mat'];
        save(filename,'averStress','t')
    end
end
%% Plot
load color_RGBK.mat
load data_from_origin.mat

Vs=[30 20 10];
nameE={'e13','e11','e11'};
idxBEM=[6 1 2];
idxFEM=[7 2 3];
n=1;
epsilon0=1e-3;
ytitles={'$${\mathrm{\bar{C}_{1313}(GPa)}}$$','$${\mathrm{\bar{C}_{1111}(GPa)}}$$',...
    '$${\mathrm{\bar{C}_{1122}(GPa)}}$$'};
ylims=[3 12 4];
lineStys={':','--','-.'};
shapeStys={'s','o','^'};

for idx1=1:3
    figure(idx1+1)
    for idx2=1:3
        filename=['./matData/V',num2str(Vs(idx2)),'_',nameE{idx1},'.mat'];
        load(filename,'averStress','t')
        plot(t,averStress(:,idxBEM(idx1)),'Color',color_RGBK{idx2},...
            'LineWidth',2.5,'LineStyle','-');
        hold on
        filename1=['./txtData/',nameE{idx1},'_v',num2str(Vs(idx2)),'.txt'];
        sFEM=readmatrix(filename1);
        plot(sFEM(1:n:end,1),sFEM(1:n:end,idxFEM(idx1))/epsilon0*1e-3,...
            'Color',color_RGBK{idx2},'LineWidth',2.5,'LineStyle',lineStys{idx2})
        hold on
        plot(AHM{idx1}(:,idx2*3-2),AHM{idx1}(:,idx2*3-1),...
            'MarkerEdgeColor',color_RGBK{4},'Marker',shapeStys{idx2},...
            'LineStyle','none','MarkerSize',11,'LineWidth',1.5)
        hold on
    end
    xlim([0 150])
    ylim([0 ylims(idx1)])
    xtitle='$${\mathrm{Time(min)}}$$';
    ytitle=ytitles{idx1};
    lgd=legend('$${\rm{{V}_{f} = 30\% \,(IGABEM)}}$$',...
        '$${\rm{{V}_{f} = 30\% \,(COMSOL)}}$$',...
        '$${\rm{{V}_{f} = 30\% \,(AHM)}}$$',...
        '$${\rm{{V}_{f} = 20\% \,(IGABEM)}}$$',...
        '$${\rm{{V}_{f} = 20\% \,(COMSOL)}}$$',...
        '$${\rm{{V}_{f} = 20\% \,(AHM)}}$$',...
        '$${\rm{{V}_{f} = 10\% \,(IGABEM)}}$$',...
        '$${\rm{{V}_{f} = 10\% \,(COMSOL)}}$$',...
        '$${\rm{{V}_{f} = 10\% \,(AHM)}}$$');
    setPlot(xtitle,ytitle,lgd)
end

%% Rmpath.
rmpath('./cfiles_New','./IGA_mesh','./Matrices_assemble','./Picture_plot','./solve','.')








