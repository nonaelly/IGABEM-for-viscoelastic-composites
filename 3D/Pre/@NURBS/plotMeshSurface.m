function plotMeshSurface(obj,plotPro)

noElementXiSpan = plotPro.NO;
noElementEtaSpan = plotPro.NO;

resolution_xi = noElementXiSpan*obj.numberElementsU + 1;
resolution_eta = noElementEtaSpan*obj.numberElementsV + 1;

uKnotVec = unique(obj.uKnot);
xiVec = zeros(resolution_xi,1);
for ixi = 1:obj.numberElementsU
    temp_pts = linspace(uKnotVec(ixi),uKnotVec(ixi+1),noElementXiSpan+1);
    xiVec((ixi-1)*noElementXiSpan + 1:ixi*noElementXiSpan+1) = temp_pts;
end

vKnotVec = unique(obj.vKnot);
etaVec = zeros(resolution_eta,1);
for ixi = 1:obj.numberElementsV
    temp_pts = linspace(vKnotVec(ixi),vKnotVec(ixi+1),noElementEtaSpan+1);
    etaVec((ixi-1)*noElementEtaSpan + 1:ixi*noElementEtaSpan+1) = temp_pts;
end

% number of distinct knot values
noKnotsU = length(uKnotVec);
noKnotsV = length(vKnotVec);

x1  = zeros(noKnotsU,resolution_eta);
y1  = zeros(noKnotsU,resolution_eta);
z1  = zeros(noKnotsU,resolution_eta);

x2  = zeros(noKnotsV,resolution_xi);
y2  = zeros(noKnotsV,resolution_xi);
z2  = zeros(noKnotsV,resolution_xi);

% NURBS curves of knot lines corresponding to xi direction

for uk=1:noKnotsU
    xi = uKnotVec(uk);
    for i=1:resolution_eta
        eta = etaVec(i);
        PhyCoords = obj.GetPhyCoords(xi,eta);
        x1(uk,i) = PhyCoords(1);
        y1(uk,i) = PhyCoords(2);
        z1(uk,i) = PhyCoords(3);
    end
end

for vk=1:noKnotsV
    eta = vKnotVec(vk);
    for i=1:resolution_xi
        xi  = xiVec(i);
        PhyCoords = obj.GetPhyCoords(xi,eta);
        x2(vk,i) = PhyCoords(1);
        y2(vk,i) = PhyCoords(2);
        z2(vk,i) = PhyCoords(3);
    end
end

hold on
plot3(x1',y1',z1','Color',plotPro.Color,'LineStyle',plotPro.LineStyle,'LineWidth',plotPro.LineWidth);
plot3(x2',y2',z2','Color',plotPro.Color,'LineStyle',plotPro.LineStyle,'LineWidth',plotPro.LineWidth);
hold on

if strcmp(plotPro.s1,'PlotCPs')
    
    plot3(obj.controlPts(:,1),obj.controlPts(:,2),obj.controlPts(:,3),'bo',...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','g',...
        'MarkerSize',7,'LineWidth',2);
    
    for i=1:obj.noPtsV
        idx = (i-1)*obj.noPtsU+1:i*obj.noPtsU;
        plot3(obj.controlPts(idx,1),obj.controlPts(idx,2),obj.controlPts(idx,3),'b--','LineWidth',1);
    end
    
    for i=1:obj.noPtsU
        idx = i:obj.noPtsU:i+obj.noPtsU*(obj.noPtsV-1);
        plot3(obj.controlPts(idx,1),obj.controlPts(idx,2),obj.controlPts(idx,3),'b--','LineWidth',1);
    end
    
end

% hold off
% axis square
axis equal
% axis tight
xlabel('x')
ylabel('y')
zlabel('z')
view(3)

end
