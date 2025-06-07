classdef NURBS1D < handle
    
    %   curve 1D
    
    properties 
        uKnot
        p
        controlPts
        weights
        eleNodeGlobal    % ��Ԫȫ�����Ӿ���
    end
    
    properties 
        numberElementsU
        numberElements
        noPtsU
        noPts
        C                % ��Ԫbezier��ȡ����
        elementNode      % ��Ԫ���Ӿ���
        elementVertex    % ��Ԫ�Ĳ�������
        elementLocation  % ��Ԫ��λ����Ϣ
    end
    
    methods 
        function obj = NURBS1D(uKnot,p,controlPts,weights)
            if nargin == 4
                %   �˴���ʾ��ϸ˵��
                obj.uKnot = uKnot;
                obj.p = p;
                obj.controlPts = controlPts;
                obj.weights = weights;
                obj.update()
            end
            
        end
    end
    
    methods % ϸ�� ���κ��� �κ������� �� ��ͼ�ĳ���
        
        function update(obj)
            % ϸ�ֺ�����֮�� Ҳ��Ҫ���и���
            obj.numberElementsU = length(unique(obj.uKnot))-1;
            obj.numberElements = obj.numberElementsU;
            obj.noPtsU = length(obj.uKnot) - obj.p - 1;
            obj.noPts = obj.noPtsU;
            
            tempC = zeros((obj.p+1),(obj.p+1),obj.numberElements);
            [C_u, ~] = bezierExtraction(obj.uKnot,obj.p);
            
            % ���㵥ԪBezier��ȡ����
            elementCounter = 0;
            for i=1:obj.numberElementsU
                elementCounter = elementCounter + 1;
                tempC(:,:,elementCounter) = C_u(:,:,i);
            end
            obj.C = tempC;
            
            elementCounter = 0;
            temp_elementVertex = zeros(obj.numberElements,2);
            temp_elementNode = zeros(obj.numberElements, (obj.p+1));
            toleq = 1e-10;
            
            for i=1:length(obj.uKnot)-1
                if abs(obj.uKnot(i+1)-obj.uKnot(i))>toleq
                    elementCounter = elementCounter + 1;
                    temp_elementVertex(elementCounter, :) = [obj.uKnot(i), obj.uKnot(i+1)];
                    % �� i-p...i ��u����
                    temp_elementNode(elementCounter,:)=i-obj.p:i;
                end
            end

            obj.elementVertex = temp_elementVertex;
            obj.elementNode = temp_elementNode;
            obj.eleNodeGlobal = temp_elementNode; % ��ʵʱ�̸�Ƭȫ�����Ӿ�����ֲ����Ӿ���һ��
            
            % ����ÿ����Ԫ�ڽ���Ԫ�б�
            neighbor = zeros(obj.numberElements,2);
            for i=1:obj.numberElementsU
                
                if i > 1
                    neighbor(i,1) = i-1; % neighbor_left
                end
                
                if i < obj.numberElementsU
                    neighbor(i,2) = i+1; % neighbor_right
                end
                
            end
            
            obj.elementLocation = neighbor;

        end
        
        p_refine_curve(obj, padd)
        
        h_refine_curve(obj, knotsForInsertion1)
        
        plotMeshCurve(obj,s1)
                
        R = Shapes(obj,xiG,ele)
        
        dRdxi = Shape1stDers(obj,xiG,ele)
                
        [dRdxi,J1] = Shape1stDersAndDRdx(obj,xiG,ele)

        [R, dRdxi] = ShapesAnd1stDers(obj,xiG,ele)

        PhyCoords = GetPhyCoords(obj,xiP)
        
        DxDxiG = GetPhyToGaussDers(obj,xiG,ele)

    end
    
end

function [C, nb] = bezierExtraction(knot,p)
% Bezier extraction

m  = length(knot)-p-1;
a  = p+1;
b  = a+1;
nb = 1;
C(:,:,1) = eye(p+1);

while b <= m
    C(:,:,nb+1) = eye(p+1);
    i=b;
    while b <= m && knot(b+1) == knot(b)
        b=b+1;
    end
    
    multiplicity = b-i+1;
    if multiplicity < p
        numerator=knot(b)-knot(a);
        for j=p:-1:multiplicity+1
            alphas(j-multiplicity)=numerator/(knot(a+j)-knot(a));
        end
        r=p-multiplicity;
        for j=1:r
            save = r-j+1;
            s = multiplicity + j;
            for k=p+1:-1:s+1
                alpha=alphas(k-s);
                C(:,k,nb)=alpha*C(:,k,nb)+(1-alpha)*C(:,k-1,nb);
            end
            if b <= m
                C(save:save+j,save,nb+1)=C(p-j+1:p+1,p+1,nb);
            end
        end
        nb=nb+1;
        if b <= m
            a=b;
            b=b+1;
        end
    elseif multiplicity==p
        if b <= m
            nb=nb+1; a=b; b=b+1;
        end
    end
end
end

