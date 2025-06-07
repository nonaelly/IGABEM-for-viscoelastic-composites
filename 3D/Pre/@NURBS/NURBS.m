classdef NURBS < handle
    
    properties 
        uKnot
        vKnot
        p
        q
        controlPts
        weights
        eleNodeGlobal     % ��Ԫȫ�����Ӿ���
        controlPts_update % Ϊ��Ҫ���п��Ƶ��غ�Ԥ�������Ա���
        weights_update    % Ϊ��Ҫ���п��Ƶ��غ�Ԥ�������Ա���
    end
    
    properties 
        numberElementsU
        numberElementsV
        numberElements
        noPtsU
        noPtsV
        noPts
        C                % ��Ԫbezier��ȡ����
        elementNode      % ��Ԫ���Ӿ���
        elementNode_update % ĳЩ��β��ӵļ���ʹ�õ�
        elementVertex    % ��Ԫ�Ĳ�������
        elementLocation  % ��Ԫ��λ����Ϣ
    end
    
    properties (Hidden)
        down_nodes
        right_nodes
        up_nodes
        left_nodes
    end
    
    methods % ��handle �̳� �� get ���� ���� Dependent ���Ե�����ֵ
        
        function obj = NURBS(uKnot,vKnot,p,q,controlPts,weights)
            if nargin == 6
                %   �˴���ʾ��ϸ˵��
                obj.uKnot = uKnot;
                obj.vKnot = vKnot;
                obj.p = p;
                obj.q = q;
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
            obj.numberElementsV = length(unique(obj.vKnot))-1;
            obj.numberElements = obj.numberElementsU*obj.numberElementsV;
            obj.noPtsU = length(obj.uKnot) - obj.p - 1;
            obj.noPtsV = length(obj.vKnot) - obj.q - 1;
            obj.noPts = obj.noPtsU * obj.noPtsV;

            % ���㵥ԪBezier��ȡ����
            [C_u, ~] = bezierExtraction(obj.uKnot,obj.p);
            [C_v, ~] = bezierExtraction(obj.vKnot,obj.q);
            tempC = zeros((obj.p+1)*(obj.q+1),(obj.p+1)*(obj.q+1),obj.numberElements);
            elementCounter = 0;
            for j=1:obj.numberElementsV
                for i=1:obj.numberElementsU
                    elementCounter = elementCounter + 1;
                    tempC(:,:,elementCounter) =  kron(C_v(:,:,j),C_u(:,:,i));
                end
            end
            obj.C = tempC;
            
            elementCounter = 0;
            temp_elementVertex = zeros(obj.numberElements,4);
            temp_elementNode = zeros(obj.numberElements, (obj.p+1)*(obj.q+1));
            toleq = 1e-10;
            for j=1:length(obj.vKnot)-1
                for i=1:length(obj.uKnot)-1
                    if (abs(obj.uKnot(i+1)-obj.uKnot(i))>toleq) && (abs(obj.vKnot(j+1)-obj.vKnot(j))>toleq)
                        elementCounter = elementCounter + 1;
                        temp_elementVertex(elementCounter, :) = [obj.uKnot(i), obj.vKnot(j), obj.uKnot(i+1), obj.vKnot(j+1)];
                        
                        tcount = 0;
                        currow = zeros(1, (obj.p+1)*(obj.q+1));
                        % �� i-p...i ��u���򣬴� j-q...j ��v����
                        for t2=j-obj.q:j
                            for t1 = i-obj.p:i
                                tcount = tcount + 1;
                                currow(tcount) = t1+(t2-1)*obj.noPtsU;
                            end
                        end
                        temp_elementNode(elementCounter,:)=currow;
                    end
                end
            end
            obj.elementVertex = temp_elementVertex;
            obj.elementNode = temp_elementNode;
            obj.eleNodeGlobal = temp_elementNode; % ��ʵʱ�̸�Ƭȫ�����Ӿ�����ֲ����Ӿ���һ��
            
            % ����ÿ����Ԫ�ڽ���Ԫ�б�
            neighbor = zeros(obj.numberElements,4);
            indexMatrix = permute(reshape(1:obj.numberElements, obj.numberElementsU, obj.numberElementsV),[2,1]);
            % permute(��ά����,[ά�������])
            % ��ά���齻��ά�����൱��ʵ���˾���ת��
            for j=1:obj.numberElementsV
                for i=1:obj.numberElementsU
                    elementIndex = indexMatrix(j,i);
                    
                    if i > 1
                        neighbor(elementIndex,4) = indexMatrix(j,i-1); % neighbor_left
                    end
                    
                    if i < obj.numberElementsU
                        neighbor(elementIndex,2) = indexMatrix(j,i+1);  % neighbor_right
                    end
                    
                    if j > 1
                        neighbor(elementIndex,1) = indexMatrix(j-1,i);  % neighbor_down
                    end
                    
                    if j < obj.numberElementsV
                        neighbor(elementIndex,3) = indexMatrix(j+1,i); % neighbor_up
                    end
                    
                end
            end
            
            obj.elementLocation = neighbor;
            
            obj.down_nodes = 1:obj.p+1;
            obj.left_nodes = 1:(obj.p+1):(1+(obj.p+1)*obj.q);
            obj.right_nodes = (obj.p+1):(obj.p+1):(obj.p+1)*(obj.q+1);
            obj.up_nodes = 1+(obj.p+1)*obj.q:(obj.p+1)*(obj.q+1);
            
%             log = LogClass.getInstance();
%             log.print('function update is called!\n');
        end
        
        p_refine_surface(obj,padd,qadd)
        
        h_refine_surface(obj, knotsForInsertion1, knotsForInsertion2)
        
        plotMeshSurface(obj,s1)
        
        R = Shapes(obj,xiG,etaG,ele)
        
        [dRdxi, dRdeta] = Shape1stDers(obj,xiG,etaG,ele)
        
        [R, dRdxi, dRdeta] = ShapesAnd1stDers(obj,xiG,etaG,ele)
        
        [R, dRdxi, dRdeta, dR2dxi2, dR2deta2, dR2dxieta] = Shapes1st2ndDers(obj,xiG,etaG,ele)
        
        PhyCoords = GetPhyCoords(obj,xiP,etaP)
        
        [nodes, elements] = sortNodesAndElemsOfEdge(obj,noEdge,FlagDir)
                
        function indexCPts = ExtractEdgeCPts(obj,noEdge)
            
            switch noEdge
                
                case 1
                    indexCPts = 1:obj.noPtsU;
                    
                case 4
                    indexCPts = 1:obj.noPtsU:obj.noPtsU*(obj.noPtsV-1)+1;
                
                case 2
                    indexCPts = obj.noPtsU:obj.noPtsU:obj.noPtsU*obj.noPtsV;
                
                case 3
                    indexCPts = obj.noPtsU*(obj.noPtsV-1)+1:obj.noPtsU*obj.noPtsV;
                
            end
            
        end
        
        function indexCPts = ExtractNextEdgeCPts(obj,noEdge)

            switch noEdge
                
                case 1
                    indexCPts = obj.noPtsU+1:2*obj.noPtsU;
                    
                case 4
                    indexCPts = 2:obj.noPtsU:obj.noPtsU*(obj.noPtsV-1)+2;
                
                case 2
                    indexCPts = (obj.noPtsU-1):obj.noPtsU:(obj.noPtsU*obj.noPtsV-1);
                
                case 3
                    indexCPts = obj.noPtsU*(obj.noPtsV-2)+1:obj.noPtsU*(obj.noPtsV-1);
                
            end
            
        end
        
        function [Uvec,Porder,CPts,Wgts] = ConstructEdge(obj,noEdge)
            
            indexCPts = obj.ExtractEdgeCPts(noEdge);
            CPts = obj.controlPts(indexCPts,:);
            Wgts = obj.weights(indexCPts);
            
            switch noEdge
                case {1,3}
                    Uvec = obj.uKnot;
                    Porder = obj.p;
                   
                case {2,4}
                    Uvec = obj.vKnot;
                    Porder = obj.q;
            end
        end
        
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
