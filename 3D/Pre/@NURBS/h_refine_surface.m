
function h_refine_surface(obj, knotsForInsertion1, knotsForInsertion2)

if ~iscolumn(knotsForInsertion1)
    knotsForInsertion1 = knotsForInsertion1';
end

if ~iscolumn(knotsForInsertion2)
    knotsForInsertion2 = knotsForInsertion2';
end

if iscolumn(obj.uKnot)
    knots1 = obj.uKnot;
else
    knots1 = obj.uKnot';
end

if iscolumn(obj.vKnot)
    knots2 = obj.vKnot;
else
    knots2 = obj.vKnot';
end

p1 = obj.p;
p2 = obj.q;

nodes = nurb2proj( obj.controlPts, obj.weights);

[knots1_new, knots2_new, nodes_new] = refine_h_surface...
    (knots1, knots2, nodes, p1, p2, knotsForInsertion1, knotsForInsertion2);

[controlPoints, PTweights] = proj2nurbs(nodes_new);

obj.controlPts = controlPoints;
obj.weights = PTweights;

if iscolumn(obj.uKnot)
    obj.uKnot = knots1_new;
else
    obj.uKnot = knots1_new';
end

if iscolumn(obj.vKnot)
    obj.vKnot = knots2_new;
else
    obj.vKnot = knots2_new';
end

obj.update() % 进行对象更新

end

function [knots1_new, knots2_new, nodes_new] = refine_h_surface(knots1, knots2, nodes, p1, p2, knotsForInsertion1, knotsForInsertion2)
    % Some useful constants
    constant_p1p1 = p1 + 1;
    constant_p2p1 = p2 + 1;
    
    numKnots1 = size(knots1, 1);
    numKnots2 = size(knots2, 1);
    numNodes1 = numKnots1 - constant_p1p1;
    numNodes2 = numKnots2 - constant_p2p1;
    numNodes  = numNodes1 * numNodes2;
    numDimensions = size(nodes, 2);
    
    numKnotInsertions1 = size(knotsForInsertion1, 1);
    numKnotInsertions2 = size(knotsForInsertion2, 1);
    numKnots1_new = numKnots1 + numKnotInsertions1;
    numKnots2_new = numKnots2 + numKnotInsertions2;
    numNodes1_new = numNodes1 + numKnotInsertions1;
    numNodes2_new = numNodes2 + numKnotInsertions2;
    
    % Check for errors
    if (numNodes ~= size(nodes, 1))
        fprintf('Error: For the given knot vector, there should be %d nodes.\n\n', numNodes);
        
        knots1_new = knots1;
        knots2_new = knots2;
        nodes_new = nodes;
        
        return;
        
    elseif (numKnotInsertions1 == 0 && numKnotInsertions2 == 0)
%         fprintf('Warning: Please specify at least one knot for insertion.\n\nNo operation has been performed.\n\n');
        
        knots1_new = knots1;
        knots2_new = knots2;
        nodes_new = nodes;
        
        return;
        
    end

    %  Knot refinement for direction 1
    if (numKnotInsertions1 > 0)
        %  Shifts to find all the nodes on the same "line"
        temp2 = (0 : (numNodes2 - 1))';
        
        shifts_for_nodes     = numNodes1     * temp2;
        shifts_for_nodes_new = numNodes1_new * temp2;
        
        clear temp2;

        %  Initialize the new knot vector and nodes array
        knots1_new = [knots1; zeros(numKnotInsertions1, 1)];
        
        clear knots1;
        
        nodes_new = zeros(numNodes1_new * numNodes2, numDimensions);
        
        for i = 1 : numNodes1
            % For the last node along direction 1, we make extra copies
            if (i < numNodes1)
                nodes_new(i + shifts_for_nodes_new, :) = nodes(i + shifts_for_nodes, :);
            else
                for j = 0 : numKnotInsertions1
                    nodes_new(i + j + shifts_for_nodes_new, :) = nodes(i + shifts_for_nodes, :);
                end
            end
        end
        
        clear nodes shifts_for_nodes;

        %  Loop over the knots for insertion
        % Search where the knot xi belongs at this index
        knotSearch_begin = constant_p1p1;
        
        for b = 1 : numKnotInsertions1
            % Knot that we will insert
            xi = knotsForInsertion1(b);
            
            % Find where to insert the knot
            for k = knotSearch_begin : (numKnots1_new - 1)
                if (knots1_new(k) <= xi && xi < knots1_new(k + 1))
                    index_knot = k + 1;
                    
                    break;
                end
            end
            
            % Shift the "end" nodes
            for i = (numNodes1 + b - 1) : -1 : index_knot
                nodes_new(i + shifts_for_nodes_new, :) = nodes_new(i - 1 + shifts_for_nodes_new, :);
            end
            
            % Update the "middle" nodes
            for i = (index_knot - 1) : -1 : (index_knot - p1)
                alpha = (xi - knots1_new(i)) / (knots1_new(i + p1) - knots1_new(i));
                
                nodes_new(i + shifts_for_nodes_new, :) = alpha  * nodes_new(i + shifts_for_nodes_new, :) + (1 - alpha) * nodes_new(i - 1 + shifts_for_nodes_new, :);
            end 
            
            % Shift the "end" knots
            for i = (numKnots1 + b) : -1 : (index_knot + 1)
                knots1_new(i) = knots1_new(i - 1);
            end

            %  Insert knot xi to the new knot vector
            knots1_new(index_knot) = xi;

            %  Update index for the next iteration
            knotSearch_begin = index_knot;
        end
        
        nodes = nodes_new;
        
        clear nodes_new;
        
    % Default action for no refinement
    else
        knots1_new = knots1;
        
        clear knots1;
    end

    %  Knot refinement for direction 2
    if (numKnotInsertions2 > 0)
        %------------------------------------------------------------------
        %  Shifts to find all the nodes on the same "line"
        %  
        %  Note that numNodes1 = numNodes1_new now. We add an extra shift
        %  of -(numNodes1_new - 1) so that we can sequentially traverse
        %  along direction 2 like we did along direction 1.
        %------------------------------------------------------------------
        shifts_for_nodes     = (0 : (numNodes1_new - 1))' - (numNodes1_new - 1);
        shifts_for_nodes_new = shifts_for_nodes;

        %  Initialize the new knot vector and nodes array
        knots2_new = [knots2; zeros(numKnotInsertions2, 1)];
        
        clear knots2;
        
        nodes_new = zeros(numNodes1_new * numNodes2_new, numDimensions);
        
        for i = 1 : numNodes2
            % For the last node along direction 1, we make extra copies
            if (i < numNodes2)
                nodes_new(numNodes1_new * i + shifts_for_nodes_new, :) = nodes(numNodes1_new * i + shifts_for_nodes, :);
            else
                for j = 0 : numKnotInsertions2
                    nodes_new(numNodes1_new * (i + j) + shifts_for_nodes_new, :) = nodes(numNodes1_new * i + shifts_for_nodes, :);
                end
            end
        end
        
        clear nodes shifts_for_nodes;

        %  Loop over the knots for insertion
        % Search where the knot xi belongs at this index
        knotSearch_begin = constant_p2p1;
        
        for b = 1 : numKnotInsertions2
            % Knot that we will insert
            xi = knotsForInsertion2(b);
            
            % Find where to insert the knot
            for k = knotSearch_begin : (numKnots2_new - 1)
                if (knots2_new(k) <= xi && xi < knots2_new(k + 1))
                    index_knot = k + 1;
                    break;
                end
            end
            
            % Shift the "end" nodes
            for i = (numNodes2 + b - 1) : -1 : index_knot
                nodes_new(numNodes1_new * i + shifts_for_nodes_new, :) = nodes_new(numNodes1_new * (i - 1) + shifts_for_nodes_new, :);
            end
            
            % Update the "middle" nodes
            for i = (index_knot - 1) : -1 : (index_knot - p2)
                alpha = (xi - knots2_new(i)) / (knots2_new(i + p2) - knots2_new(i));
                
                nodes_new(numNodes1_new * i + shifts_for_nodes_new, :) = alpha * nodes_new(numNodes1_new * i + shifts_for_nodes_new, :) + (1 - alpha) * nodes_new(numNodes1_new * (i - 1) + shifts_for_nodes_new, :);
            end
            
            % Shift the "end" knots
            for i = (numKnots2 + b) : -1 : (index_knot + 1)
                knots2_new(i) = knots2_new(i - 1);
            end

            %  Insert knot xi to the new knot vector
            knots2_new(index_knot) = xi;

            %  Update index for the next iteration
            knotSearch_begin = index_knot;
        end
        
    % Default action for no refinement
    else
        knots2_new = knots2;
        nodes_new = nodes;
        
        clear knots2 nodes;
    end
end

function projcoord = nurb2proj( controlPoints, PTweights)
%--------------------------------------------------------------
% function projcoord = nurb2proj(nob, controlPoints, weights)
% transform NURBS data into projective coordinates
% INPUT:
% controlPoints: vector of control points (1 per row)
% PTweights :    : column vector of weights
% OUTPUT:
% projcoord    : matrix with projective coordinates
%--------------------------------------------------------------
projcoord = controlPoints;
for i=1:size(controlPoints,1)
    projcoord(i,:) = projcoord(i,:)*PTweights(i);
end
projcoord = [projcoord, PTweights];
end

function [controlPoints, PTweights] = proj2nurbs(projcoord)
%--------------------------------------------------------------
% function [controlPoints, weightVector] = proj2nurbs(projcoord)
% transform projective coordinates into NURBS data
% INPUT:
% projcoord    : matrix with projective coordinates
% OUTPUT:
% controlPoints: vector of control points (1 per row)
% PTweightsVector : column vector of weights
%--------------------------------------------------------------
dimension     = size(projcoord,2);
PTweights       = projcoord(:,dimension);
controlPoints = projcoord(:,1:dimension-1);

for i=1:size(PTweights,1)
    controlPoints(i,:) = controlPoints(i,:)* 1/(PTweights(i));
end
end


