
function h_refine_curve(obj, knotsForInsertion1)

if ~iscolumn(knotsForInsertion1)
    knotsForInsertion1 = knotsForInsertion1';
end

if iscolumn(obj.uKnot)
    knots1 = obj.uKnot;
else
    knots1 = obj.uKnot';
end

p1 = obj.p;

nodes = nurb2proj( obj.controlPts, obj.weights);

[knots1_new, nodes_new] = refine_h_curve(knots1, nodes, p1, knotsForInsertion1);

[controlPoints, PTweights] = proj2nurbs(nodes_new);

obj.controlPts = controlPoints;
obj.weights = PTweights;

if iscolumn(obj.uKnot)
    obj.uKnot = knots1_new;
else
    obj.uKnot = knots1_new';
end

obj.update() % 进行对象更新

end

function [knots_new, nodes_new] = refine_h_curve(knots, nodes, p, knotsForInsertion)
    % Some useful constants
    constant_pp1 = p + 1;
    
    numKnots = size(knots, 1);
    numNodes = numKnots - constant_pp1;
    
    numKnotInsertions = size(knotsForInsertion, 1);
    numKnots_new = numKnots + numKnotInsertions;
%   numNodes_new = numNodes + numKnotInsertions;
    
    % Check for errors
    if (numNodes ~= size(nodes, 1))
        fprintf('Error: For the given knot vector, there should be %d nodes.\n\n', numNodes);
        
        knots_new = knots;
        nodes_new = nodes;
        
        return;
        
    elseif (numKnotInsertions == 0)
%         fprintf('Warning: Please specify at least one knot for insertion.\n\nNo operation has been performed.\n\n');
        
        knots_new = knots;
        nodes_new = nodes;
        
        return;
        
    end
    
    
    %----------------------------------------------------------------------
    %  Knot refinement for direction 1
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %  Initialize the new knot vector and nodes array
    %----------------------------------------------------------------------
    knots_new = [knots; zeros(numKnotInsertions, 1)];
    
    clear knots;
    
    nodes_new = [nodes; repmat(nodes(numNodes, :), numKnotInsertions, 1)];
    
    clear nodes;
    
    
    %----------------------------------------------------------------------
    %  Loop over the knots for insertion
    %----------------------------------------------------------------------
    % Search where the knot xi belongs at this index
    knotSearch_begin = constant_pp1;
    
    for b = 1 : numKnotInsertions
        % Knot that we will insert
        xi = knotsForInsertion(b);
        
        % Find where to insert the knot
        for k = knotSearch_begin : (numKnots_new - 1)
            if (knots_new(k) <= xi && xi < knots_new(k + 1))
                index_knot = k + 1;
                
                break;
            end
        end
        
        % Shift the "end" nodes
        for i = (numNodes + b - 1) : -1 : index_knot
            nodes_new(i, :) = nodes_new(i - 1, :);
        end
        
        % Update the "middle" nodes
        for i = (index_knot - 1) : -1 : (index_knot - p)
            alpha = (xi - knots_new(i)) / (knots_new(i + p) - knots_new(i));
            
            nodes_new(i, :) = alpha * nodes_new(i, :) + (1 - alpha) * nodes_new(i - 1, :);
        end
        
        % Shift the "end" knots
        for i = (numKnots + b) : -1 : (index_knot + 1)
            knots_new(i) = knots_new(i - 1);
        end
        
        
        %------------------------------------------------------------------
        %  Insert knot xi to the new knot vector
        %------------------------------------------------------------------
        knots_new(index_knot) = xi;
        
        
        %------------------------------------------------------------------
        %  Update index for the next iteration
        %------------------------------------------------------------------
        knotSearch_begin = index_knot;
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


