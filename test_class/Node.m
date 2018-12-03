% ================================================================== %
classdef Node
    properties
        id 
        state % 0: not visited; 1: visited; 2: boundary
        faces % Face id list that shared this node
        coor % 3D coordinates
    end
end
% ================================================================== %