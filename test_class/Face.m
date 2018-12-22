% ================================================================== %
classdef Face
    properties
        id
        state % 0: not visited; 1: visited; 2: boundary
        norm % 3D norm vector
        nodes % surrounded nodes id(# = 3)
        neighbors % surrounded faces id
        center % coordinates of the center point
        ajac_faces % those faces that share only one node with this face
        all_neighbors
    end
end
% ================================================================== %