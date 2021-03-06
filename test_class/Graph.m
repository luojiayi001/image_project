% ================================================================== %
classdef Graph
    properties
        nodes 
        faces
    end
    methods
        % Count the total number of visited faces
        function r = Sum_visited_faces(obj)
            r = 0;
            for i = 1 : size(obj.faces, 2)
                if obj.faces{i}.state > 0
                    r = r + 1;
                end
            end
        end
        % Create the visited faces vector
        function r = Visited_faces(obj)
            r = [];
            for i = 1 : size(obj.faces, 2)
                if obj.faces{i}.state ~= 0
                    r(i) = 1;
                    
                else
                    r(i) = 0;
                end
            r = r';
            end
        end
        % Create the boundary faces vector
        function r = Boundary_faces(obj)
            r = [];
            for i = 1 : size(obj.faces, 2)
                if obj.faces{i}.state == 2
                    r(i) = 1;
                else
                    r(i) = 0;
                end
            r = r';
            end
        end
        % Create the boundary faces vector
        function vertex = Extract_v(obj)
            vertex = [];
            for i = 1:size(obj.nodes, 2)
                 vertex = [vertex; obj.nodes{i}.coor];
            end
            vertex = vertex';
        end
        function face = Extract_f(obj)
            face = [];
            for i = 1:size(obj.faces, 2)
                face = [face; obj.faces{i}.nodes];
            end
            face = face';
        end

        % Create the boundary faces vector
        function idx = Error_detect(obj)
            idx = 0;
            for i = 1 : size(obj.faces, 2)
                if obj.faces{i}.state == 1
                    for j = 1 : size(obj.faces{i}.neighbors)
                        state_temp = obj.faces{obj.faces{i}.neighbors(j)}.state;
                        if state_temp == 0
                            idx = idx + 1;
                        end
                    end
                end
            end
        end
        
        % Extract the vertex and faces from a graph
        function new_g = Get_small(obj, state)
            new_g = Graph;
            face_ref = [];
            node_ref = [];
            for i = 1 : size(obj.faces, 2)
                idx = 0;
                if obj.faces{i}.state == state%1
                    face_ref(end + 1) = obj.faces{i}.id;
                    for j = 1 : 3
                        if (size(find(node_ref == obj.faces{i}.nodes(j)), 2) == 0)
                            node_ref(end + 1) = obj.faces{i}.nodes(j);
                        end
                    end
                end
            end
            size(node_ref);
            size(face_ref);
            for i = 1 : size(node_ref, 2)
                temp = Node;
%                 temp.id = i;
                temp.id = obj.nodes{node_ref(i)}.id;
                temp.coor = obj.nodes{node_ref(i)}.coor;
                for n = 1 : size(obj.nodes{node_ref(i)}.faces, 2)
                    if size(find(face_ref == obj.nodes{node_ref(i)}.faces(n))) == 1
                        temp.faces(end + 1) = find(face_ref == obj.nodes{node_ref(i)}.faces(n));
                    end
                end
                new_g.nodes{end + 1} = temp;
            end
            idx = 0;
            for i = 1 : size(obj.faces, 2)
                if obj.faces{i}.state == state%1
                    idx = idx + 1;
                    temp = Face;
%                     temp.id = idx;
                    temp.id = obj.faces{i}.id;
%                     temp.state = 1;
                    temp.state = state;
                    temp.norm = obj.faces{i}.norm;
                    temp.center = obj.faces{i}.center;
                    for j = 1 : size(obj.faces{i}.nodes, 2)
                        temp.nodes(j) = find(node_ref == obj.faces{i}.nodes(j));
                    end
                    for m = 1 : size(obj.faces{i}.neighbors, 2)
                        if size(find(face_ref == obj.faces{i}.neighbors(m)), 2) == 1
                            temp.neighbors(end + 1) = find(face_ref == obj.faces{i}.neighbors(m));
                        end
                    end
                    new_g.faces{end + 1} = temp;
                end
            end
        end
    end
end
% ================================================================== %