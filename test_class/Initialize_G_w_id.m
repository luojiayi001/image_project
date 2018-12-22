% This function is aimed to create a subgraph with the same id
function g = Initialize_G_w_id(vertex, vid, faces, fid, state)
% Denote V = No. of vertex, F = No. of faces
% vertex: 3 x V matrix for coordinates (x, y, z) of V vertices
% faces:  3 x F matrix for vertex index (I1, I2, I3) of F faces
face_rings = compute_face_ring(faces); % face_rings{i}: a list of adjacent faces of ith face 
face_normals = zeros(size(faces));     % face_normals(:, i): 3 x 1 face normal vector of ith face
face_centers = zeros(size(faces));     % face_centers(:, i): 3 x 1 centroid of ith face
% Compute face normals and face centers
for i = 1 : size(faces, 2)
    face = faces(:, i);
    v1 = vertex(:, face(1));
    v2 = vertex(:, face(2));
    v3 = vertex(:, face(3));
    normal = cross(v2 - v1, v3 - v1); % assume a counter-clockwise order convention of vertices, then the cross product should also be counter-clockwise, v2-v1-->v3-v1, right-hand principle
    normal = normal / norm(normal); % normalize vector
    face_normals(:, i) = normal;
    face_centers(:, i) = (v1 + v2 + v3) / 3;
end
% Create a Graph class
g = Graph;
for i = 1:size(vertex, 2)
    temp = Node;
    temp.state = 0; % Initialize it as new node 
%     temp.id = i;
    temp.id = vid(i);
    temp.coor = vertex(:,i)';
    g.nodes{end + 1} = temp;
end

for i = 1:size(faces, 2)
    temp = Face;
    temp.state = state; % Initialize it as new face
%     temp.id = i;
    temp.id = fid(i);
    temp.neighbors = face_rings{i};
    temp.norm = face_normals(:,i);
    temp.center = face_centers(:,i);
    temp.nodes = faces(:,i)';
    g.faces{end + 1} = temp;
    
    % Add ajacent faces to node class
    for j = 1:3
        g.nodes{faces(j,i)}.faces(end + 1) = i;
    end
end

 % Add adj_faces to the face obj
 for i = 1:size(faces, 2)
     temp = [];
     for k = 1:3
         vertex_id = g.faces{i}.nodes(k);
         temp = union(temp, g.nodes{vertex_id}.faces);
     end
     g.faces{i}.all_neighbors = setdiff(temp, i);
%      g.faces{i}.all_neighbors = temp;
     g.faces{i}.ajac_faces = setdiff(temp, [i, g.faces{i}.neighbors]);
 end
end 