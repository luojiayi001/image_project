% Libraries:
% Graph & Mesh: https://www.mathworks.com/matlabcentral/fileexchange/5355-toolbox-graph 
% Queue/Stack (queue line96 bug fix): https://www.mathworks.com/matlabcentral/fileexchange/28922-list-queue-stack
addpath(genpath('toolbox_graph'));
addpath(genpath('datastructure'));
clear options;

clc;
close all;
clear;
name = 'all_particle'; %'single_particle';
options.name = name; % useful for displaying
[vertex,faces] = read_mesh(name);
% A = triangulation2adjacency(faces, vertex);

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
    temp.id = i;
    temp.coor = vertex(:,i)';
    g.nodes{end + 1} = temp;
end

for i = 1:size(faces, 2)
    temp = Face;
    temp.state = 0; % Initialize it as new face
    temp.id = i;
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

% Breadth-First Search
q = CQueue();                          % queue for BFS, storing the face index, i.e. ith face
%visited = zeros(size(faces, 2), 1);    % array for recording visited nodes. 1-visited, 0-unvisited
%boundary = zeros(size(faces, 2), 1);   % array for labelling boundary faces. 1-boundary face, 0-non-boundary face
threshold1 = 0.6;                      % local criteria for convexity/curvature (adjustable)
threshold2 = 0.001;                      % global criteria for centroid facing (adjustable)

seed = 2388; % seed = randi(size(faces, 2)); % or push a random starting face index
q.push(seed); center_sum = zeros(3, 1); object_count = 0;
g.faces{seed}.state = 1;
% stop = 0; sign = 11404;
s = 1; % Track how many visited vertices
while q.isempty() ~= 1
%     stop = stop + 1;
%     if stop == sign
%         break;
%     end
    curr_id = q.pop();
    curr_normal = - g.faces{curr_id}.norm; % flip sign, now normal points to the interior of a particle
    curr_center = g.faces{curr_id}.center;
    
    % Update the centroid
    center_sum = center_sum + curr_center;
    object_count = object_count + 1;
    centroid = center_sum / object_count;
    
    % Check if a neighboring face is object OR boundary face
    neighbors = face_rings{curr_id};
    %neighbors = g.faces{curr_id}.neighbors;
    for f = 1 : length(neighbors)
        adj_id = neighbors(f);
        adj_normal = - g.faces{adj_id}.norm; % flip sign, now normal points to the interior of a particle
        adj_center = g.faces{adj_id}.center;
        if g.faces{adj_id}.state == 0
            g.faces{adj_id}.state = 1;
            s = s + 1;
            % Criteria 1: Local criteria (between a face and its adjacent faces)
            local_center_diff = adj_center - curr_center;
            local_center_diff = local_center_diff / norm(local_center_diff); % normalize
            local_measure = dot(curr_normal, local_center_diff + adj_normal);
            % Criteria 2: Global criteria (between a face and the global centroid)
            global_center_diff = adj_center - centroid;
            % global_measure = dot(adj_normal, global_center_diff); % Attempt 1: dot product
            global_center_diff = global_center_diff / norm(global_center_diff);
            global_measure = dot(-adj_normal, global_center_diff); % Attempt 2: projection (global_center_diff is normalized, so |b|*cos(theta) = projection)
            if s < 11404/10
                test_t = -0.5;
            else
                test_t = 0.3;
            end
            if local_measure > threshold1 && global_measure > test_t % global_measure < threshold2
                % object face (push to queue)
                q.push(adj_id);
            else
                % boundary face (not push to queue)
                % label its adjacency as well
                boundary_neighbors = g.faces{adj_id}.neighbors;
                for b = 1 : length(boundary_neighbors)
                    boundary_id = boundary_neighbors(b);
                    if g.faces{boundary_id}.state == 0 
                        g.faces{boundary_id}.state = 2;
                        s = s + 1;
                    end
                end
            end
        end
    end
end

% Display the mesh
% Per-face coloring: Non-visited face-0 (black), object vertex-0.3 (red), boundary vertex-1.0 (white)
visited = g.Visited_faces(); % Create the visited vector for plotting
boundary = g.Boundary_faces(); % Create the boundary vector for plotting
objects = visited == 1 & boundary == 0;
boundaries = boundary == 1;
face_colors = zeros(size(faces, 2), 1);
face_colors(objects) = 0.3;
face_colors(boundaries) = 1.0;

new_g = g.Filter_particle();
visited1 = new_g.Visited_faces()';
boundary1 = new_g.Boundary_faces()';
new_faces = zeros(3 , 1);
for i = 1 : size(new_g.faces)
    new_faces(end + 1) = new_g.faces{i}.nodes;
DISPLAY_HHH = false;
if DISPLAY_HHH
    % My own display settings
    clf;
    h = patch('vertices',vertex','faces',faces','FaceVertexCData', face_colors, 'FaceColor', 'flat');
    lighting phong;
    camproj('perspective');
    axis square; 
    axis off;
    cameramenu;
    axis tight;
    axis equal;
    shading interp;
    camlight;
    % camzoom(1.0);
else
    % Use toolbox display settings
    options.face_vertex_color = face_colors;
    clf;
    plot_mesh(vertex, faces, options);
end
hold on;
% id = 2388; % plot normal for a specific face
% quiver3(face_centers(1,id), face_centers(2,id), face_centers(3,id), face_normals(1,id), face_normals(2,id), face_normals(3,id),  'r'); % draw normal
shading faceted; % display edges % shading interp; % display smooth surface 
colormap hot;  % colormap jet(256);
caxis([0 1]); % fix colormap range
% particle = Filter_particle(g);
% Per-vertex coloring (not recommended)
% objects = faces(:, visited == 1 & boundary == 0);
% boundaries = faces(:, boundary == 1);
% objects = unique(objects); % per-vertex coloring, so we should give it a list of vertex. 
% boundaries = unique(boundaries);
% colors = zeros(size(vertex, 2), 1);
% colors(objects) = 0.5;
% colors(boundaries) = 1.0;
% options.face_vertex_color = colors;
% clf;
% plot_mesh(vertex, faces, options);
% shading interp; 
% colormap jet(256);
% caxis([0 1]); % fix colormap range

% Consider first calculate the principal direction of all the visited
% faces, and project the current faces onto a plane along the principal
% direction. Then, for the uncertain faces, we do the same projection and
% check if its normal direction points towards or away from the 2D
% centroid, and decide whether this face is a whole on the surface, or the
% curvature between different particles. This could be a solution for
% identifying if a face points to the "inside" of the body.

% Thinking:
% 1. the center difference vector should be normalized. (otherwise if the 
% magnitude of center difference vector is much smaller than the normalized 
% normal vector, the normal vector will dominate the resultant direction)
% 2. the threshold of the dot product should be further investigated in
% terms of curvature limit, looks like 0.5~0.6 an promising range
% 3. currently the boundary does not close, maybe more restrictions should
% be applied to close the boundary: (1) use the dynamically updated "trial 
% center" -- reliable but expensive (2) label adjacent faces of a boundary
% face as boundary faces as well -- cheap but vulnerable