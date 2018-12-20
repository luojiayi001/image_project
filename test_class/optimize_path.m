
% bound_E: n by 2 matrix corresponds to the faces relationship within the boundary
% global_E: n by 2 matrix corresponds to the faces relationship within the global domain
% V: the coordinates for all the nodes
% threshold: Use subtraction for threshold (absolute diff in distance)
% In order to avoid local optimization, two points of which the distance between should be larger than threshold
function [point_1, point_2] = optimize_path(bound_E, global_E, V, face_rings, threshold, num_check) 
% Overall structure
point_1 = 0;
point_2 = 0; % Let algorithm search within the visited path
iter = 0;
while (point_1 == 0 || point_2 == 0) && iter < num_check % Another point futher enough will be elected to be the new start point
    visited = zeros(size(bound_E, 1), 1);
    iter = iter + 1;
    point_1 = get_corner_point(bound_E, global_E, V, face_rings, threshold, visited, start);
    point_2 = get_corner_point(bound_E, global_E, V, face_rings, threshold, 1 - visited, start);
    new_start = randi(size(bound_E, 1)); 
    while cal_dist(new_start, start) < threshold
        new_start = randi(size(bound_E, 1)); 
    end
    start = new_start;
end
end

function [point_id] = get_corner_point(bound_E, global_E, V, face_rings, threshold, visited, start)
point_id = 0;
q = CQueue(); 
q.push(start);
while(q.isempty() ~= 1)
    dest = q.pop();
    [cost,g_path] = dijkstra(V,global_E,[start],[dest]);
    [cost,b_path] = dijkstra(V,bound_E,[start],[dest]);
    dist_g = cal_dist(g_path);
    dist_b = cal_dist(b_path);
    if dist_b - dist_g > threshold
        point_id = dest;
        break;
    else
        neighbors = face_rings(dest);
        for j = 1:length(neighbors)
            visited(neighbors(j)) = 1;
            q.push(neighbors(j));
        end
    end      
end
end

function dist = cal_dist(path)
cur = path(1);
dist = 0;
for i = 1:length(path)
    dist = dist + (sum(cur - path(i)).^2)^0.5;
end
end