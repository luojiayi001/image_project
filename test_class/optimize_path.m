function count = optimize_path(E_local, E_global, V_all, boundary_points, threshold, g) 

for j = 1: size(boundary_points, 2)
    start_point = boundary_points(j);
    for i = 1: size(boundary_points, 2)
        [c_l, p_l] = dijkstra(V_all,E_local, start_point, boundary_points(i));
        [c_g, p_g] = dijkstra(V_all,E_global, start_point, boundary_points(i));
        if abs(c_l - c_g) >= threshold && c_l ~= Inf
            count = count + 1;
            % change the local path's faces state back to unvisited
            for m = 1: size(p_l, 2)
                g.faces{m}.state = 0;
            end
            
            % change the global path's faces state to visited
            for m = 1: size(p_g, 2)
                g.faces{m}.state = 2;
            end
        end
    end
end