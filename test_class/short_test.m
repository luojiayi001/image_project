start_point = boundary_points(5);
for i = 1: size(boundary_points, 2)
    [c_l, p_l] = dijkstra(V_all,E_local, start_point, boundary_points(i));
    [c_g, p_g] = dijkstra(V_all,E_global, start_point, boundary_points(i));
    if abs(c_l - c_g) >= 2 && c_l ~= Inf
        p_l
        p_g
    end
end