function [C, g] = optimize_path(E_local, E_global, V_all, boundary_points, threshold, g) 
% count: how many local paths have been changed
    C = 0;

for j = 1: 1% size(boundary_points, 2)
    start_point = boundary_points(10);
    for i = 1: size(boundary_points, 2)
%         1
%         start_cor = V_all(start_point);
%         end_cor = V_all(boundary_points(i));
%         dist = (sum(start_cor - end_cor).^2)^0.5;
%         if dist <= 10
            [c_l, p_l] = dijkstra(V_all,E_local, start_point, boundary_points(i));
            [c_g, p_g] = dijkstra(V_all, E_global, start_point, boundary_points(i));
            if abs(c_l - c_g) >= threshold && c_l ~= Inf
                C = C + 1;
                % change the local path's faces state back to unvisited
                for m = 1: size(p_l, 2)
                    g.faces{p_l(m)}.state = 1;
                end

                % change the global path's faces state to visited
                for n = 1: size(p_g, 2)
                    g.faces{p_l(n)}.state = 2;
                end
            end
        end
    end
end