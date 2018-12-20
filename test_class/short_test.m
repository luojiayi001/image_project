    n = 7; V = 10*rand(n,2)
    I = delaunay(V(:,1),V(:,2));
    J = I(:,[2 3 1]); E = [I(:) J(:)]
    [costs,paths] = dijkstra(V,E)