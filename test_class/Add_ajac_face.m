function ret = Add_ajac_face(g, curr_id, ID)
    ret = -1;
    if intersect(g.faces{curr_id}.nodes, g.faces{ID}.nodes) == 1
        g.faces{ID}.ajac_faces = union(g.faces{ID}.ajac_faces, curr_id);
    end
end