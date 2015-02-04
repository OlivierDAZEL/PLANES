function vec_edge=normal_edge(coord)

vec_edge=coord(:,2)-coord(:,1);
vec_edge=vec_edge/norm(vec_edge);
vec_edge=[0 1;-1 0]*vec_edge;
