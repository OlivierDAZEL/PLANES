[i_temp,j_temp,s_temp] = find(vh);
ll=length(i_temp);
i_h_limp=[i_h_limp;index_p(i_temp)'];
j_h_limp=[j_h_limp;index_p(j_temp)'];
v_h_limp(1:end+ll,element_num_mat(ie))=[v_h_limp(:,element_num_mat(ie));s_temp];
[i_temp,j_temp,s_temp] = find(vq);
ll=length(i_temp);
i_q_limp=[i_q_limp;index_p(i_temp)'];
j_q_limp=[j_q_limp;index_p(j_temp)'];
v_q_limp(1:end+ll,element_num_mat(ie))=[v_q_limp(:,element_num_mat(ie));s_temp];

