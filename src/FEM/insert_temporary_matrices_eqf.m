[i_temp,j_temp,s_temp] = find(vh);
ll=length(i_temp);
i_h_eqf=[i_h_eqf;index_p(i_temp)'];
j_h_eqf=[j_h_eqf;index_p(j_temp)'];
v_h_eqf(1:end+ll,element_num_mat(ie))=[v_h_eqf(:,element_num_mat(ie));s_temp];
[i_temp,j_temp,s_temp] = find(vq);
ll=length(i_temp);
i_q_eqf=[i_q_eqf;index_p(i_temp)'];
j_q_eqf=[j_q_eqf;index_p(j_temp)'];
v_q_eqf(1:end+ll,element_num_mat(ie))=[v_q_eqf(:,element_num_mat(ie));s_temp];

