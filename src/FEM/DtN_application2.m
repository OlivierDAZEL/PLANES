if (nb_R+nb_T)>0
        A_2(nb_dof_FEM+nb_R+nb_T,nb_dof_FEM+nb_R+nb_T)=0;
        F_2(nb_dof_FEM+nb_R+nb_T,1)=0;
end
    
    DtN_porous_R=0;
    DtN_porous_T=0;
    DtN_elas_R=0;
    DtN_elas_T=0;
    
    for ie=1:nb_loads  
        typ=loads(ie,4);
        switch typ
            case {10}
                DtN_porous_R=1;
                excitation_10000
            case {11}
                DtN_elas_R=1;
                excitation_20000
            case {21}
                DtN_elas_T=1;
                excitation_21000
            case {12}
                DtN_porous_R=1;
                excitation_50000
            otherwise
                disp('Unknown load')
                stop
        end
    end
    
 
    
    if (nb_media.eqf*DtN_porous_R)~=0
        %    DtN Biot 1998
        if (nb_R~=0)
            F_2(nb_dof_FEM+nb_Bloch_waves+1)=period;
            for i_R=1:nb_R
                 A_2(nb_dof_FEM+i_R,nb_dof_FEM+i_R)=A_2(nb_dof_FEM+i_R,nb_dof_FEM+i_R)-period;
            end
        end
    end
   
        if (nb_media.limp*DtN_porous_R)~=0
        %    DtN Biot 1998
        if (nb_R~=0)
            F_2(nb_dof_FEM+nb_Bloch_waves+1)=period;
            for i_R=1:nb_R
                 A_2(nb_dof_FEM+i_R,nb_dof_FEM+i_R)=A_2(nb_dof_FEM+i_R,nb_dof_FEM+i_R)-period;
            end
        end
    end
    
    
    
    
if (nb_media.pem98*DtN_porous_R)~=0
        %    DtN Biot 1998
        if (nb_R~=0)
            F_2(nb_dof_FEM+nb_Bloch_waves+1)=period;
            for i_R=1:nb_R
                 A_2(nb_dof_FEM+i_R,nb_dof_FEM+i_R)=A_2(nb_dof_FEM+i_R,nb_dof_FEM+i_R)-period;
            end
        end
    end  
    
    
    
    if (nb_media.pem01*DtN_porous_R)~=0
        %    DtN Biot 2001
        if (nb_R~=0)
            F_2(nb_dof_FEM+1)=period;
            for i_R=1:nb_R
                A_2(nb_dof_FEM+i_R,nb_dof_FEM+i_R)=A_2(nb_dof_FEM+i_R,nb_dof_FEM+i_R)-period;
            end
        end
    end
    
    
    
    
    
    if (nb_media.elas*DtN_elas_R)~=0
        %    DtN Biot elas
        if (nb_R~=0)
            F_2(nb_dof_FEM+nb_Bloch_waves+1)=F_2(nb_dof_FEM+nb_Bloch_waves+1)+period*(1i*k_z)/(rho_0*omega^2);
            for i_R=1:nb_R
                A_2(nb_dof_FEM+i_R,nb_dof_FEM+i_R)=A_2(nb_dof_FEM+i_R,nb_dof_FEM+i_R)-...
                    period*(1i*vec_k_z(i_R))/(rho_0*omega^2);
            end
        end
        if (nb_T~=0)
            for i_T=1:nb_T
                A_2(nb_dof_FEM+nb_R+i_T,nb_dof_FEM+nb_R+i_T)=A_2(nb_dof_FEM+nb_R+i_T,nb_dof_FEM+nb_R+i_T)+...
                    period*(1i*vec_k_z(i_T))/(rho_0*omega^2);
            end
        end
    end