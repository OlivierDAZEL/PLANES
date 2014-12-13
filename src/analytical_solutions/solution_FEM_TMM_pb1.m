% Calcul de la solution analytique
        k_parallel=k_air*sin(theta);
        k_perp=sqrt(k_eq^2-k_parallel^2);
        k_z=k_air*cos(theta);
        %[R A B C D]
        x=0;
        M(1,4)=-1;
        M(1,5)=1;
        x=-dy;
        M(2,2)=-k_perp*exp(-j*k_perp*x);
        M(2,3)= k_perp*exp( j*k_perp*x);
        M(2,4)= k_z*exp(-j*k_z*x);
        M(2,5)=-k_z*exp( j*k_z*x);
        M(3,2)=K_eq_til*k_eq^2*exp(-j*k_perp*x);
        M(3,3)=K_eq_til*k_eq^2*exp( j*k_perp*x);
        M(3,4)=-K_0*k_air^2*exp(-j*k_z*x);
        M(3,5)=-K_0*k_air^2*exp( j*k_z*x);
        x=-dy-delta_y_MMT;
        FF(4,1)=-k_z*exp(-j*k_z*x);
        M(4,1)= k_z*exp( j*k_z*x);
        M(4,2)= k_perp*exp(-j*k_perp*x);
        M(4,3)=-k_perp*exp( j*k_perp*x);
        FF(5,1)=K_0*k_air^2*exp(-j*k_z*x);
        M(5,1)= K_0*k_air^2*exp( j*k_z*x);
        M(5,2)=-K_eq_til*k_eq^2*exp(-j*k_perp*x);
        M(5,3)=-K_eq_til*k_eq^2*exp( j*k_perp*x);
        
        X=M\FF;
        R_analytic(i_theta)=-X(1);
        abs_analytic(i_theta)=1-abs(R_analytic(i_theta))^2;