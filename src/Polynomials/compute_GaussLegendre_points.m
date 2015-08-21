function Gauss_points= compute_GaussLegendre_points(N)

Gauss_points.nb=N;
L_n=Legendre_polynomial(N);
Gauss_points.xi=roots(L_n(end:-1:1))';
L_n_p=derive_polynom(L_n);
Gauss_points.w=zeros(1,N);

for ii=1:N
    Gauss_points.w(ii)=0;
    for jj=1:length(L_n_p)
        Gauss_points.w(ii)=Gauss_points.w(ii)+L_n_p(jj)*Gauss_points.xi(ii)^(jj-1);
    end
    Gauss_points.w(ii)=2/((1-Gauss_points.xi(ii)^2)*Gauss_points.w(ii)^2);
end


