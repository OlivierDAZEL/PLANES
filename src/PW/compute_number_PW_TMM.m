for ii=1:nb_layers
    switch floor(multilayer(ii).mat/1000)
        case 1
            n_w(ii)=4;
        case {0 2 3}
            n_w(ii)=2;
        case {4 5}
            n_w(ii)=6;
    end
end
nb_amplitudes=sum(n_w);


if termination~=0
    nb_amplitudes=nb_amplitudes+1;
end
