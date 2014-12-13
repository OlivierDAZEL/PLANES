if(nb_frequencies>0)
    if nb_frequencies==1
        vec_freq=freq_min;
    else
        vec_freq=linspace(freq_min,freq_max,nb_frequencies);
    end
else
    vec_freq=logspace(log10(freq_min),log10(freq_max),abs(nb_frequencies));
end