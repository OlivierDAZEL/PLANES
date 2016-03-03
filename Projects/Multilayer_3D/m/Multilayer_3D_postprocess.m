fid_out=fopen(name.file_PW3D,'a');
switch termination(i_m)
    case 0
        fprintf(fid_out,'%12.8f\t%12.8f\t%12.8f\t%12.8f\n',frequency.vec(i_f),abs_PW_3D(i_f),real(rflx_PW_3D(i_f)),imag(rflx_PW_3D));
    case 1
        fprintf(fid_out,'%12.8f\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n',frequency.vec(i_f),abs_PW_3D(i_f),real(rflx_PW_3D(i_f)),imag(rflx_PW_3D),TL_PW_3D(i_f));
    otherwise
        INCORRECT_TERMINATION_FOR_MULTILAYER
end
fclose(fid_out);
