function f=display_multilayer(m)

for ii=1:size(m,2)
    disp(['#####################']);
    disp(['Layer #' num2str(ii)]);
   disp(['Material =' num2str(m(ii).mat)]);
   disp(['thickness=' num2str(m(ii).d)]);
end
disp(['#####################']);
f=0;