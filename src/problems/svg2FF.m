

itemp=1;
x1(1)=temp(itemp);
itemp=itemp+1;
x1(2)=temp(itemp);

itemp=itemp+1;
v1(1)=x1(1)+temp(itemp);
itemp=itemp+1;
v1(2)=x1(2)+temp(itemp);

itemp=itemp+1;
v2(1)=x1(1)+temp(itemp);
itemp=itemp+1;
v2(2)=x1(2)+temp(itemp);

itemp=itemp+1;
x2(1)=x1(1)+temp(itemp);
itemp=itemp+1;
x2(2)=x1(2)+temp(itemp);

x=x1(1).*(1-t).^3+3*v1(1).*t.*(1-t).^2+3*v2(1).*t.^2.*(1-t)+x2(1).*t.^3;
y=x1(2).*(1-t).^3+3*v1(2).*t.*(1-t).^2+3*v2(2).*t.^2.*(1-t)+x2(2).*t.^3;

delta_x=(x(2:end)-x(1:end-1));
delta_y=(y(2:end)-y(1:end-1));
l2=sqrt(delta_x.^2+delta_y.^2);


i_border=i_border+1;
n_border(i_border)=sign_boundary*ceil(sum(l2));


text_file=['border b' num2str(i_border) '(t=0,1){x=' num2str(x1(1)/1000) '*(1.0-t)^3+3.0*' num2str(v1(1)/1000) '*t*(1.0-t)^2+3*' num2str(v2(1)/1000) '*t^2*(1-t)+' num2str(x2(1)/1000) '*t^3 ;y=' num2str(x1(2)/1000) '*(1.0-t)^3+3.0*' num2str(v1(2)/1000) '*t*(1.0-t)^2+3*' num2str(v2(2)/1000) '*t^2*(1-t)+' num2str(x2(2)/1000) '*t^3;label=3;}'];
fprintf(fid,'%s\n',text_file);
plot(x,y)

for ii=2:floor(length(temp)/6)
    x1=x2;
    itemp=itemp+1;
    v1(1)=x1(1)+temp(itemp);
    itemp=itemp+1;
    v1(2)=x1(2)+temp(itemp);
    
    itemp=itemp+1;
    v2(1)=x1(1)+temp(itemp);
    itemp=itemp+1;
    v2(2)=x1(2)+temp(itemp);
    
    itemp=itemp+1;
    x2(1)=x1(1)+temp(itemp);
    itemp=itemp+1;
    x2(2)=x1(2)+temp(itemp);
    
    
    x=x1(1).*(1-t).^3+3*v1(1).*t.*(1-t).^2+3*v2(1).*t.^2.*(1-t)+x2(1).*t.^3;
    y=x1(2).*(1-t).^3+3*v1(2).*t.*(1-t).^2+3*v2(2).*t.^2.*(1-t)+x2(2).*t.^3;
    
    delta_x=(x(2:end)-x(1:end-1));
    delta_y=(y(2:end)-y(1:end-1));
    l2=sqrt(delta_x.^2+delta_y.^2);
    
    
    
    i_border=i_border+1;
    n_border(i_border)=sign_boundary*ceil(sum(l2));
    text_file=['border b' num2str(i_border) '(t=0,1){x=' num2str(x1(1)/1000) '*(1.0-t)^3+3.0*' num2str(v1(1)/1000) '*t*(1.0-t)^2+3*' num2str(v2(1)/1000) '*t^2*(1-t)+' num2str(x2(1)/1000) '*t^3 ;y=' num2str(x1(2)/1000) '*(1.0-t)^3+3.0*' num2str(v1(2)/1000) '*t*(1.0-t)^2+3*' num2str(v2(2)/1000) '*t^2*(1-t)+' num2str(x2(2)/1000) '*t^3;label=3;}'];
    fprintf(fid,'%s\n',text_file);
    
    plot(x,y)
end