function Q=multiply_polynom_2D(P1,P2)

Q=zeros(size(P1,1)+size(P2,1)-1,size(P1,2)+size(P2,2)-1);
for x1=1:size(P1,1)
    for y1=1:size(P1,2)
        for x2=1:size(P2,1)
            for y2=1:size(P2,2)
                Q(x1+x2-1,y1+y2-1)=Q(x1+x2-1,y1+y2-1)+P1(x1,y1)*P2(x2,y2);
            end
        end
        
    end
end

