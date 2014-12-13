function f=uxy(node)

f(1:2:2*length(node))=(3*(node-1)+1);
f(2:2:2*length(node))=(3*(node-1)+2);