function f=centre_element(ie)

global vcor
global kconec

f=mean(vcor(kconec(ie,:),1:2))';