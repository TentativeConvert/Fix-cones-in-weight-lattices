restart
loadPackage "WeylGroups"
viewHelp "WeylGroups"
R = rootSystemE(6)
s1    = reduce(R,{1})
s1min = minimalRepresentative(s1 % P)
s1minn = minimalRepresentative(P % s1)
s1*FundamentalWeight
s1min*FundamentalWeight
s1minn*FundamentalWeight

a1 = simpleRoot(R,2)
a1_0
a1_1
a1_2
b0 = weight(R,{0,0,0,0,0,1})
b1 = reduce(R,{6})*b0
b2 = reduce(R,{5})*b1
b3 = reduce(R,{4})*b2
b41 = reduce(R,{3})*b3
b42 = reduce(R,{2})*b3
rank R

restart
loadPackage "WeylGroupsExtras"
check "WeylGroupsExtras"


    
