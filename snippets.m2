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
R = rootSystemE(6)
omega1 = fundamentalWeight(R,1)
omega6 = fundamentalWeight(R,6)
omega3 = fundamentalWeight(R,3)
omega5 = fundamentalWeight(R,5)
w0 = longestWeylGroupElement(R)
W = {};
for i from 0 to coxeterLength(longestWeylGroupElement(R)) do (
    W = W|listWeylGroupElements(R,i)
    )

orbit1 = {}
orbit6 = {}
orbit3 = {}
orbit5 = {}
orbit16 = {}
orbit35 = {}
for w in W do (
    orbit1 = orbit1|{w*omega1};
    orbit3 = orbit3|{w*omega3};
    orbit5 = orbit5|{w*omega5};
    orbit6 = orbit6|{w*omega6};
    orbit16 = orbit16|{w*(omega1 + omega6)};
    orbit35 = orbit35|{w*(omega3 + omega5)};
    )
orbit1 = unique(orbit1)
orbit3 = unique(orbit3)
orbit5 = unique(orbit5)
orbit6 = unique(orbit6)
orbit16 = unique(orbit16)
orbit35 = unique(orbit35)

dominants16 = {}
for omega in orbit16 do (    
    if isInClosedChamber(omega) then (	  
	dominants16 = dominants16|{omega}
	)
    )
dominants16

dominants = method(TypicalValue => List)
dominants(List) := (L) -> (
    newL = {};	  
    for omega in L do (    
    	if isInClosedChamber(omega) then (	  
	    newL = newL|{omega}
	    )    
    	);	  
    return newL;
    )

dominants35 = dominants(orbit35)
dominants16 = dominants(orbit16)


   
    	
