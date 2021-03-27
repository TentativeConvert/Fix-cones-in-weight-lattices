restart
loadPackage "WeylGroupsExtras"
R = rootSystemE(6)
omega1 = fundamentalWeight(R,1)
omega6 = fundamentalWeight(R,6)
omega3 = fundamentalWeight(R,3)
omega5 = fundamentalWeight(R,5)
w0 = longestWeylGroupElement(R)

assert(omega1 == -w0*omega6) -- omega1 & omega6 are dual
assert(omega3 == -w0*omega5) -- omega3 & omega5 are dual
    
computeWeylGroup = method(TypicalValue => List);
computeWeylGroup(RootSystem) := (R) -> (
    print("Computing Weyl group ...");
    w0 := longestWeylGroupElement(R);
    l  := coxeterLength(w0);
    W := {};
    for i from 0 to l do (
    	print("... "|toString(round(i/l*100.0))|" %");
    	W = W|listWeylGroupElements(R,i);
    	);
    print("... Done.");    
    return W;
    );

W = computeWeylGroup(R);    

orbit = method(TypicalValue => List)
orbit(Weight) := (omega) -> (
    myorbit := {};
    for w in W do (
	myorbit = myorbit|{w*omega}
	);
    return unique(myorbit);
    );

print("Computing orbit of omega_1 ...");
orbit1 = orbit(omega1);
print("Computing orbit of omega_3 ...");
orbit3 = orbit(omega3);
print("Computing orbit of omega_5 ...");
orbit5 = orbit(omega5);
print("Computing orbit of omega_6 ...");
orbit6 = orbit(omega6);
print("... Done.");

#orbit1 -- 27
#orbit3 -- 216
#orbit5 -- 216
#orbit6 -- 27

dominants = method(TypicalValue => List)
dominants(List) := (L) -> (
    newL := {};	  
    for omega in L do (    
    	if isInClosedChamber(omega) then (	  
	    newL = newL|{omega}
	    )    
    	);	  
    return newL;
    )

dominants1 = dominants(orbit1) -- should be just omega1
dominants3 = dominants(orbit3) -- should be just omega3
dominants5 = dominants(orbit5) -- should be just omega5
dominants6 = dominants(orbit6) -- should be just omega6


productorbit = method(TypicalValue => List)
productorbit(List,List) := (L1,L2) -> (
    newL := {};
    for omega1 in L1 do (
	for omega2 in L2 do (
	    newL = newL|{omega1+omega2}
	    )
	);
    return newL;
    )	  

productorbit16 = productorbit(orbit1,orbit6)    
productorbit35 = productorbit(orbit3,orbit5)    

dominants16 = dominants(productorbit16)
dominants35 = dominants(productorbit35)

tally(dominants16)
-- S(ω1)S(ω6) = S(ω1+ω6) + 6 S(ω2) + 27

tally(dominants35)
-- S(ω3)S(ω5) = S(ω3+ω5) + 4 S(ω1+ω2+ω6) + 10 S(ω5+ω6) +  10 S(ω1+ω3) + 15 S(2ω2) + 18 S(ω4) + 32 S(ω1+ω6) + 60 S(ω2) + 216
   
    	

