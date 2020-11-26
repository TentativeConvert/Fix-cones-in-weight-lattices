restart
loadPackage "WeylGroupsFixed"

angle = method(TypicalValue => RR)
angles = method(TypicalValue => Table)

angle(RootSystem,ZZ,ZZ) := (R,a,b) -> (
    sc := scalarProduct(R,a,b)/(sqrt(scalarProduct(R,a,a))*sqrt(scalarProduct(R,b,b)));
    return acos(sc)/pi*180
    )    	 
angles(RootSystem) := (R) -> (
    return table(toList(1..rank R),toList(1..rank R), (i,j) -> format(3,3,angle(R,i,j)))
)

R = rootSystemA(2); netList angles(R)
R = rootSystemA(9); netList angles(R)
R = rootSystemB(9); netList angles(R)
R = rootSystemC(9); netList angles(R)
R = rootSystemD(9); netList angles(R)
R = rootSystemF4; netList angles(R)
R = rootSystemE(6); netList angles(R)
R = rootSystemE(7); netList angles(R)
R = rootSystemE(8); netList angles(R)
