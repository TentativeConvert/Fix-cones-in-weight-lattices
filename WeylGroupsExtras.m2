newPackage(
	"WeylGroupsExtras",
	Headline => "Some auxiliarly methods for the package WeylGroups",
	Version => "0.0.1",
	Date => "2020",
        Authors => {{Name => "Marcus Zibrowius"}},
	PackageExports => {"WeylGroups"},	   
	DebuggingMode => true)
-- Put here the name of functions that should be visible to users
export{
    "fundamentalWeight","isInClosedChamber","levyComponent","otherComponent","weylOrbit",
    "prettyfy", "texify"
}


fundamentalWeight = method(TypicalValue => Weight)
fundamentalWeight(RootSystem,ZZ) := (R,i) -> (
    if i < 1 then error "fundamental weights cannot have index < 1";
    if i > rank R then error "fundamental weights cannot have index > rank of root system";
    vec := toList(apply(1..rank(R), j -> if i==j then 1 else 0));
    return weight(R,vec)
    )    


isInClosedChamber = method(TypicalValue => Boolean)
isInClosedChamber(Weight) := (w) -> (
    for i from 0 to length(entries(w))-1 do (
	if w_i < 0 then return false
	);
    return true
    )
isInClosedChamber(Parabolic,Weight) := (P,w) -> (
    -- no error handling if weight and parabolic belong to different root systems! --
    for i in toList(P) do (
	if w_(i-1) < 0 then return false
	);
    return true
    )	    	     


levyComponent = method(TypicalValue => Vector)
levyComponent(RootSystem,Parabolic,Weight) := (R,P,omega) -> (
    (ZZ^(rank R))^(apply(sort toList P, i -> i-1))*(vector entries omega)    
    )
otherComponent = method(TypicalValue => Vector)
otherComponent(RootSystem,Parabolic,Weight) := (R,P,omega) -> (
    (ZZ^(rank R))^(apply(sort toList (set(1..rank R) - P), i -> i-1))*(vector entries omega)
    )    


weylOrbit = method(TypicalValue => List)
weylOrbit(RootSystem,Weight) := (R,w) -> (
    if not isInClosedChamber(w) then error "computation of Weyl orbit currently only implemented starting from /heighest/ weight in the orbit";
    orbitlevels := {{w}};
    -- commented lines with "reflections" only here for debugging purposes
    -- reflections := {};
    while last orbitlevels != {} do (
	nextlevel := {};	
	-- nextreflections := {};
    	for w in last orbitlevels do (	 
    	    for i from 1 to length(entries(w)) do (
	    	if w_(i-1) > 0 then (
		    -- nextreflections = append(nextreflections,i);
		    nextlevel = append(nextlevel,reduce(R,{i})*w)
		    )
	    	)    
    	    );   
	-- reflections = append(reflections,nextreflections);
	orbitlevels = append(orbitlevels,unique(nextlevel))	   
    	);
    -- print reflections;
    return flatten(orbitlevels);
)    	 
-- with R = F4:
-- benchmark "o = weylOrbit(R,fundamentalWeight(R,4))" --> 0.0064
-- benchmark "o = weylOrbit(R,fundamentalWeight(R,3))" --> 0.033
-- benchmark "o = weylOrbit(R,fundamentalWeight(R,2))" --> 0.033
-- benchmark "o = weylOrbit(R,fundamentalWeight(R,1))" --> 0.0064
----------------------------------------------------------------
-- total: 0.14


--------------------------------------------------
-- printing methods:
prettyfy = method(TypicalValue => String)
texify = method(TypicalValue => String)

prettyfy(Parabolic) := (P) -> (
    return toString sort toList(P)
    )
texify(Parabolic) := (P) -> (
    s := toString sort toList(P);
    s = replace("[{]","\\{",replace("[}]","\\}",s));
    return "\\("|s|"\\)"
    )
prettyfy(DynkinType) := (type) -> (
    s := "";
    for component in toList(type) do s = s|toString(component_0)|toString(component_1)|" x ";
    if length(s) > 0 then s=substring(s,0,length(s)-3); -- cut off last "\\times"
    return s
    )	 
texify(DynkinType) := (type) -> (
    s := "";
    for component in toList(type) do s = s|toString(component_0)|"_"|toString(component_1)|" \\times ";
    if length(s) > 0 then s=substring(s,0,length(s)-7); -- cut off last "\\times"
    return "\\("|s|"\\)"
    )

----------------------------------------------------------------------------------------------------

beginDocumentation() -- very sensitive to indentation
doc ///
Key
  WeylGroupsExtras
Headline
  WeylGroupsExtras
Description
 Text
  An auxiliary package that slightly extends the package "WeylGroups" by Baptiste Calmès and  Viktor Petrov.
  Developed for research related to the manuscript:
  Tobias Hemmert &  Marcus Zibrowius, The Witt rings of many flag varieties are exterior algebras.
///
TEST ///
R = rootSystemE(6)
print fundamentalWeight(R,6)
///
TEST ///
R = rootSystemA(6)
P = parabolic(R,set(4,1))
lc = omega -> levyComponent(R,P,omega)
oc = omega -> otherComponent(R,P,omega)

omega1 = new Weight from vector{-1,-2,-3,-4,-5,-6}
omega2 = new Weight from vector{10,20,30,40,50,60}

lc(omega1)
lc(omega2)
oc(omega1)
oc(omega2)
///
TEST ///
R = rootSystemE(7)
P = parabolic(R,set{7,1,2,3})

texify(P)
prettyfy(P)
t = dynkinType(rootSystem(R,P))
texify(t)
prettyfy(t)
///
