restart
loadPackage "WeylGroupsExtras"
--check "WeylGroupsExtras"
loadPackage "Auxiliary"
--check "Auxiliary"

--------------------------------------------------
-- administration of input & output
RootSystemPair = new Type of MutableHashTable  -- needed for sorting only
rootSystemPair = method(TypicalValue => RootSystemPair)
rootSystemPair(RootSystem,Parabolic) := (R,P) -> (new RootSystemPair from {
	"rootsystem" => R,
	"parabolic" => P
	}
    )
RootSystemPair ? RootSystemPair := (RSP1,RSP2) -> (
    L := (RSP) -> (
	R := RSP#"rootsystem";
	P := RSP#"parabolic";
	return {#dynkinType(rootSystem(R,P))}|{sort toList(dynkinType(rootSystem(R,P)))}|{sort toList P}
	);
    return L(RSP1) ? L(RSP2)
    )

ResultLine = new Type of MutableHashTable
resultLine = method(TypicalValue => ResultLine)
resultLine(RootSystem,Parabolic) := (R,P) -> (new ResultLine from {
	    "type" => dynkinType(rootSystem(R,P)),
    	    "parabolic" => P
	    })

writeOut = method()
writeOut(RootSystem,List) := (R,results) -> (
    file := ("results_"|prettyfy(dynkinType(R))|".tex") << endl;
    print "--------------------------------------------------------------------------------------------------";	   
    print (pad("parabolic",30)|pad("Dynkin type", 50)|pad("(single cell)",13)|pad("(free)",7)|pad("(orbit basis)",20));
    print "--------------------------------------------------------------------------------------------------";	   
    for result in results do (
    	assert(result#"free1" == result#"free2");
	strFree    = (if result#"free1" then "\\Free" else "\\no");
	strSingleCell = (if result#"singlecell" then "\\SingleCell" else "\\no");
        
	if member("orbitbasis",keys result) then (
	    strOrbitBasis = (if result#"orbitbasis" then "\\OrbitBasis" else "\\no");
	    ) else strOrbitBasis = "???";
	
	if result#"singlecell" then (
	    assert(result#"free1");
	    strFree = "("|strFree|")";
	    if member("orbitbasis",keys result) then assert(result#"orbitbasis");
	    strOrbitBasis = "("|strOrbitBasis|")";
	    ) else (	 
	    strFree = " "|strFree|" ";
	    if result#"free1" then strOrbitBasis = " "|strOrbitBasis|" " else strOrbitBasis = "("|strOrbitBasis|")";    
	    );
    	
	print (pad(prettyfy(result#"parabolic"),30)|pad(prettyfy(result#"type"),50)|pad(strSingleCell,13)|pad(strFree,7)|pad(strOrbitBasis,20));
	print "--------------------------------------------------------------------------------------------------";	   
	file << pad(texify(result#"parabolic"),30)|" & "|pad(texify(result#"type"),50)|" & "|pad(strSingleCell,13)|" & "|pad(strFree,7)|" & "|pad(strOrbitBasis,20)|" \\\\" << endl
	);	 
    file << close
)

divideWeightsByTwoIfPossible = method(TypicalValue => List) 
divideWeightsByTwoIfPossible(RootSystem,List) := (R,L) -> (
    -- auxiliary method that divides each weight in a list of weights by two, if possible
    for omega in L do (
	if substitute(vector entries omega,ZZ/2) == 0_((ZZ/2)^(rank R)) then (
	    L = delete(omega, L);
	    L = append(L, weight(R,(entries omega)//2))
	    )
	);
    return L
    )
 
distinguishedGeneratorsOfFixedPointMonoid = method(TypicalValue => List)
distinguishedGeneratorsOfFixedPointMonoid(RootSystem,Parabolic) := (R,P) -> (
    -- distinguished generators of FixedPointMonoid
    -- linearly independent, and generate the real fixed point cone 
    wP := longestWeylGroupElement(R,P);    
    distinguishedGenerators := apply(toList(P), i -> fundamentalWeight(R,i));
    distinguishedGenerators = unique(apply(distinguishedGenerators, omega -> omega-wP*omega));
    return (divideWeightsByTwoIfPossible(R,distinguishedGenerators))
    )

furtherGeneratorsOfFixedPointMonoid = method(TypicalValue => List)
furtherGeneratorsOfFixedPointMonoid(RootSystem,Parabolic) := (R,P) -> (    
    -- potential further integral generators of FixedPointMonoid
    wP := longestWeylGroupElement(R,P);    
    PSgeq2 := delete(null,apply(powerSet(P), S -> if #S >= 2 then S else null));                 -- subsets of P with at least 2 elements
    furtherGenerators := apply(PSgeq2, S -> sum(apply(toList(S), i -> fundamentalWeight(R,i)))); -- linear combinations of at least two fundamentalWeights w_i with i in P
    furtherGenerators = unique(apply(furtherGenerators, omega -> omega-wP*omega));
    return  divideWeightsByTwoIfPossible(R,furtherGenerators)
    )

symmetrizedFundamentalWeights = method(TypicalValue => List)
symmetrizedFundamentalWeights(RootSystem) := (R) -> (
    wR := longestWeylGroupElement(R);    
    fundamentalWeights := toList(apply(1..rank R, j -> fundamentalWeight(R,j)));
    return (unique(apply(fundamentalWeights, w -> if -wR*w == w then w else w -wR*w)))
)

checkIfFixedPointMonoidIsFree1 = method(TypicalValue => Boolean)
checkIfFixedPointMonoidIsFree1(RootSystem,Parabolic) := (R,P) -> (
    -- check for freenees via ZZ/2 reduction
    wP := longestWeylGroupElement(R,P);
    lc := omega -> levyComponent(R,P,omega);
    oc := omega -> otherComponent(R,P,omega);
    M := map(ZZ^(rank R - #P),ZZ^0,0);
    for i from 1 to rank R do (
    	omega = fundamentalWeight(R,i);
    	if lc(omega) == lc(-wP*omega) then M = M|matrix{oc(omega-wP*omega)}
    	);
    return isStandardSubspace(kernel substitute(M,ZZ/2))
    )

checkIfFixedPointMonoidIsFree2 = method(TypicalValue => Boolean)
checkIfFixedPointMonoidIsFree2(RootSystem,Parabolic) := (R,P) -> (
    -- check for freeness more directly
    result := true;
    distGens    := distinguishedGeneratorsOfFixedPointMonoid(R,P);
    furtherGens := furtherGeneratorsOfFixedPointMonoid(R,P);
    if #distGens>0 and #furtherGens>0 then (
    	M1 := transpose matrix(apply(distGens, omega -> entries omega));
    	M2 := transpose matrix(apply(furtherGens, omega -> entries omega));
	result = (image(M1) == image(M1|M2))
	);
    return result
)

checkIfFixedPointMonoidSatisfiesSingleCell = method(TypicalValue => Boolean)
checkIfFixedPointMonoidSatisfiesSingleCell(RootSystem,Parabolic) := (R,P) -> (
    result := true;
    distGens  := distinguishedGeneratorsOfFixedPointMonoid(R,P);	
    for alpha in toList(positiveRoots(R)) do (
    	signs := apply(distGens, omega -> (eval(R,omega,alpha) ? 0));
    	if member(symbol <, signs) and member(symbol >, signs) then (
	    result = false;
	    break
    	    )
	);
    return result
    )

checkIfFixedPointMonoidSatisfiesOrbitBasis = method(TypicalValue => Boolean)
checkIfFixedPointMonoidSatisfiesOrbitBasis(RootSystem,Parabolic) := (R,P) -> (
    wP := longestWeylGroupElement(R,P);
    -- find selfdual weights in FixedPointMonoid
    selfDualWeightsInOrbits := {};
    for omega in symmetrizedFundamentalWeights(R) do (
	omegasOrbit := weylOrbit(R,omega);
	omegasSelfDualWeights := {};
	for e in omegasOrbit do (
	    if isInClosedChamber(P,e) and -wP*e==e then omegasSelfDualWeights = append(omegasSelfDualWeights,e)
	    );
	selfDualWeightsInOrbits = append(selfDualWeightsInOrbits,omegasSelfDualWeights)
	);
    -- check condition for getting an exterior algebra:
    result := true;	   
    for omega in distinguishedGeneratorsOfFixedPointMonoid(R,P) do (
    	if not member(toList{omega},selfDualWeightsInOrbits) then (
	    result = false;
	    break;
	    )    
    	);	
    return result
)

--------------------------------------------------
-- computational core:  analysing the situation for a specific parabolic
analyseParabolic = method(TypicalValue => ResultLine)    
analyseParabolic(RootSystem,Parabolic) := (R,P) -> (
    result := resultLine(R,P);
    result#"singlecell" = checkIfFixedPointMonoidSatisfiesSingleCell(R,P);
    result#"free1" = checkIfFixedPointMonoidIsFree1(R,P);
    result#"free2" = checkIfFixedPointMonoidIsFree2(R,P);
    if rank(R) < 8 then result#"orbitbasis" = checkIfFixedPointMonoidSatisfiesOrbitBasis(R,P);
    return result
    )

--------------------------------------------------
-- the MAIN routine

main = method()
main(List) := (listOfRootSystems) -> (
    for R in listOfRootSystems do (
    	listOfRootSystemPairs = sort(apply(powerSet(set(1..rank R)), S -> rootSystemPair(R,parabolic(R,S)))); 
    	results = {};
    	for RSP in listOfRootSystemPairs do ( 
       	    print (prettyfy(dynkinType(RSP#"rootsystem"))|"  "|prettyfy(RSP#"parabolic")|"  ...   ");
    	    results = append(results,analyseParabolic(RSP#"rootsystem",RSP#"parabolic"))
    	    
	    );
    	print "Printing ...";
    	writeOut(R,results)
    	)
    )    

-- main({rootSystemA(7),rootSystemB(7),rootSystemC(7),rootSystemD(7)})
main({rootSystemG2,rootSystemF4,rootSystemE(6),rootSystemE(7),rootSystemE(8)})

