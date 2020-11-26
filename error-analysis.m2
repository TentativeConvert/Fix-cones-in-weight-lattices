restart
loadPackage "WeylGroupsExtras"
loadPackage "Auxiliary"
R = rootSystemC(5)
P = parabolic(R,set{2,3,4})
rcritical = addRoots(R,simpleRoot(R,4),addRoots(R,simpleRoot(R,4),simpleRoot(R,5)))
rootCoefficients(R,rcritical)

-- from main.m2:
print (prettyfy(dynkinType(R))|"  "|prettyfy(P)|"  ...   ");
result := resultLine(R,P);
wP := longestWeylGroupElement(R,P);
        
--------------------------------------------------
-- check for freeness via ZZ/2-reduction
lc = omega -> levyComponent(R,P,omega);
oc = omega -> otherComponent(R,P,omega);
M = map(ZZ^(rank R - #P),ZZ^0,0);
for i from 1 to rank R do (
    omega = fundamentalWeight(R,i);
    if lc(omega) == lc(-wP*omega) then M = M|matrix{oc(omega-wP*omega)}
    );
print isStandardSubspace(kernel substitute(M,ZZ/2))

--------------------------------------------------
-- check for freeness more directly
divideByTwoIfPossible = method(TypicalValue => List);  -- auxiliary method that divides each weight in a list of weights by two, if possible
divideByTwoIfPossible(RootSystem,List) := (R,L) -> (
    for omega in L do (
	if substitute(vector entries omega,ZZ/2) == 0_((ZZ/2)^(rank R)) then (
	    L = delete(omega, L);
	    L = append(L, weight(R,(entries omega)//2))
	    )
	);
    return L;
    );
-- distinguished generators (potential basis) of L_\Theta:
distinguishedGenerators = apply(toList(P), i -> fundamentalWeight(R,i));
distinguishedGenerators = unique(apply(distinguishedGenerators, omega -> omega-wP*omega));
distinguishedGenerators = divideByTwoIfPossible(R,distinguishedGenerators);

matrix toList(apply(1..5, i -> toList(apply(1..5, j -> eval(R,fundamentalWeight(R,i),simpleRoot(R,j))))))
matrix toList(apply(1..5, i -> toList(apply(1..5, j -> scalarProduct(R,fundamentalWeight(R,i),simpleRoot(R,j))))))

rootCoefficients(R,rcritical)
e1 = distinguishedGenerators_0
e2 = distinguishedGenerators_1

scalarProduct(R,e1,simpleRoot(R,4))
scalarProduct(R,e1,simpleRoot(R,5))
scalarProduct(R,e1,rcritical)

scalarProduct(R,e2,simpleRoot(R,4))
scalarProduct(R,e2,simpleRoot(R,5))
scalarProduct(R,e2,rcritical)

eval(R,e1,rcritical)
eval(R,e2,rcritical)
