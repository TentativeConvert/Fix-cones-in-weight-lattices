-- Tobias' program:
restart
loadPackage "WeylGroups"
R=rootSystemF4
fundweights={weight(R,{1,0,0,0}),weight(R,{0,1,0,0}),weight(R,{0,0,1,0}),weight(R,{0,0,0,1})}
w1=longestWeylGroupElement(R)
a=coxeterLength(w1)
L={}
for i from 0 to a do L=L|listWeylGroupElements(R,i)
-- benchmark with R = F4:
-- benchmark "for i from 0 to a do L=L|listWeylGroupElements(R,i)"  -- 0.76
#L

S={1,4}
P=parabolic(R,set(S))
w0 = longestWeylGroupElement(R,P)
-- [...]
    
X={}
for fw in fundweights do (X=append(X,unique(delete(null,apply(L,
		    i->if isMinimalRepresentative(P,i) and -w0*i*fw==i*fw then (i*fw))))));
X

