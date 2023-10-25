[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5742827.svg)](https://doi.org/10.5281/zenodo.5742827)

The code in this repository complements the following research manuscript:

[HZ] Tobias Hemmert & Marcus Zibrowius, *The Witt rings of many flag varieties are exterior algebras*

## Prerequisites

The code is written in Macaulay2 and requires the `WeylGroups` package, written by Baptiste Calmès and Viktor Petrov. 
It has been tested with [Macaulay2 version 1.22](https://github.com/Macaulay2/M2/releases/tag/release-1.22), which includes `WeylGroups` package version 0.5.3.
It should also run with [Macaulay2 version 1.18](https://github.com/Macaulay2/M2/releases/tag/release-1.18) and `WeylGroups` package version 0.5.2.  Older versions of the `WeylGroups` package contain a bug that will render the results incorrect.  To see which version of the package you have, you can type `readPackage "WeylGroups"` in Macaulay2.  


## Executing the code

The computations necessary to complete the proof of [HZ, Proposition 3.3] can be run by executing all code in `main.m2`.  The code there uses functions from the `WeylGroups` package and the two small auxiliary packages `WeylGroupsExtra` and `Auxiliary` provided here.  The code is currently set up to peform computations for all Dynkin diagrams Σ of exceptional types (E6, E7, E8, G2, F4).  This can easily be changed by editing the very last line of `main.m2`.

The proof of [HZ, Proposition 3.3] only requires the verification of the condition *single cell* for certain pairs (Σ, H).  In addition, the code also checks whether the fix point moniod written as $\overline{\mathcal C}(H)^{[H]}$ in [HZ] is a *free* abelian monoid (in two different ways), and whether the pair (Σ, H) satisfies the weaker condition *orbit basis*.  However, because of the associated computational costs, the verification of this last condition is restricted to root systems Σ of rank < 8.  To change this behaviour, remove the conditional in the line

```
  if rank(R) < 8 then result#"orbitcondition" = checkIfLSatisfiesOrbitCondition(R,P);
```

in `main.m2`.

## Viewing the results

The results of the computations are written to tex files (`results_G2.tex`, `results_F4.tex`, ...).  To view them, compile the auxiliary file `ViewResults.tex` also provided here.  To display results for other than the exceptional types, the contents of `ViewResults.tex` need to be edited in an obvious way.  For reference, results for exceptional types are already included in the folder `results`.

## Interpretation

The notation in `ViewResults.tex` follows [HZ]. The numbering of simple roots follows the conventions of Nicolas Bourbaki, *Lie groups and Lie algebras 4-6* (see plates at the end of the book).  The conditions *single cell* and *orbit basis* are explained in [HZ, §1: Overview].  The condition *free* signifies whether the fixed point monoid $\overline{\mathcal C}(H)^{[H]}$ is free abelian.  Within (the comments in) the code, the fixed point monoid $\overline{\mathcal C}(H)^{[H]}$ is referred to as `FixedPointMonoid`.

