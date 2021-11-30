
The code in this repository complements the following research manuscript:

[HZ] Tobias Hemmert & Marcus Zibrowius, *The Witt rings of many flag varieties are exterior algebras*

## Prerequisites

The code is written in Macaulay2 (version 1.18) and requires the `WeylGroups` package, version >= 0.5.2, written by Baptiste Calmès and Viktor Petrov.  To see which version of the package you have, you can type `readPackage "WeylGroups"` in Macaulay2.  Version 0.5.2 was added to the Macaulay2 repositories in November 2021.

## Executing the code

The computations necessary to complete the proof of [HZ, Proposition 3.3] can be run by executing all code in `main.m2`.  The code there uses functions from the `WeylGroups` package and the two small auxiliary packages `WeylGroupsExtra` and `Auxiliary` provided here.  The code is currently set up to peform computations for all Dynkin diagrams Σ of exceptional types (E6, E7, E8, G2, F4).  This can easily be changed by editing the very last line of `main.m2`.

## Viewing the results

The results of the computations are written to tex files (`results_G2.tex`, `results_F4.tex`, ...).  To view them, compile the auxiliary file `ViewResults.tex` also provided here.  To display results for other than the exceptional types, the contents of `ViewResults.tex` need to be edited in an obvious way.  For reference, results for types A3, A5, B5, C5, D5, E6, E7, E8, F4, G2 are already included in the folder `results`.

## Interpretation

The notation in `ViewResults.tex` follows [HZ]. The numbering of simple roots follows the conventions of Nicolas Bourbaki, *Lie groups and Lie algebras 4-6* (see plates at the end of the book).  The conditions *single cell* and *orbit basis* are explained in [HZ, §1: Overview].  The condition *free* signifies whether the fixed point monoid $\overline{\mathcal C}(H)^{[H]}$ is free abelian.

