


The code in this repository complements the following research manuscript:

[HZ] Tobias Hemmert & Marcus Zibrowius, *The Witt rings of many flag varieties are exterior algebras*

It is written in Macaulay2 (version 1.18) and requires the `WeylGroups` package (version >= 0.5.2) written by Baptiste Calmès and Viktor Petrov.  To see which version of the package you have, you can type `readPackage "WeylGroups"` in Macaulay2.

The computations necessary to complete the proof of [HZ, Proposition 3.3] can be run by executing the code on `main.m2`.  The code there uses functions from the `WeylGroups` package and dthe two auxiliary packages `WeylGroupsExtra` and `Auxiliary` provided here.  The code is currently set up to peform computations of Dynkin diagrams Σ of exceptional types (G_2, F_4, E_6, E_7, E_8).  This can easily be changed by editing the very last line of `main.m2`.










