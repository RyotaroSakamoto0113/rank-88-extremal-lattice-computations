# Exact number-field computations

The class-number and unit-group calculations use PARI/GP `bnfcertify`. The reproduction driver runs `compute/L/src/field_arithmetic.gp`, `small_ideals.gp`, and `nebe_class_number.gp`. Their recorded outputs are in the separate certificate archive under `L/results/fields/`, `L/results/small_ideals/`, and `L/results/nebe_class_number/`.

|Field|Computed result|Defining data|
|---|---|---|
|F=Q(sqrt(-23))|class group C3; units ±1|polynomial x²−x+6|
|K=Q(zeta23)|class group C3|degree22, discriminant −23²¹, polynomial Phi23|
|Kplus=Q(zeta23+zeta23^-1)|class number1; unit-signature rank11, narrow class number1|degree11, discriminant23¹⁰; real cyclotomic field|
|Q(sqrt(-23),zeta5)|class group C3; relative unit norm onto totally positive units|degree8, discriminant23⁴5⁶; its real subfield has unit-signature rank3 and norm-image rank1 modulo squares|
|Q(sqrt(-7),zeta13)|class group C7|degree24; real subfield class number1; Hasse unit index2; comparison in Remark 3.15|

A class group or a list of fundamental units is not itself a statement of relative unit-norm surjectivity. That inference also needs a norm/signature argument beyond the class-group and unit-group output. Explicit cyclotomic ideal identities needed for the lattice construction, and the concrete semilinear beta satisfying `2 beta bar(beta)=1`, also have an independent direct check in `verify_beta.gp`, in addition to the full class-group and unit computations.
