-- Authors: Nida Obatake, Elise Walker
-- Last Updated: Jan 2, 2022
-- This code is supplemental to the paper:
--     "Newton-Okounkov bodies of chemical reaction systems"
--     


-- Example 5.3
-- Example 6.5.3 from Perez Millan's Thesis
-- This computes the birationally invariant self-intersection index for Example 6.5.3

------------------------------------------------------------------------------------------
restart
needsPackage "SubalgebraBases"
needsPackage "Polyhedra"
needsPackage "Cremona"
needsPackage "PHCpack"
needs "./Procedures.m2"

-- Chemical reaction network in Example 6.5.3, where 's' is a grading variable
------------------------------------------------------------------------------
R = QQ[x_1, x_2, s]
system = {x_1*((x_1-1)^2+(x_2-2)^2)*s, x_2*((x_1-1)^2+(x_2-2)^2)*s}


-----Computing a Newton-Okounkov Body-----
------------------------------------------
-- Choose a vector space basis such that each polynomial 'system' is in the vector space
-- Here, we select trinomials from each polynomial in 'system'
basisV = createTrinomials system --basis has 4 elements

-- In this case, 'basisV' is a SAGBI basis already, which we verify in the following:
(verifySagbi basisV)#"isSAGBI"  

-- Next we compute the Newton-Okounkov body as the convex hull of the semigroup generators
semigroupGenerators = for i in basisV list flatten exponents leadTerm i
semigroupGenerators = for i in semigroupGenerators list drop(i/(last i), -1)
P = convexHull transpose matrix semigroupGenerators

-- The volume of the Newton-Okounkov body P is 1
NObodyVol = volume P


-----Computing the Degree of the Kodaira Map-----
-------------------------------------------------
-- Homogenize polynomials so grading is captured by total degree
dehomBasisV = for i in basisV list sub(i, s=>1_R);
basisVHomToSameDegree = homToSameDegree(dehomBasisV, s)
-- Create Kodaira map
S = QQ[y_1..y_(#dehomBasisV)]
phi = map(R, S, basisVHomToSameDegree)
kodairaDegree = degreeMap phi --returns 2

-----Computing the Index of the Lattice-----
--------------------------------------------
latticeIndex = computeIndexSNF(basisV) -- returns 1

-----Computing the Birationally Invariant Self-Intersection Index-----
----------------------------------------------------------------------
2! * kodairaDegree * NObodyVol / latticeIndex -- biii of 4

-----Computing the Base Locus-----
----------------------------------
complexR = CC[drop(gens R, -1)];
complexBasisV = for i in basisV list sub(sub(i, s=>1_R), complexR);
baseLocus = solveSystem(complexBasisV)
#baseLocus --3 solutions


-----Comparing to Mixed Volume-----
-----------------------------------
-- The mixed volume of a general system pulled from vector space V is:
generalSystem = makeRandomSystem(complexBasisV);
mixedVolume(generalSystem) --returns 8
-- The mixed volume of our original system is:
mixedVolume(for i in system list sub(sub(i, s=>1_R), complexR)) --returns 4


------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------



------------------------------------------------------------------------------------------
-- There are several different choices for the vector space, and different choices give
-- different Newton-Okounkov bodies, as seen below.
------------------------------------------------------------------------------------------



------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------
restart
needsPackage "SubalgebraBases"
needsPackage "Polyhedra"
needsPackage "Cremona"
needsPackage "PHCpack"
needs "./Procedures.m2"

-- Chemical reaction network in Example 6.5.3, where 's' is a grading variable
------------------------------------------------------------------------------
R = QQ[x_1, x_2, s]
system = {x_1*((x_1-1)^2+(x_2-2)^2)*s, x_2*((x_1-1)^2+(x_2-2)^2)*s}


-----Computing a Newton-Okounkov Body-----
------------------------------------------
-- Choose a vector space basis such that each polynomial 'system' is in the vector space
-- Here, we select binomials from each polynomial in 'system'
basisV = createBinomials system 

-- In this case, 'basisV' is not a SAGBI basis already, which we verify in the following:
(verifySagbi basisV)#"isSAGBI" 
-- We compute a SAGBI basis by the following:
Sag = sagbi(basisV, Limit=> 20) --computes SAGBI basis
SB = flatten entries gens Sag   -- SAGBI basis has 8 elements 

-- Next we compute the Newton-Okounkov body as the convex hull of the semigroup generators
semigroupGenerators = for i in SB list flatten exponents leadTerm i
semigroupGenerators = for i in semigroupGenerators list drop(i/(last i), -1)
P = convexHull transpose matrix semigroupGenerators

-- The volume of the Newton-Okounkov body P is 1
NObodyVol = volume P


-----Computing the Degree of the Kodaira Map-----
-------------------------------------------------
-- Homogenize polynomials so grading is captured by total degree
dehomBasisV = for i in basisV list sub(i, s=>1_R);
basisVHomToSameDegree = homToSameDegree(dehomBasisV, s)
-- Create Kodaira map
S = QQ[y_1..y_(#dehomBasisV)]
phi = map(R, S, basisVHomToSameDegree)
kodairaDegree = degreeMap phi --returns 1


-----Computing the Birationally Invariant Self-Intersection Index-----
----------------------------------------------------------------------
2! * kodairaDegree * NObodyVol -- biii of 6


-----Computing the Base Locus-----
----------------------------------
complexR = CC[drop(gens R, -1)];
complexBasisV = for i in basisV list sub(sub(i, s=>1_R), complexR);
baseLocus = solveSystem(complexBasisV)
#baseLocus --1 solution, (0,0)


-----Comparing to Mixed Volume-----
-----------------------------------
-- The mixed volume of a general system pulled from vector space V is:
generalSystem = makeRandomSystem(complexBasisV);
mixedVolume(generalSystem) --returns 8
-- The mixed volume of our original system is:
mixedVolume(for i in system list sub(sub(i, s=>1_R), complexR)) --returns 4
