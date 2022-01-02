-- Authors: Nida Obatake, Elise Walker
-- Last Updated: Jan 2, 2022
-- This code is supplemental to the paper:
--     "Newton-Okounkov bodies of chemical reaction systems"
--     

-- Example 5.4
-- From Example 4.2 in Boros, Craciun, Yu (2019)
-- This computes the birationally invariant self-intersection index for Example 4.2

------------------------------------------------------------------------------------------
restart
needsPackage "SubalgebraBases"
needsPackage "Polyhedra"
needsPackage "Cremona"
needsPackage "PHCpack"
needs "./Procedures.m2"

-- Chemical reaction network in Example 4.2
-------------------------------------------
R = QQ[x, y, s]
p1 = x^2*y^2 + x^2 + y^2 + 1 -5*x*y;
p2 = 1;
p3 = x*y;
p4 = x;
p5 = y;

system = {p1*(p2-p3), p1*(p4-p5)}

-----Computing a Newton-Okounkov Body-----
------------------------------------------
-- Choose a vector space basis such that each polynomial 'system' is in the vector space
-- Here, 's' is a grading variable
basisV = {s*p1*p2, s*p1*p3, s*p1*p4, s*p1*p5}

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
S = QQ[z_1..z_(#dehomBasisV)]
phi = map(R, S, basisVHomToSameDegree)
kodairaDegree = degreeMap phi --returns 1

-----Computing the Index of the Lattice-----
--------------------------------------------
latticeIndex = computeIndexSNF(basisV) -- returns 1


-----Computing the Birationally Invariant Self-Intersection Index-----
----------------------------------------------------------------------
2! * kodairaDegree * NObodyVol / latticeIndex -- biii of 2


-----Solving a General System-----
----------------------------------
complexR = CC[drop(gens R, -1)]
complexBasisV = for i in basisV list sub(sub(i, s=>1_R), complexR)
generalSystem = makeRandomSystem(complexBasisV);
sol = solveSystem(generalSystem) 
#sol 


-----Computing the Base Locus-----
----------------------------------
baseLocus = solveSystem(complexBasisV)
#baseLocus 


-----Comparing to Mixed Volume-----
-----------------------------------
-- The mixed volume of a general system pulled from vector space V is:
generalSystem = makeRandomSystem(complexBasisV);
mixedVolume(generalSystem) --returns 18
-- The mixed volume of our original system is:
mixedVolume(for i in system list sub(sub(i, s=>1_R), complexR)) --returns 18


------------------------------------------------------------------------------------------
