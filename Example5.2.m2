-- Authors: Nida Obatake, Elise Walker
-- Last Updated: Jan 1, 2022
-- This code is supplemental to the paper: 
--     "Newton-Okounkov bodies of chemical reaction systems"
--     

-- Example 5.2
-- This computes the birationally invariant self-intersection index for a simple network

needsPackage "SubalgebraBasesDevo"
needsPackage "Polyhedra"
needsPackage "PHCpack"
needsPackage "Cremona"
needsPackage "Bertini"
needs "./Procedures.m2"

-- Step 0: 
-- Choose monomial term order. We use the default of grvlex
for i in 1..5 do k_i = random(QQ, Height => 1000); 
R = QQ[x,y]
system = { k_1 * (x^2 + y^2) - k_2 * x + k_3,
           k_1 * (x^2 + y^2) - k_4 * y + k_5 };

-- Step 1:
-- Choose a vector space basis s.t. each polynomial from 'system' is in the vector space
basisV = { 1 , x, y, x^2 + y^2};

-- Step 2: 
-- Compute a finite Khovanskii basis for R_V
-- In this case, basisV is already a finite Khovanskii basis, verified with verifySagbi
Rs = QQ[gens R | {s}]
basisVs = for i in basisV list sub(i, Rs)*s;
(verifySagbi basisVs)#"isSAGBI" --returns true

-- Step 3: 
-- Compute the normalized volume of the Newton-Okounkov body associated to R_V
normalizedNObodyVolume = computeNObodyVol(basisVs)

-- Step 4:
-- Compute the degree of the Kodaira map
-- We use the degreeMap function from the Cremona.m2, which requires a homogeneous ring.
dehomBasisV = for i in basisVs list sub(i, s=>1_Rs);
basisVHomToSameDegree = homToSameDegree(dehomBasisV, s);
S = QQ[z_1..z_(#dehomBasisV)];
kodairaMap = map(Rs, S, basisVHomToSameDegree);
kodairaDegree = degreeMap kodairaMap --returns 1

-- Step 5:
-- Compute the index of the lattice 
-- We show that the index of the entire semigroup is 1, which shows the index of G_0 is 1
latticeIndex = computeIndexSNF(basisVs) --returns 1

-- Step 6:
-- Compute the birationally invariant self-intersection index
kodairaDegree * normalizedNObodyVolume / latticeIndex -- biii of network, returns 2


------------------------------------------------------------------------------------------
Comparing to Mixed Volume
------------------------------------------------------------------------------------------
-- Comparing to mixed volume of a general system from V
CCR = CC[gens R];
CCbasisV = for i in basisV list sub(i, CCR);
randomSystemFromV = makeRandomSystem(CCbasisV);
mixedVolume(randomSystemFromV) --returns 4

-- Comparing to mixed volume of the original system
CCsystem = for i in system list sub(i, CCR);
mixedVolume(CCsystem)  --returns 4
