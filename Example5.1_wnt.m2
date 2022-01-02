-- Authors: Nida Obatake, Elise Walker
-- Last Updated: Jan 2, 2022
-- This code is supplemental to the paper: 
--     "Newton-Okounkov bodies of chemical reaction systems"
--     

-- Motivating Example, Summarized in Example 5.1
-- Gross, Harrington, Rosen, Sturmfels
-- https://arxiv.org/pdf/1502.03188.pdf
-- This computes the birationally invariant self-intersection index for WNT pathway model

needsPackage "SubalgebraBasesDevo"
needsPackage "Polyhedra"
needsPackage "PHCpack"
needsPackage "Cremona"
needsPackage "Bertini"
needs "./Procedures.m2"


--Step 0: Choose monomial term order. We use the default of grvlex
-- WNT pathway model with random rational coefficients and conservation laws
for i in 1..31 do k_i = random(QQ, Height => 1000); 
for i in 1..5  do c_i = random(QQ, Height => 1000);
R = QQ[x_1..x_19]; --default term order of grvlex
system = {-k_3*x_2*x_4+k_1*x_1+(-k_2-k_26)*x_2+k_27*x_3+(k_4+k_5)*x_14,
          -k_14*x_3*x_6+k_26*x_2-k_27*x_3+(k_15+k_16)*x_15,
          -k_6*x_5*x_8-k_28*x_5+k_29*x_7+k_5*x_14+k_7*x_16,
          -k_14*x_3*x_6-k_20*x_6*x_11+k_15*x_15+k_19*x_17+(k_21+k_22)*x_19,
          -k_17*x_7*x_9+k_28*x_5-k_29*x_7+k_16*x_15+k_18*x_17,
          -k_9*x_4*x_10+(-k_13-k_30)*x_10+k_31*x_11+k_10*x_18+k_12,
          -k_20*x_6*x_11-k_24*x_11*x_12+k_30*x_10+(-k_23-k_31)*x_11+k_25*x_13+k_21*x_19,
          k_24*x_11*x_12-k_25*x_13,
          k_3*x_2*x_4+(-k_4-k_5)*x_14,
          k_14*x_3*x_6+(-k_15-k_16)*x_15,
          k_6*x_5*x_8+(-k_7-k_8)*x_16,
          k_17*x_7*x_9+(-k_18-k_19)*x_17,
          k_9*x_4*x_10+(-k_10-k_11)*x_18,
          k_20*x_6*x_11+(-k_21-k_22)*x_19,
          x_1+x_2+x_3+x_14+x_15-c_1,
          x_4+x_5+x_6+x_7+x_14+x_15+x_16+x_17+x_18+x_19-c_2,
          x_8+x_16-c_3,
          x_9+x_17-c_4,
          x_12+x_13-c_5
          };

-- Step 1:
-- Choose a vector space basis s.t. each polynomial from 'system' is in the vector space
basisV = {-k_3*x_2*x_4+ (-k_2-k_26)*x_2+ (k_4+k_5)*x_14,
          k_1*x_1+k_27*x_3,
          -k_14*x_3*x_6+ -k_27*x_3 +(k_15+k_16)*x_15,
          k_26*x_2,
          -k_6*x_5*x_8-k_28*x_5+k_7*x_16,
          k_29*x_7+k_5*x_14,
          -k_14*x_3*x_6-k_20*x_6*x_11+k_15*x_15,
          k_19*x_17+(k_21+k_22)*x_19,
          -k_17*x_7*x_9-k_29*x_7+k_18*x_17,
          k_28*x_5+k_16*x_15,
          -k_9*x_4*x_10+(-k_13-k_30)*x_10+k_10*x_18,
          k_31*x_11+k_12,
          (-k_23-k_31)*x_11-k_24*x_11*x_12+k_25*x_13,
          -k_20*x_6*x_11+k_21*x_19+k_30*x_10,
          k_24*x_11*x_12-k_25*x_13,
          k_3*x_2*x_4+(-k_4-k_5)*x_14,
          k_14*x_3*x_6+(-k_15-k_16)*x_15,
          k_6*x_5*x_8+(-k_7-k_8)*x_16,
          k_17*x_7*x_9+(-k_18-k_19)*x_17,
          k_9*x_4*x_10+(-k_10-k_11)*x_18,
          k_20*x_6*x_11+(-k_21-k_22)*x_19,
          x_2+x_14,
          x_3+x_15,
          x_1-c_1,
          x_4+x_6+x_15+x_19+x_14+x_18,
          x_5+x_16+x_7+x_17-c_2,
          x_8+x_16-c_3,
          x_9+x_17-c_4,
          x_12+x_13-c_5
          };

-- Step 3: Compute a finite Khovanksii basis for R_V
-- first homogenize as in R_V
Rs = QQ[gens R |{s}];
basisVs = for i in basisV list sub(i, Rs)*s;
-- Check if basis is a SAGBI basis
(verifySagbi basisVs)#"isSAGBI" --returns false 
SB = flatten entries gens sagbi(basisVs, Limit=>20) --32 generators, 6 binom, 6 six-term
(verifySagbi SB)#"isSAGBI" --returns true 


-- Step 6:
-----Computing a Newton-Okounkov Body-----
------------------------------------------
NObodyVol = computeNObodyVol(SB) --returns 32 for the normalized volume of the NObody


-- Step 4:
-----Computing the Degree of the Kodaira Map-----
-------------------------------------------------
dehomBasisV = for i in basisVs list sub(i, s=>1_Rs);
basisVHomToSameDegree = homToSameDegree(dehomBasisV, s);
-- Create Kodaira map
S = QQ[y_1..y_(#dehomBasisV)];
phi = map(Rs, S, basisVHomToSameDegree);
kodairaDegree = degreeMap phi --returns 1


-- Step 5: 
-----Computing the Index of the Lattice-----
--------------------------------------------
latticeIndex = computeIndexSNF(SB) -- returns 1


-- Step 7:
-----Computing the Birationally Invariant Self-Intersection Index-----
----------------------------------------------------------------------
kodairaDegree * NObodyVol / latticeIndex -- biii of WNT network, returns 32


----------------------------------------------------------------------
----------------------------------------------------------------------
Discussion
----------------------------------------------------------------------
----------------------------------------------------------------------


-----Solving a General System-----
----------------------------------
-- The intersection index should be the number of solutions to a general 
-- system drawn from V.
-- This computes the number of solutions to a general system drawn from V, 
-- which is 32. This is equal to the intersection index, as expected.
--
complexR = CC[drop(gens Rs, -1)]
complexBasisV = for i in dehomBasisV list sub(i, complexR)
generalSystem = makeRandomSystem(complexBasisV);
sol = solveSystem(generalSystem) 
#sol --32 solutions


-----Computing the Base Locus-----
----------------------------------
-- One should check there are no real positive solutions in the base locus of V.
-- By inspection, the base locus for V is empty for general choices of coefficients
baseLocus = bertiniZeroDimSolve(complexBasisV_{11,12,14,28}, AffVariableGroup=>{x_11, x_12, x_13})
#baseLocus --empty


-----Comparing to Mixed Volume-----
-----------------------------------
-- The intersection index is tighter than both the mixed volume of a general system
-- drawn from V, as well as the mixed volume of the original chemical reaction system.
-- The mixed volume of a general system pulled from vector space V is:
mixedVolume(generalSystem) --returns 60
-- The mixed volume of our original chemical reaction system is:
mixedVolume(for i in system list sub(i, complexR)) --returns 56
