-- Authors: Nida Obatake, Elise Walker
-- Last Updated: Dec 12, 2021
-- This code is supplemental to the paper: 
--     "Newton-Okounkov bodies of chemical reaction systems"
--     (include arxiv link)

-- This contains supporting functions for computing the Newton-Okounkov body bound
--   of chemical reaction networks.

-- This code requires the following packages:
-- Polyhedra


computeNObodyVol = (subBasis) -> (
-- Input: subBasis, a subalgebra basis with grading variable the last variable in ring
-- Output: normalized volume of NObody
   semigroupGenerators = for i in subBasis list flatten exponents leadTerm i;
   semigroupGenerators = for i in semigroupGenerators list drop(i/(last i), -1);
   NewtonOkounkovBody = convexHull transpose matrix semigroupGenerators;
   return (numgens (ring first subBasis) -1 )! * (volume NewtonOkounkovBody)   
   )


computeIndexSNF = (sagbiBasis) -> (
-- Computes the index of S(V,nu) in ZZ^d+1, NOT the index S(V,nu)_0 in ZZ^d
-- If the index of S(V,nu) in ZZ^d+1 is 1, then the index of S(V,nu)_0 in ZZ^d is 1
-- Input: sagbiBasis, a subalgebra basis with grading variable the last variable in ring
-- Output: an integer, index of S(V,nu) in ZZ^d+1
   sagbiBasisVertices := for i in sagbiBasis list flatten exponents leadTerm i;
   sagbiBasisVertices = for i in sagbiBasisVertices list drop(i/(last i), -1);
   latticeGenerators := transpose matrix sagbiBasisVertices;
   D := smithNormalForm( latticeGenerators, ChangeMatrix=>{false, false});
   diagEntries := for i from 0 to min(numcols D, numrows D) - 1 list D_(i,i);
   nonzeroDiagEntries := select(diagEntries, i->not i == 0);
   return product nonzeroDiagEntries
   )

createTrinomials = (system) -> (
-- Input: system, a list of polynomials defining a chemical reaction system
-- Output: a list of trinomials such that system is in the span of the trinomials
--         note: polynomials in output may have 3 or fewer terms 
   monomialSystem := for i in system list terms i;
   trinomialV := {};
   newLast := {};
   for i in monomialSystem do (
      pairSystem := pack(i, 3);
      trinomialBasis := for j in pairSystem list sum(j);
      trinomialV = append(trinomialV, trinomialBasis);
      );   
   return flatten trinomialV   
   )

createBinomials = (system) -> (
-- Input: system, a list of polynomials defining a chemical reaction system
-- Output: a list of binomials such that system is in the span of the binomials
--         note: polynomials in output may have 2 or 1 term(s) 
   monomialSystem := for i in system list terms i;
   binomialV := {};
   newLast := {};
   for i in monomialSystem do (
      pairSystem := pack(i, 2);
      binomialBasis := for j in pairSystem list sum(j);
      binomialV = append(binomialV, binomialBasis);
      );   
   return flatten binomialV   
   )

makeCoef = (a, b) -> (
-- Input: a and b are positive integers
-- Output: returns one list of #a lists of length #b of random complex numbers
    apply(a,(j->(apply(b, (i ->(x:=random(sub(0,RR),sub(1,RR));cos(2*pi*x)+sin(2*pi*x)*ii))))))
    )

makeRandomSystem = (basisV) -> (
-- Input: basisV is a list of functions which are a basis for a vector space
-- Output: a random square system drawn from basisV, coefficients are complex
    numberEquations := length gens ring first basisV;
    coef := makeCoef(numberEquations, length basisV);
    randomSystem := for i in coef list sum apply(i, basisV, (j,k)-> j*k);
    return randomSystem 
    )

homToSameDegree = (polynomialList, s) -> (
-- Input: polynomialList, homogenization variable s where all inputs belong to same ring
---       s should not appear in polynomialList
-- Output: polynomials homogenized with s to the same degree (degree chosen minimally)
    homogList := for i in polynomialList list homogenize(i, s);
    maxDegree := max flatten for i in homogList list degree i;
    degToMax := for i in homogList list (maxDegree - (degree i)#0);
    return for i from 0 to #homogList -1 list s^(degToMax#i)*homogList#i
   )

