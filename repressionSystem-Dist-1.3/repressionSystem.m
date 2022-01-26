(* Mathematica Raw Program *)

generateRandomParams[regulatorsMask_, seed_]:=
	Module[{bVector,cMatrix,numGenes,i,j},
		
		numGenes = Length[regulatorsMask];
		bVector = ConstantArray[0.0,numGenes];
		cMatrix = ConstantArray[0.0,{numGenes,numGenes}];
		SeedRandom[seed];
		Do[
			
			bVector[[i]] = RandomReal[{.1,1}];
				Do[
					If[i != j && regulatorsMask[[i]] == 1,cMatrix[[i,j]] = RandomReal[]];
		
				,{j,numGenes}];
		,{i,numGenes}];

	{bVector,cMatrix}
	]
(* genericEqns generates system of equations for n mutually repressing proteins. They are
   generic in that the parameters remain as indexed symbols b[j] and c[j,i], except that
   c[i,i] is replaced by zero, implementing the ban on auto-repression.*)
genericEqns[n_]:=
	Module[{eqns=ConstantArray[0,n],i,j},
	Clear[b];Clear[m];Clear[c];
		Do[
		eqns[[i]] = (b[i]m[i]Product[If[i!=j,m[j]c[j,i]+1,1],{j,n}]-1.0==0);
			
		,{i,n}];
	eqns
	]
(* mutantEqns takes a genePresenceMask indicating which genes were deleted in the cells 
   in which these mRNA measurements were taken. It generates a set of generic equations 
   and then replaces the mRNA level variables for the deleted genes with zeros, ensuring 
   that any actual measurments made on those deleted genes are ignored and any solver 
   doesn't have to deal with them.
*)
mutantEqns[genePresenceMask_]:=
Module[{eqns,numGenes = Length[genePresenceMask],i,newEqns},
	eqns = genericEqns[numGenes];
	newEqns = {};
	Do[
		If[genePresenceMask[[i]] == 0, m[i] = 0, AppendTo[newEqns,eqns[[i]]]];
	
	,{i,numGenes}];

newEqns
]

(* mutantExpVars returns a list of variables representing the expression levels of
   genes that remain in the genotype described by genePresenceMask. *)
mutantExpVars[genePresenceMask_]:=
Module[
{i,vars = {}}, 
Do[
	If[genePresenceMask[[i]] == 1,AppendTo[vars,m[i]]];
	,{i,Length[genePresenceMask]}];
	vars
]
(* expressionEqns replaces the symbolic parameters of the generic equations with actual values
   from a parameter matrix, returning a final set of equations with fixed parameters but variables
   for mRNA levels. Solving these for the mRNA variables yields a set of expression levels that
   are consistent with the parameters.*)
expressionEqns[eqns_, {bVector_, cMatrix_}]:=
Module[{i,j},
	Do[
		b[i] = bVector[[i]];
		Do[
			c[j,i] = cMatrix[[j,i]]
		,{j,Length[cMatrix[[i]]]}]
	,{i,Length[bVector]}];
	eqns
]
(* expressionProfile finds a set of steady state mRNA expression levels that are consistent with
   a set of expression equations and a genotype as indicated by genePresenceMask. To do this, you
   will use the built in function FindRoot to find a solution to the set of simultaneous equations
   you have constructed with the function expressionEqns. IMPORTANT: FindRoot allows you to specify
   a starting, minimum, and maximum value for each unknown (in this case the mRNA levels). You need
   to do this, making use of the bounds specified in the Introduction section of the assignment notebook. 
   
   The result of FindRoot is a list of replacement rules for expression variables. These are then used 
   to construct a vector of length n containing the expression levels for genes that are present and 
   zeros for genes that are deleted.*)
expressionProfile[expressionEqns_, genePresenceMask_]:=
Module[{rules={},roots},
	
	Map[(AppendTo[rules,{#,1,0,100}])&, mutantExpVars[genePresenceMask]];
	
	roots =	FindRoot[expressionEqns,rules];
	
	roots
]
(* expressionMatrix takes parameters and a list of genePresenceMasks, one for each strain that we
   want an expression profile for. For each mask, it makes a set of equations and calls
   expressionProfile with those equations and the mask. This returns a list of rules that is
   applied to the expression level variables to produce actual values, which are returned.
   The return value is a matrix in which each row represents the expression levels for one
   strain (as described by one genePresenceMask). *)
expressionMatrix[params_, genePresenceMasks_]:=
	(*This ensures that b, c, m are undefined here and in all functions called from here.*)
Module[{i,matrix = {},mVals,possibleVals},	
	Block[{b, c, m},
		Map[(mVals = expressionProfile[expressionEqns[mutantEqns[#], params],#];
			possibleVals = Table[m[i],{i,Length[#]}];
			
			AppendTo[matrix,possibleVals/.mVals])&
			
			,genePresenceMasks];	
		
	];
	matrix
]	
(*Part 4: Parameter inference.*)
	
(* inferParameters takes a list of expression profiles and a 
   matched list of genePresenceMasks, one for each strain that we have an expression 
   profile for. For each profile-mask pair, it makes a set of equations with the c parameters
   as unknown variables. It is important that the the right-hand-side of each equation be zero.
   A vector of the left-hand-sides of these equations is then dotted with itself to make a
   sum of squared errors, where the error is the deviation from zero. NMminize is then used to 
   find values of the c[i,j] parameters that are consistent with the equations, or as close to
   consistent as possible. We must use NMiminize instead of NSolve because the number of 
   equations may be larger than the number of variables and NSolve won't accept over-determined 
   system of equations *)

(* If the level of regulator j is 0 in all experiments then c[j,i] will never appear in the equations 
   and therefore can't be included in the list of variables to solve for.*)

inferParameters[regulatorsMask_, geneExpressionProfiles_, genePresenceMasks_]:=
	(*Make sure b, c, m are undefined here and in all functions called from here.*)
	Module[{paramsList,eqns,i,j,Bs,Cs,bVals,cVals,newEqns,sumLeft,ans,allB,allC},
	Block[{b, c, m},
	paramsList = generateRandomParams[regulatorsMask,150];
	Bs = paramsList[[1]];
	Cs = paramsList[[2]];
	(*get b and c values that are going to need to be found by nMinimize*)
	bVals = Table[If[Bs[[i]]!=0,b[i],##&[]],{i,Length[Bs]}];
	cVals = Table[If[Cs[[i,j]]!=0,c[i,j],##&[]],{i,Length[Cs]},{j,Length[Cs[[i]]]}];
	
	eqns = Table[parameterEqns[mutantEqns[genePresenceMasks[[i]]],geneExpressionProfiles[[i]]],{i,Length[geneExpressionProfiles]}];

	newEqns = eqns;
	(*get rid of == 0*)
	Do[newEqns[[i,j]] = eqns[[i,j,1]],{i,Length[eqns]},{j,Length[eqns[[i]]]}];
	
	(*sum the squared errors of left sides of eqns*)
	sumLeft = Total[Map[parameterSumSquaredError[#]&,newEqns]];
	Print[sumLeft];
	ans = NMinimize[{sumLeft,
		Apply[And,Map[1 >= # >= 0 &,Flatten[{bVals,cVals}]]]},
		
		Flatten[{bVals,cVals}]];
	
	allB = Table[If[Bs[[i]]!=0,b[i],0.0],{i,Length[Bs]}];
	allC = Table[If[Cs[[i,j]]!=0,c[i,j],0.0],{i,Length[Cs]},{j,Length[Cs[[i]]]}];

	 {allB,allC}/.ans[[2]]
	 ]
	 
	]
(* parameterEqns replaces the symbolic mRNA expression levels of the generic equations with actual values
   from a gene expression matrix, returning a final set of equations with fixed mRNA levels but variables
   for parameters. Solving these for the parameter variables yields a set of parameters that
   are consistent with the expression levels.*)

parameterEqns[eqns_, geneExpressionProfile_]:=
Module[{i},
Do[m[i] = geneExpressionProfile[[i]],{i,Length[geneExpressionProfile]}];

eqns

]
parameterSumSquaredError[eqns_]:=
Module[{},
Dot[eqns,eqns]



]
(*SECTION Some auxiliary functions to support testing parameter inference.*)

(* Generate a random set of parameters, generate expression profiles for those 
   parameters in wild-type cells and cells with all genes deleted one at a time.
   Then use these expression levels to infer the parameters. Compare the inferred 
   parameters to the ones used to generate the expression levels. If the 
   difference is less than or equal to 10^-4 replace it by zero. Return the resulting matrix of differences. 
   
   It should contain nothing but zeros.*)		 

testParameterInference[n_, seed_]:=
	With[{regulatorsMask=ConstantArray[1, {n}],
		  genePresenceMasks=allSinglesAndWTPresenceMatrix[n]},
		With[{params=generateRandomParams[regulatorsMask, seed]},
			With[{expressionMatrix=expressionMatrix[params, genePresenceMasks]},
				With[{diff=params-inferParameters[regulatorsMask,
					 	                  		  expressionMatrix,
					 	     			          genePresenceMasks]},
				Round[diff, 0.0001]
				]]]]

allSinglesAndWTPresenceMatrix[n_]:=
	Join[ConstantArray[1,{1,n}],
		 ConstantArray[1, {n,n}] - IdentityMatrix[n]]
		 
