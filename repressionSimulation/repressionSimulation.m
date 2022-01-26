(* Wolfram Language package *)

nSolveRegVect[synthRtVect_, repressionMatrix_, mDegRtConstVect_, translationRtConstVect_, pDegRtConstVect_, initMRNAConcsVect_, initProteinConcsVect_] := 
 Module[{numGenes = Length[synthRtVect], eqns, funs, solns,i,j},
  (* Use Table to define a list of numGenes equations for m' and numGenes equations for p'. Set 'eqns' to be the list of equations.*)
  (* Also use Table to list the functions you want NDSolve to solve for. See docs on NDSolve. set 'funs' to be list of functions.*)
  (* My implementation of the missing code is about 10 lines, including extra linebreaks to make it readable. Your formatting may vary,
     but this gives you a rough idea of the scope of the assignment. *)
 Clear[m];Clear[p];

 eqns = Table[{
 	m[i]'[t] == synthRtVect[[i]]*Product[1/(1+repressionMatrix[[j,i]]*p[j][t]),{j,numGenes}]-mDegRtConstVect[[i]]*m[i][t],
	m[i][0]==initMRNAConcsVect[[i]],
	p[i]'[t] == translationRtConstVect[[i]]*m[i][t] - pDegRtConstVect[[i]]*p[i][t],
	p[i][0]==initProteinConcsVect[[i]] 
 }
 ,{i,numGenes}];

 eqns = Flatten[eqns];
 	funs = Table[{m[i],p[i]},{i,1,numGenes}];
 	funs = Flatten[funs];
  solns = NDSolve[eqns, funs, {t, 0, 40}];
  Print[eqns];
  Print[Plot[Evaluate[Flatten[Table[{m[i][t], p[i][t]}, {i, numGenes}], 1] 
  	                          /. solns], 
  	        {t, 0, 10}, 
  	        PlotRange -> All]];
  solns]