 sampleLabels[binCounts_, oldLabels_,thetas_, hyperParams_] :=
 Module[{newLabels,C1,P1,C2,P2,faceProduct1,faceProduct2,i,j},
 
 newLabels = oldLabels;




 Do[
 	C1 = BinCounts[newLabels,{1,3}][[1]];
 	C2 = BinCounts[newLabels,{1,3}][[2]];
 	
 	faceProduct1 = Times@@(Table[safeExponentiate[thetas[[1]][[j]],binCounts[[i]][[j]]],{j,Length[thetas[[1]]]}]);
 	
 	P1 = ((C1 + hyperParams[[1]] - 1)/(C1+ C2 + hyperParams[[1]] + hyperParams[[2]] - 1))*faceProduct1; 

	faceProduct2=Times@@(Table[safeExponentiate[thetas[[2]][[j]],binCounts[[i]][[j]]],{j,Length[thetas[[2]]]}]);
	
	P2 = ((C2 + hyperParams[[2]] - 1)/(C1 + C2 + hyperParams[[1]] + hyperParams[[2]] - 1))*faceProduct2;
	
	
 	newLabels[[i]] = 1 + RandomVariate[BernoulliDistribution[P2/(P1+P2)]];		
 				,{i,Length[oldLabels]}];
 
 Return[newLabels];
 ]
 	 