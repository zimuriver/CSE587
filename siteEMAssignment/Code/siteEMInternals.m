(* Put your code for sequencePosteriors, updateMotifPrior, and updatePFMCounts,
   along with any other functions you write to implement them, in this file.
   Don't forget that you can use your implementation of sitePosterior as part of your
   implementation of sequencePosteriors or, if you prefer, you can request a correct
   version from the TAs. *)
sequencePosteriors[sequence_, sitePrior_, siteProbs_, backgroundProbs_] :=
	Module[{posteriors,backgroundPrior,motifLength,i},
	   backgroundPrior = 1 - sitePrior;
	   motifLength = Length[siteProbs];
	   posteriors = Table[sitePosterior[Take[sequence,{i,i+motifLength-1}],sitePrior,backgroundPrior,siteProbs,backgroundProbs], {i ,Length[sequence] - motifLength+1}];
	  
	   posteriors
	]
	
updateMotifPrior[normalizedPosteriors_] := 
	Module[{motifCount,backgroundCount},
		
	    motifCount = Total[normalizedPosteriors,Infinity];
	   
		backgroundCount = Total[1-normalizedPosteriors,Infinity];
		
		motifCount/(motifCount+backgroundCount)
	]
	
updatePFMCounts[motifLength_, input_, normalizedPosteriors_, motifPseudocounts_, erasers_] :=
	Module[{currentErasers,sequence,i,j,k,newMotifCounts},
		
		newMotifCounts = Table[motifPseudocounts,motifLength];
		
		currentErasers = erasers;
		
	    If[erasers == False,currentErasers = input/input]; (*If no erasers, set them to 1*)
		Table[
	    	Table[
	  
	      sequence = Take[input[[j]],{i,i+motifLength-1}];
	    
	    
	      Do[newMotifCounts[[k,sequence[[k]]]] += (normalizedPosteriors[[j,i]]*currentErasers[[j,i+k-1]]),{k,motifLength}];
	    
	    
	   ,{i,Length[input[[j]]]-motifLength+1}]; 
	   			,{j,Length[input]}];
	   
	   Return[newMotifCounts];
	    
	]

(* Put your code for sitePosterior here. *)
sitePosterior[sequence_, sitePrior_, backgroundPrior_, siteProbs_, backgroundProbs_] := 
	Module[{probOfSites,probOfBackground},
 
 	
 	probOfSites = Map[siteProbs[[#]][[sequence[[#]]]]&, Range[Length[sequence]]];
 	probOfBackground = Map[backgroundProbs[[sequence[[#]]]]&,Range[Length[sequence]]];
 	
 	
 	(Apply[Times,probOfSites]*sitePrior)/(Apply[Times,probOfSites]*sitePrior + Apply[Times,probOfBackground]*backgroundPrior)
    
	  
]