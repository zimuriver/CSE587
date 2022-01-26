(* Mathematica package *)


(* This file contains code for two "hacks" described in Bailey & Elkan (1993), posterior
   window normalization and erasers. Neither of those are part of the EM algorithm. Posterior
   windown normalization is needed to prevent convergence on repetitive motifs like "GGGGGGGG".
   Erasers are needed to find more than one motif in a given set of sequences. Note that these
   hacks are needed because the probability model behind MM doesn't quite match the input.*)
   
(* normalizePosteriorWindows is the top level function for normalizing the posteriors according
   to Bailery and Elkin's algorithm. The inputs are motifLength and posteriors. posteriors is a 
   list containing one list for each input sequence, containing the posterior probabilities of
   motif occurrences at each position in the sequence. The length of the list is the length of the 
   sequence minus motifLength plus 1, since that is the maximum number of overlapping motif instances
   that can occur in a sequence.
*)
normalizePosteriorWindows[posteriors_, motifLength_]:=
	Map[rowNormalizePosteriorWindows[#, motifLength] &, posteriors]

(* rowNormalizePosteriorWindows divides the posterior for a motif starting at a given position
   by the sum of the posteriors for motifs that overlap that position, but only if the sum is
   > 1.0.*)
rowNormalizePosteriorWindows[posteriorsRow_, motifLength_]:=
	MapThread[If[#2 > 1.0,
		         #1 / #2,
		         #1] &, 
		      {posteriorsRow,
		       rowWindowMaxes[rowWindowSums[posteriorsRow, motifLength], motifLength]}]

(* For each entry in the list posteriorsRow, rowWindowSums calculates a sum of it and the previous motifLength-1 
   entries. If there are fewer than motifLength-1 previous entries, all existing previous entries are used. Put
   differently, it computes a running sum of elements, including at most motifLength elements in the sum.
   
   To calculate the running sum for a given position, it starts with the previous running sum and adds in 
   the one new element for the current position and subtracts out the one element of the previous sum that
   is not included in the current sum. This keeps the asymptotic running time constant, independent of
   motifLength. However, the savings is probably not large for a typical motifLength of 8. And doing it this
   way propogates any errors do to limited precision. On balance I think it's worth doing it this way, but
   the alternative of recalculating each time could be reasonable, too.
   *)
rowWindowSums[posteriorsRow_, motifLength_]:=
	Block[{sum=0.},
		Reap[Do[sum=sum+posteriorsRow[[i]];
		       If[i > motifLength,
		   	      sum = sum - posteriorsRow[[i - motifLength]]];
		       Sow[sum], 
		   {i, 1, Length[posteriorsRow]}]][[2,1]]
		]

(* For each entry in the list sumsRow, rowWindowMaxes calculates a max of it and the motifLength-1 subsequent
   entries. If there are fewer than motifLength-1 subsequent entries, all existing subsequent entries are used.
   
   To avoid doing work proportional to motifLength, a running max is kept and updated for each position. Most of
   the time this will only require constant work per position in the sumsRow, but occasionally the running max
   has to be recalculated from scratch.
   *)
rowWindowMaxes[sumsRow_, motifLength_]:=
	Module[{rowLength=Length[sumsRow],
		    max=0.,
		    maxesList={}},
		maxesList=Reap[
			Do[If[sumsRow[[i]] > max,
		      max = sumsRow[[i]],
		      (*sumsRow[[i]] is not the max of the new window, but if sumsRow[[i + motifLength]] != max
		        then you know that whatever was the max in the old window is max in the new window.*)
		      If[i + motifLength <= rowLength,
		      	 If[max == sumsRow[[i + motifLength]],
		      	    (*You don't know which element is the max after you drop sumsRow[[i + motifLength]]
		      	      from the window over which you are maxing, so recalculate the max from scratch*)
		      	    max = 0;
		      	    Do[If[sumsRow[[j]] > max,
		      	 	     max=sumsRow[[j]]],
		      	 	  {j, i, i + motifLength - 1 }]]]];
		       Sow[max],
		      {i, rowLength, 1, -1}]][[2,1]];
		Reverse[maxesList]]


updateErasers[erasers_, posteriors_, motifLength_]:=
	MapThread[rowUpdateErasers[#1, #2, motifLength] &, {erasers, posteriors}]

(* The length of an erasers row is the same as the length of the input row (and greater than the length of the
   corresponding posteriors row by the motifLength-1.
   
   To avoid doing work proportional to motifLength at every step, a running product is kept and updated at 
   each step by multiplying by a new element and dividing by an old element. See the comment for rowWindowSums.
   The exception is when you would encounter a divide by zero (and hence the current product has a zero factor
   and is therefore equal to zero). The the running product is multiplied by a zero factor, all the information
   that would allow incremental calculation is wiped out and the product must be recalculated from scratch
   whenever the zero factor drops out of the window.
*)
   
rowUpdateErasers[erasersRow_, posteriorsRow_, motifLength_]:=
	Module[{posteriorsLength=Length[posteriorsRow],
		    product=1.},
		Reap[
		Do[If[i<=posteriorsLength,
			  product=product*(1.-posteriorsRow[[i]])];
		   If[i > motifLength,
		   	 If[posteriorsRow[[i - motifLength]]!=1.0,
    		    product = product / (1. - posteriorsRow[[i - motifLength]]),
    		    (*If the factor you're trying to divide out was zero then you have to recalculate the
    		      product for the entire new window.*)
    		    product=1.0;
    		    Do[product=product*(1. - posteriorsRow[[j]]),
    		    	{j, i-motifLength+1, Min[i, posteriorsLength]}]]];
		   Sow[erasersRow[[i]]*product], 
		   {i, 1, Length[erasersRow]}]][[2, 1]]
	]