(* Mathematica Raw Program *)

(* MM sets things up and calls MMIteration once for each motif to be found. If more than one motif is to be found,
   MM takes care of updating the erasers after each call to MMIteration. Whether the number of motifs to be found is one
   or some larger number, MM always returns a list of results, one per motif. Each result is list with three parts:
   the final values for {motifPrior, PFM, posteriors}. *)
MM[fastaFile_, motifLength_, pseudocountWeight_, maxIterations_, accuracy_, includeReverseStrand_:False, numMotifsToFind_:1, outputFile_:"MMOutput"]:=
	Module[{input, inputLengths, erasers, baseCounts, backgroundFreqs, motifPseudocounts, randomInitialPFM, initialMotifPrior, results={}},
		(*Initialization*)
		input=readInput[fastaFile, motifLength, includeReverseStrand];
		inputLengths=Map[Length, input];
		erasers=
			If[numMotifsToFind == 1,
			   False,
			   Map[ConstantArray[1.0, #] &, inputLengths]];
		baseCounts=nucleotideCounts[input];
		(*backgroundFreqs represents the frequency with which each nucleotide occurs in the
		  promoter sequences of which the input sequences are a sample, at sites not generated 
		  by the motif. These should be the vast majority sites in the promoter regions. 
		  There is a strong prior that no nucleotides should be excluded from such "generic" 
		  promoter sites. In fact, the lowest nucleotide frequency should be within roughly 
		  a factor of two of the greatest, so significant pseudocounts are justified. *) 
		backgroundFreqs=(baseCounts + {1., 1., 1., 1.}) / (Total[baseCounts] + 4.);
		motifPseudocounts=pseudocountWeight*backgroundFreqs;
		(*In the published paper they iterate through many combinations of initializations
		  for the PWM and the prior, or expected number of motif occurrences. They say this
		  is important and is the primary determinant of running time. Not implemented here.*)
		randomInitialPFM = initializePFM[motifLength];
		(* The prior probability of finding a motif at a given site is one over the average length
		   of the input sequences, so the expectataion based on this prior is one motif per input. *)
		initialMotifPrior=Length[inputLengths] / Total[inputLengths];
		Do[AppendTo[results, MMIteration[input, motifLength, motifPseudocounts, initialMotifPrior, randomInitialPFM, backgroundFreqs, 
			                             maxIterations, accuracy, erasers]];
		   PutAppend[i, results[[i]], erasers, outputFile];
		   Print["PFM number ", i, " consensus: ", PFMConsensus[results[[-1,2]]]];
		   Print[prettyPrintPFM[results[[i,2]]]];
		   (* To update the erasers after finding a motif, you need the posteriors from the last iteration. *)
		   If[numMotifsToFind > 1,
			  erasers=updateErasers[erasers, results[[-1, 3]], motifLength]],
		   {i, numMotifsToFind}];
		results
		]
		
(* MMIteration returns one motif. It makes use of the erasers to avoid finding previously reported motifs. 
   MMIteration calls updateProbs repeatedly until either convergence (the PWMs are haven't changed much
   since the last iteration) or the maximum number of iterations is reached. *)
MMIteration[input_, motifLength_, motifPseudocounts_, initialMotifPrior_, randomInitialPFM_,backgroundFreqs_, 
	        maxIterations_, accuracy_, erasers_]:=
	Module[{PFM=randomInitialPFM, oldPFM, motifPrior=initialMotifPrior, posteriors,iter=0},
		Do[oldPFM=PFM;
			iter += 1;
			
		   {motifPrior, PFM, posteriors}=
			 updateProbs[input, backgroundFreqs, PFM, motifPrior, motifPseudocounts, erasers];
		   If[Total[Map[Abs, oldPFM-PFM, {2}], 2] < accuracy,
		   	  Print[iter];Return[]],
		   {maxIterations}];
		{motifPrior, PFM, posteriors}
	]

(* Total count of each nucleotide, summing over all input sequences. *)
nucleotideCounts[input_]:=
	Total[Map[BinCounts[#, {Range[1, 5]}] &, input]]
		
(* This random initialization ensures that the initial probabilities (before normalization)
    of the different nucleotides don't differ by more than a factor of 2.*)
initializePFM[motifLength_]:=
    Map[Normalize[#, Total] &, 
       RandomReal[{0.5, 1.0},{motifLength, 4}]]

(* UpdateProbs is where one full iteration of the EM algorithm actually happens. *)
updateProbs[input_, backgroundFreqs_, oldPFM_, oldMotifPrior_, motifPseudocounts_, erasers_]:=
 	Module[{motifLength, unNormalizedPosteriors, normalizedPosteriors, newMotifPrior, newPFM},
		motifLength=Length[oldPFM];
		(* You have to write the code for sequencePosteriors. They are the posteriors calculated according to EM.
		    They are only unnormalized in the sense that the window normalization "hack" has not yet been applied. *)
		unNormalizedPosteriors = Map[sequencePosteriors[#, oldMotifPrior, oldPFM, backgroundFreqs] &, input];
		(* normalizePosteriorWindows is provided for you in fixUps.m *)
		normalizedPosteriors = normalizePosteriorWindows[unNormalizedPosteriors, motifLength];
		(* You have to write the code for updateMotifPrior. The motif prior is the probability of finding a motif
		    in a given position, prior to considering the sequence at that position. *)
		newMotifPrior=updateMotifPrior[normalizedPosteriors];
		(* You have to write the code for updatePFMCounts. This is the maximum likelihood estimate of the PFM
		    based on the expected counts of the nucleotides at each position in the PFM. *)
		newPFM = Map[Normalize[#, Total] &, 
			         updatePFMCounts[motifLength, input, normalizedPosteriors, motifPseudocounts, erasers]];
		{newMotifPrior, newPFM, normalizedPosteriors}
		]