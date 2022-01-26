(* This file contains auxiliary functions for the HMM labs. You do not need to make any changes to it.*) 

readFasta[fastaFile_]:=
	Map[Characters, Import[fastaFile]] /. {"A"->1, "C"->2, "G"->3, "T"->4}

(* readHMM takes a HMM text file and outputs an HMM object containing the state 
 names, the alphabet the transition matrix, and emission matrix of the HMM

INPUT
 - file = path to an HMM file. An example HMM file can be found in humanMalaria.hmm.
 The format of the file is:
 - state names (one line, names in double quotes)
 - initial state probabilities (one line, one value per state)
 - state transtion matrix (as many rows as there are states, with each row
   containing the probabilities of transitioning from one state into each other state.
   Each row must sum to one.
 - alpabet (one line, names in double quotes)
 - state emission matrix (as many rows as there are states, with each row
   containing the probability of emitting one alphabet letter. Each row must sum to one).
   NOTE that after this is read in, it is transposed so that each row contains the emission
   probabilities for a given alphabet symbol and the dimensions are number of alphabet symbols
   (rows) by number of states (columns). 
   
OUTPUT
An HMM Object, which is a new symbol. To access the components, call it on
the appropriate string (see code).
*)

readHMMFile[filePath_] := 
	Module[{fileContents = Import[filePath, "CSV"],
		    numberOfStates, hmmObject},
		hmmObject["states"]=fileContents[[1]];
		numberOfStates = Length[hmmObject["states"]];
		hmmObject["initialStateProbs"] = fileContents[[2]];
		hmmObject["transitionMatrix"] = fileContents[[3 ;; 3 + numberOfStates - 1]];
		hmmObject["alphabet"] = fileContents[[3 + numberOfStates]]; 
		hmmObject["emissionMatrix"]=Transpose[fileContents[[4 + numberOfStates;; -1]]];
		hmmObject
	  ]   
 
(* checkHMMValidity returns True if its argment is a valid HMM object; otherwise it returns False.
    errorTolerance can be provided in case this is an empirically estimated HMM, where probabilities
    may not sum to exactly one.*)	  
checkHMMValidity[hmmObject_, errorTolerance_:0.00001] :=
	Module[{numberOfStates = Length[hmmObject["states"]],
		    numberOfAlphabetSymbols = Length[hmmObject["alphabet"]],
		    validHMM = True},
	If[Not[AllTrue[hmmObject["states"], StringQ]],
	   Print["Invalid HMM: The list of HMM states", hmmObject["states"], "is not a list of strings"];
	   validHMM = False];
    	If[Not[Length[hmmObject["initialStateProbs"]] == numberOfStates],
    	   Print["Invalid HMM: The list of initial state probabilities", 
    	   	     hmmObject["initialStateProbs"], 
    	   	     "is not the same length as the list of states", 
    	   	     numberOfStates];
    	   	validHMM = False];
    	If[Not[AllTrue[hmmObject["initialStateProbs"], NumberQ]], 
    	   Print["Invalid HMM: initial state probabilities", 
    	   	     hmmObject["initialStateProbs"],
    	   	     "are not all numeric"];
    	   validHMM = False];
	If[Not[checkProbsList[hmmObject["initialStateProbs"], errorTolerance]],
	       Print["Invalid HMM: initial state probabilities", 
    	   	          hmmObject["initialStateProbs"],
    	   	          "are not all numbers between zero and one that sum to 1, within the error tolerance."];
    	   validHMM = False];
    If[Dimensions[hmmObject["transitionMatrix"]] != {numberOfStates, numberOfStates},
	   Print["Invalid HMM: the transition matrix,",
	    	     hmmObject["transitionMatrix"],
	    	     "does not have dimensions of {numberOfStates, numberOfStates} ", 
      	     {numberOfStates, numberOfStates}]; 
    	   validHMM = False];
    If[Not[AllTrue[hmmObject["transitionMatrix"], checkProbsList[#, errorTolerance] &]],
	       Print["Invalid HMM: the rows of the transition matrix,",
	       	     hmmObject["transitionMatrix"],
	       	     "are not all numbers between zero and one that sum to 1, within the error tolerance."];
    	   validHMM = False]; 
    If[Not[AllTrue[hmmObject["alphabet"], StringQ]],	       	     
	   Print["Invalid HMM: The list of alphabet symbols", hmmObject["alphabet"], "is not a list of strings"];
	   validHMM = False];
    If[Dimensions[hmmObject["emissionMatrix"]] != {numberOfAlphabetSymbols, numberOfStates},
	   Print["Invalid HMM: the emission matrix,",
	       	  hmmObject["emissionMatrix"],
	       	  "does not have dimensions of {numberOfAlphabetSymbols, numberOfStates}", 
	       	  {numberOfAlphabetSymbols, numberOfStates}]; 
	       validHMM = False];
    	If[Not[AllTrue[Transpose[hmmObject["emissionMatrix"]], checkProbsList[#1, errorTolerance] &]],
	       Print["Invalid HMM: the rows of the transposed emission matrix,",
	       	     Transpose[hmmObject["emissionMatrix"]],
	       	     "are not all numbers between zero and one that sum to 1, within the error tolerance."];
 	       validHMM = False];
 	validHMM
]

	
checkProbsList[probsList_, errorTolerance_:0.00001] := 
  AllTrue[probsList, NumberQ] &&
  AllTrue[probsList, 0.0 <= # <= 1.0 &] &&	
  Abs[Total[probsList] - 1] <= errorTolerance
  
printHMM[hmmObject_]:= (
  Print["States: ", MatrixForm[hmmObject["states"]]];
  Print["Transition Matrix: ",	MatrixForm[hmmObject["transitionMatrix"]]];
  Print["Initial state probabilities: ",	 MatrixForm[hmmObject["initialStateProbs"]]];
  Print["Alphabet: ",	MatrixForm[hmmObject["alphabet"]]];
  Print["Emission Matrix (transposed relative to file): "];
  Print[MatrixForm[Prepend[Transpose[Prepend[Transpose[hmmObject["emissionMatrix"]], 
  	                                         hmmObject["alphabet"]]],
  	                       Prepend[hmmObject["states"], ""]]]];
  Print["This HMM is not guaranteed to be valid. To check it, use checkHMMValidity."]
)

calculateAccuracy[predictedStatesList_, trueStatesList_] :=
	Count[MapThread[#1 == #2 &,{predictedStatesList, trueStatesList}], True]

(*myRound is a hack to get around a problem with rounding numbers in Mathematica.Please don't bother to try to understand it.*)
myRound[x_, n_] :=
  N[IntegerPart[Round[x, 10^-n]*10^n]
                            /10^n];