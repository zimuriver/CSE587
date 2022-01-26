(* Wolfram Language package *)

(*parser for the MSA file into matrix, return a 2D array of the sequences*)
(**)
parsingMSA[msaFile_]:=
	Module[{input, seqs},
		input = StringSplit[Import["sequence.txt"]][[;; -2]];
		seqs = Map[StringPartition[#, 1] &, input];
		seqs /. {"A" -> 1, "C" -> 2, "G" -> 3, "T" -> 4}
	]

(*takes in a number of match states and the model alphabet length, and returns a skeleton profile HMM*)
skeletonPHMM[numMatch_, alphabetLength_]:=
	Module[{initialInsertEmissions, initialTransitions, emissions, transitions},
		initialInsertEmissions = ConstantArray[1, alphabetLength];
		initialTransitions = ConstantArray[1, 5];
		emissions = Table[1, numMatch, 2, alphabetLength];
		transitions = Table[1, numMatch, 7];
		transitions[[-1, 3]] = 0;
		transitions[[-1, 7]] = 0;
		{initialInsertEmissions, initialTransitions, emissions, transitions}
	]
