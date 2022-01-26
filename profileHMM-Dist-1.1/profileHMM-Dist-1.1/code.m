<<tools`
 
createPHMM[msaFile_]:=
	Module[{sequences, matches, matchStates, numMatch, pHMM, i},
		(*read in the msaFile*)
		sequences = parsingMSA[msaFile];
		(*determine which positions in the sequences are match states*)
		matches = Select[Map[Select[#, StringQ] &, Transpose[sequences]], Length[#] < Length[sequences]/2 &];
		matchStates = Sort[DeleteDuplicates[Flatten[Map[Position[Map[Select[#, StringQ] &, Transpose[sequences]], #] &, matches]]]];
		numMatch = Length[matchStates];
		(*get skeleton PHMM, assuming alphabet is DNA (4 letters)*)
		pHMM=skeletonPHMM[numMatch, 4];
		(*loop through all sequences, updating the PHMM counts*)
		For[i=1, i<=Length[sequences], i++,
			pHMM=updatePHMM[sequences[[i]], matchStates, pHMM];
			Print[sequences[[i]]];
		];
		(*normalize the PHMM counts to return PHMM probabilities*)
		normalizePHMM[pHMM]
	]

(*
Reminder about the PHMM structure:
Insert0 state emission probabilities (list length = length of alphabet)
Transition probabilities for B -> M1, B -> I0, B -> D1; I0 -> M1, I0 -> I0
A 3 D emission array, 
	where the first level corresponds to the index of the states (so if we have 10 match states, we expect length of 10 for the emission array), 
	the second level corresponds to the state type (Insert, or Match, so we expect length of 2 for the emission subarrays), 
	and the last level corresponds to the alphabet emission probabilities (so for the DNA alphabet, we expect length of 4 for the emission sub - subarray)
A 2 D transition array, 
	where the first level corresponds to the index of the the states (so if we have 10 match states, we expect length of 10 for the transition array), 
	and the second level corresponds to the transition probabilities to the next states 
		(so we expect 7 values for each subarray in the transition array corresponding to 
		Mk -> Mk + 1, Mk -> Ik, Mk -> Dk + 1; Ik -> Mk + 1, Ik -> Ik; Dk -> Mk + 1, Dk -> Dk + 1)

update PHMM pseudocounts based on given sequence
*)
updatePHMM[sequence_, matchStates_, pHMM_]:=
	Module[{pHMMIndex, currentState, sequenceIndex,
		initialInsertEmissions, initialTransitions, emissions, transitions,firstMatchIndex = matchStates[[1]],i,prevState},
		
		{initialInsertEmissions, initialTransitions, emissions, transitions}=pHMM;
	
		Do[
			If[MemberQ[matchStates,i];,
				initialInsertEmissions[sequence[[i]]]++;
				Break[]]
			,{i,Length[sequence]}];
			
	Do[		
		If[checkDeletion[sequence,matchStates,i],initialTransitions[[3]]++;];
		
		If[checkInsertion[sequence,matchStates,i], initialTransitions[[2]]++;
			If[checkInsertion[sequence,matchStates,i+1], initialTransitions[[5]]++;,initialTransitions[[4]]++;];
		];
		If[checkMatch[sequence,matchStates,i], initialTransitions[[1]]++;];
	,{i,firstMatchIndex}]; 	

		(*
		Walk through given sequence.
		Note that transitions up to the first match state will update the initial transitions list of the PHMM
		All following transitions will update the 2D transitions array
		*)
		currentState = matchStates[[1]]; 
		prevState = 1;
		Do[
		(*check each state and update the emissions matrix*)	
			If[checkMatch[sequence,matchStates,i], emissions[[currentState,2,sequence[[i]]]]++;];
			If[checkInsertion[sequence,matchStates,i],emissions[[currentState,1,sequence[[i]]]]++;];
			
			(*TRANSITIONPART*)
			If[i<Length[sequence],
				(*check for gaps and adjust state check accordingly*)
				If[sequence[[prevState+1]] != "-" || MemberQ[matchStates,prevState+1] || MemberQ[matchStates,i], prevState = i];
				(*check each transition*)  (*MemberQ is probably pretty inefficient here but it works fine for these smaller sequences and i couldnt think of another way*)
				If[checkMatch[sequence,matchStates,prevState] && checkMatch[sequence,matchStates,i+1], transitions[[currentState,1]]++];
				If[checkMatch[sequence,matchStates,prevState] && checkInsertion[sequence,matchStates,i+1], transitions[[currentState,2]]++];
				If[checkMatch[sequence,matchStates,prevState] && checkDeletion[sequence,matchStates,i+1], transitions[[currentState,3]]++;];
				If[checkInsertion[sequence,matchStates,prevState] && checkMatch[sequence,matchStates,i+1], transitions[[currentState,4]]++];
				If[checkInsertion[sequence,matchStates,prevState] && checkInsertion[sequence,matchStates,i+1], transitions[[currentState,5]]++];
				If[checkDeletion[sequence,matchStates,prevState] && checkMatch[sequence,matchStates,i+1], transitions[[currentState,6]]++];
				If[checkDeletion[sequence,matchStates,prevState] && checkDeletion[sequence,matchStates,i+1], transitions[[currentState,7]]++];
			
			];
			If[i<Length[sequence],If[MemberQ[matchStates,i+1],currentState++;]];
		,{i,1,Length[sequence]}];
		
		(*Update the transitions into the end state (Mk+1)*)
		If[checkMatch[sequence,matchStates,Length[sequence]],transitions[[Length[transitions],1]]++;]; 
		
		(*return updated pseudocounts*)
		{initialInsertEmissions, initialTransitions, emissions, transitions}
	]
	
	
(*Normalize the parts of the PHMM to get probability values instead of counts
Keep in mind that probabilities should not add up to one for every sublist, but rather for:
	transitions out of the same state
	emissions from the same state
*)
normalizePHMM[{initialInsertEmissions_, initialTransitions_, emissions_, transitions_}]:=
	Module[{normalizedInitialInsertEmissions, normalizedInitialTransitions, normalizedEmissions, normalizedTransitions,i,j},
	
	normalizedInitialInsertEmissions = initialInsertEmissions/Total[initialInsertEmissions];
	normalizedInitialTransitions = initialTransitions;
	normalizedInitialTransitions[[1;;3]] = initialTransitions[[1;;3]]/Total[initialTransitions[[1;;3]]];
	normalizedInitialTransitions[[4;;5]] = initialTransitions[[4;;5]]/Total[initialTransitions[[4;;5]]];
	normalizedEmissions = emissions;
	Do[
		Do[
			normalizedEmissions[[i,j]] /= Total[normalizedEmissions[[i,j]]];
			,{j,Length[emissions[[i]]]}];,{i,Length[emissions]}];
	normalizedTransitions = transitions;	
	Do[	
	normalizedTransitions[[i,1;;3]] = transitions[[i,1;;3]]/Total[transitions[[i,1;;3]]];
	normalizedTransitions[[i,4;;5]] = transitions[[i,4;;5]]/Total[transitions[[i,4;;5]]];
	normalizedTransitions[[i,6;;7]] = transitions[[i,6;;7]]/Total[transitions[[i,6;;7]]];
			,{i,Length[transitions]}];
		(*Return the normalized PHMM*)
		{normalizedInitialInsertEmissions, normalizedInitialTransitions, normalizedEmissions, normalizedTransitions}
	]
	
	(*check for different states using their conditions*)
	checkInsertion[sequence_,matchStates_,i_] := sequence[[i]] != "-" && !MemberQ[matchStates,i];
	checkDeletion[sequence_,matchStates_,i_] := sequence[[i]] == "-" && MemberQ[matchStates,i];
	checkMatch[sequence_,matchStates_,i_] := sequence[[i]] != "-" && MemberQ[matchStates,i];
	