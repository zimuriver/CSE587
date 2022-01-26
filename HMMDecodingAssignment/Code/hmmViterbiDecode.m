(* decode
INPUT 
- observationSeq, in the form output by readFasta in tools.m
- hmm, in form output by readHMM in tools.m.
See tools.m for details.
     
OUTPUT
- stateSeq = a list containing the most likely state sequence, e.g., {h,h,m,m,m}.  
The length of this sequence equal to the length to the input observationSeq.  

IMPLEMENTATION
I recommend the following:
1) Implement the Viterbi table as a matrix that is transposed relative to the way
   it is shown in class and course notes. In other words, the rows correspond to
   observations and the columns to states.
2) After computing the Viterbi probabilities for each observation, normalize
   them by dividing by their sum. The matrix entries are no longer the Viterbi 
   probabilities, but the entries for each observation remain proportional
   to the Viterbi probabilities. The only thing you will use this matrix for is picking
   the state with the greatest entry for each observation, so scaling them all by a
   constant factor won't change the result. The benefit of doing this is that you
   avoid numerical underflow.
*)

viterbiDecode[observationSeq_, hmm_] :=
  Module[{numericStateSeq},
  (* Put your code here. *)
  numericStateSeq=traceback[buildMatrix[observationSeq,hmm],hmm];
  Map[hmm["states"][[#]]&,numericStateSeq]
  (* Return the sequence of state names corresponding to the Viterbi decode path. *)

 ]

buildMatrix[observationSeq_, hmm_] := 
	Module[{numberOfObservations, numberOfStates, viterbiMatrix, observationIndex},
	  (* Put your code for bulding the Viterbi matrix here. To give you an idea of what to expect,
         my code is 13 lines. But use as many lines as you need to make your code clear and readable. *)
         numberOfObservations=Length[observationSeq];
         numberOfStates=  Length[hmm["states"]];
         viterbiMatrix={};  
         (* calculating the Vertibi probability for begin *)
         AppendTo[viterbiMatrix,hmm["emissionMatrix"][[observationSeq[[1]]]]*hmm["initialStateProbs"]];
         (* normalize for each list in vertibi matrix to prevent getting too small after Times *)
         If[Total[viterbiMatrix[[-1]]]==0,viterbiMatrix[[-1]]=viterbiMatrix[[-1]]/0
         	,viterbiMatrix[[-1]]=Normalize[viterbiMatrix[[-1]],Total]];
         (* calculating the Vertibi probability for each state at that time after begin,
         find bi(Ot) for current state i, find the max for the product of aj,i the transition probs from different states to current state=i and previous vertibi probs *)
         Do[observationIndex=observationSeq[[iterationNumberForObservation]];
         AppendTo[viterbiMatrix,
         hmm["emissionMatrix"][[observationIndex]]*
         Map[Max[Times[#,viterbiMatrix[[iterationNumberForObservation-1]]]]&,Transpose[hmm["transitionMatrix"]]]
         ];If[Total[viterbiMatrix[[-1]]]==0,viterbiMatrix[[-1]]=viterbiMatrix[[-1]]/0,viterbiMatrix[[-1]]=Normalize[viterbiMatrix[[-1]],Total]]
         ,{iterationNumberForObservation,2,numberOfObservations}];
         
         (* Print[viterbiMatrix];*)
	 (* Return the Viterbi matrix. *)
	 viterbiMatrix] 

traceback[viterbiMatrix_, hmm_] :=
  (* Put your code for tracing back through the Viterbi matrix here. To give you an idea of what to expect,
     my code is 10 lines. But use as many lines as you need to make your code clear and readable. *)
     Module[{finalStateIndex,lengthOfStateSeq,stateSeq,sourcesForMax,maxCalled,maxIndex},
     	(* find the final state of the most likely state seq *)
     	lengthOfStateSeq=Length[viterbiMatrix];
     	maxCalled=Max[viterbiMatrix[[-1]]];
     	finalStateIndex=Min[Flatten[Position[viterbiMatrix[[-1]],maxCalled]]];
     	stateSeq={finalStateIndex};
     	(* To find the mostlikely previous state for the most likely state seq,
     	find the vertibi probs for previous state, 
     	find the transition probs for previous different states to current state, 
     	find the max source from previous state,break tie by chose smaller index *)
     	Do[
     	 sourcesForMax=Times[Transpose[hmm["transitionMatrix"]][[stateSeq[[1]]]],viterbiMatrix[[iterationNumberForStateSeq]]];
     	 (* Print[sourcesForMax]; *)
     	 maxCalled=Max[sourcesForMax];
     	 maxIndex=Min[Flatten[Position[sourcesForMax,maxCalled]]];
     	 PrependTo[stateSeq,maxIndex]
     	,{iterationNumberForStateSeq,Reverse[Range[lengthOfStateSeq]][[2;;-1]]}];
	 (* Return the list of states. *)
	 stateSeq
	]	   
 
	