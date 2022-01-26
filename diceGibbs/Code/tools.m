(* Wolfram Language package *)

stackStats[labelsStack_, thetasStack_]:=
	(Print["Label means"];
	 Print[Round[Mean[labelsStack] - 1, 0.01]];
	 Print["thetaStack"];
	 Do[
	 	Print["die Type ", dieType];
	 Print[Grid[Prepend[Map[Round[{Quantile[#, 0.05], Quantile[#, 0.5], Quantile[#, 0.95]}, 0.01] &, Transpose[thetasStack[[All, dieType]]]],
	 	               {"5%", "Median", "95%"}]]],
	 {dieType, 1, 2}	]
	 )

(* This is the right thing to do when calculating the likelihood of a sample from a multinomial distribution
   and the exponent represents the number of times an outcome was seen. If an outcome was never seen, it doesn't
   affect the likelihood, even its probability is zero. *)
safeExponentiate[base_, exponent_]:=
	If[TrueQ[N[exponent] == 0.0],
	   1,
	   base^exponent]

diceSample[type2Prob_, type1FaceProbs_, type2FaceProbs_, draws_, rollsPerDraw_] :=
 Module[{drawTypes= 1 + RandomVariate[BernoulliDistribution[type2Prob],
 	                                  draws]},
 	Map[RandomChoice[{type1FaceProbs, type2FaceProbs}[[#]] -> Range[Length[type1FaceProbs]],
 		              rollsPerDraw] &, drawTypes]
 	]