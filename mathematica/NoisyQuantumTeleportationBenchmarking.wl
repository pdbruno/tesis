(* ::Package:: *)

BeginPackage["NoisyQuantumTeleportationBenchmarking`"];


ForAllChannelsAndDistances::usage =
        "ForAllChannelsAndDistances[f] takes a function and evaluates it for all possible pairs of noisy channels and distance functions. 
		The last argument of f is a function to calculate the average distance of teleportation"
GetChannelLabel::usage = "GetChannelLabel[ch] returns a string representing the channel's description"
GetChannelOutputFantasyName::usage = "GetChannelOutputFantasyName[ch] returns a string representing the channel's output id"

GetDistanceLabel::usage = "GetDistanceLabel[d] returns a string representing the distances's description"
GetDistanceThreshold::usage = "GetDistanceThreshold[d] returns the distances's classical threshold"
GetDistanceOutputFantasyName::usage = "GetDistanceOutputFantasyName[d] returns a string representing the distances's output id"

ADC::usage = "Symbol that represents the Amplitude Damping Channel"
MADC::usage = "Symbol that represents the Mirrored Amplitude Damping Channel"
DC::usage = "Symbol that represents the Depolarizing Channel"
PDC::usage = "Symbol that represents the Phase Damping Channel"
Channels::usage = "List with all the noisy channels"

TrDist::usage = "Trace Distance function for qubits in Bloch Representation"
Fid::usage = "Fidelity function for qubits in Bloch Representation"
Aff::usage = "Affinity function for qubits in Bloch Representation"
WootersDist::usage = "Wooters Distance function for qubits in Bloch Representation"
Distances::usage = "List with all the distance functions"

AvgDistanceOfTeleportation::usage = "Average Distance Of Teleportation Formula"
QuantumCertified::usage = "Decides whether a value for average distance of teleportation can be quantum certified for a given distance measure"

ResourceNegativityWithChannels::usage = "Gives the negativity of a Bell pair exposed to chA and chB"
ChannelsMakeResourceSeparable::usage = "Returns whether a Bell pair exposed to chA and chB is separable"


Begin["`Private`"]


(*Funciones generales*)
Fano[\[Rho]_]:={{ 2Re[\[Rho][[2,3]]+\[Rho][[1,4]]], 2Im[\[Rho][[2,3]]\[Minus]\[Rho][[1,4]]], 2Re[\[Rho][[1,3]]\[Minus]\[Rho][[2,4]]]},
{\[Minus]2Im[\[Rho][[2,3]]+\[Rho][[1,4]]], 2Re[\[Rho][[2,3]]\[Minus]\[Rho][[1,4]]],\[Minus]2Im[\[Rho][[1,3]]\[Minus]\[Rho][[2,4]]]},{ 2Re[\[Rho][[1,2]]\[Minus]\[Rho][[3,4]]],\[Minus]2Im[\[Rho][[1,2]]\[Minus]\[Rho][[3,4]]], \[Rho][[1,1]]\[Minus]\[Rho][[2,2]]\[Minus]\[Rho][[3,3]]+\[Rho][[4,4]]}};
DensityMatrix[s_] := Outer[Times, s, s];
ExpectationValue[psi_, m_] := ConjugateTranspose[psi] . m . psi;
ToBloch[psi_] := {ExpectationValue[psi, PauliMatrix[1]], ExpectationValue[psi, PauliMatrix[2]], ExpectationValue[psi, PauliMatrix[3]]};


(*Definici\[OAcute]n de los estados de Bell*)
ket[0] = {1, 0};
ket[1] = {0, 1};
bell[x_, y_] := Flatten[(KroneckerProduct[ket[0], ket[y]] + (-1)^x KroneckerProduct[ket[1], ket[1-y]])/Sqrt[2]];
bellStates = Flatten[Table[bell[x, y], {x, 0, 1}, {y, 0, 1}], 1];


(*C\[AAcute]lculo de matrices w*)
w = Table[Fano[DensityMatrix[bell]], {bell, bellStates}];


(*C\[AAcute]lculo de autovectores de matrices de Pauli*)
pauliEigenvectors = (ToBloch[Normalize[#]])& /@ Flatten[Table[Eigenvectors[PauliMatrix[i]], {i, 3}], 1];


(*Definici\[OAcute]n de las funciones de distancia*)
TrDist[r_,s_]:=Sqrt[(r-s) . (r-s)]/2;
GetDistanceLabel[TrDist] ^= "Trace Distance";
GetDistanceThreshold[TrDist] ^= N[8(11\[Minus]2\[Sqrt]10)/81];
GetDistanceOutputFantasyName[TrDist] ^= "trace_distance";
QuantumCertified[TrDist, avgDistTelep_] ^:= avgDistTelep < GetDistanceThreshold[TrDist];

Fid[r_,s_]:=(1+r . s)/2;
GetDistanceLabel[Fid] ^= "Fidelity";
GetDistanceThreshold[Fid] ^= N[2/3];
GetDistanceOutputFantasyName[Fid] ^= "fidelity";
QuantumCertified[Fid, avgDistTelep_] ^:= avgDistTelep > GetDistanceThreshold[Fid];

Aff[r_,s_]:=(r . s+(1+Sqrt[1-r . r]) (1+Sqrt[1-s . s]))/((Sqrt[1+Sqrt[r . r]]+Sqrt[1-Sqrt[r . r]]) (Sqrt[1+Sqrt[s . s]]+Sqrt[1-Sqrt[s . s]]));
GetDistanceLabel[Aff] ^= "Affinity";
GetDistanceThreshold[Aff] ^= N[\[Sqrt]5/3];
GetDistanceOutputFantasyName[Aff] ^= "affinity";
QuantumCertified[Aff, avgDistTelep_] ^:= avgDistTelep > GetDistanceThreshold[Aff];

WootersDist[r_,s_] := ArcCos[Sqrt[Fid[r, s]]];
GetDistanceLabel[WootersDist] ^= "Wooters Distance";
GetDistanceThreshold[WootersDist] ^= 0.589;
GetDistanceOutputFantasyName[WootersDist] ^= "wooters_distance";
QuantumCertified[WootersDist, avgDistTelep_] ^:= avgDistTelep < GetDistanceThreshold[WootersDist];

Distances = { TrDist, Fid, Aff, WootersDist };


Channels = { ADC, MADC, DC, PDC };


GetAffineDecompositionM[ADC[p_]] ^:= DiagonalMatrix[{Sqrt[1-p],Sqrt[1-p],1-p}];
GetAffineDecompositionC[ADC[p_]] ^:= {0,0,p};
GetChannelLabel[ADC] ^= "Amplitude Damping Channel";
GetChannelOutputFantasyName[ADC] ^= "ADC";

GetAffineDecompositionM[MADC[p_]] ^:= DiagonalMatrix[{Sqrt[1-p],Sqrt[1-p],1-p}];
GetAffineDecompositionC[MADC[p_]] ^:= {0,0,-p};
GetChannelLabel[MADC] ^= "Mirrored Amplitude Damping Channel";
GetChannelOutputFantasyName[MADC] ^= "MADC";

GetAffineDecompositionM[DC[p_]] ^:= DiagonalMatrix[{1-p,1-p,1-p}];
GetAffineDecompositionC[DC[p_]] ^:= {0,0,0};
GetChannelLabel[DC] ^= "Depolarizing Channel";
GetChannelOutputFantasyName[DC] ^= "DC";

GetAffineDecompositionM[PDC[p_]] ^:= DiagonalMatrix[{Sqrt[1-p],Sqrt[1-p],1}];
GetAffineDecompositionC[PDC[p_]] ^:= {0,0,0};
GetChannelLabel[PDC] ^= "Phase Damping Channel";
GetChannelOutputFantasyName[PDC] ^= "PDC";


GetCorrelationMatrixForNoises[chA_, chB_] := GetAffineDecompositionM[chA] . w[[1]] . GetAffineDecompositionM[chB] + Outer[Times, GetAffineDecompositionC[chA], GetAffineDecompositionC[chB]];

FanoForm[chA_, chB_] := {GetAffineDecompositionC[chA], GetAffineDecompositionC[chB], GetCorrelationMatrixForNoises[chA, chB]};
GetRA[{rA_, _, _}] := rA;
GetRB[{_, rB_, _}] := rB;
GetR[{_, _, r_}] := r;


BindFreeVars = Function[{pA, pB}, #]&;


AssociateChannelsPairsTo[f_] := Association[Table[{chA, chB} -> f[chA, chB], {chA, Channels},{chB, Channels}]];


fanoForms = AssociateChannelsPairsTo[{chA, chB} |->  BindFreeVars[FanoForm[chA[pA], chB[pB]]]];


GetBellDiagonalResourceStateForNoises[chA_[pA_], chB_[pB_]] := 
	(Join[{1}, Diagonal[GetR[fanoForms[{chA, chB}][pA, pB]]]] . Table[KroneckerProduct[PauliMatrix[i],PauliMatrix[i]],{i,0,3}])/4;

normalizeEigenvectors[{eigenvalues_, eigenvectors_}] := {eigenvalues, eigenvectors * 1/Sqrt[2]};

eigensystems = AssociateChannelsPairsTo[{chA, chB} |-> BindFreeVars[normalizeEigenvectors[FullSimplify[Eigensystem[GetBellDiagonalResourceStateForNoises[chA[pA], chB[pB]]], 0<=pA <=1 && 0<=pB <=1]]]];


(* If one of the ev is already a bell state, return it. 
	Else, try to do a linear combination (adding them all, maybe a little rough). 
	If the linear combination isn't a bell state, supose any choice is valid, so return the first bell state *)
	
(* Possibly refactor using better suiting WL functions to improve preformance/readability (ojo con el Which) *)
FindMatchingBellState[maxEigenvectors_] := SelectFirst[maxEigenvectors, MemberQ[bellStates, #]&, If[MemberQ[bellStates, Total[maxEigenvectors]], Total[maxEigenvectors], bellStates[[1]]]];

GetEigenvectorsOfMaxEigenvalue[eigenvalues_, eigenvectors_] := eigenvectors[[Flatten[Position[eigenvalues, Max[eigenvalues]]]]];

GetMaxEigenvector[{eigenvalues_, eigenvectors_}] := 
	FindMatchingBellState[GetEigenvectorsOfMaxEigenvalue[eigenvalues, eigenvectors]];

GetLMax[eigensystem_] := FirstPosition[bellStates, GetMaxEigenvector[eigensystem]][[1]];


RotOptFid[i_, lMax_]:= w[[i]] . w[[lMax]];
p[i_, t_, fano_]:=(1 + t . (w[[i]] . GetRA[fano])) / 4;
tBob[i_, t_, fano_, lMax_]:=1 / (4 * p[i, t, fano]) * RotOptFid[i, lMax] . (GetRB[fano] + Transpose[w[[i]] . GetR[fano]] . t);


Score[t_, fanoForm_, d_, lMax_] := Sum[p[i, t, fanoForm] * d[t, tBob[i, t, fanoForm, lMax]], {i, 1, 4}];

FanoFormForChannel[chA_[pA_], chB_[pB_]] := fanoForms[{chA, chB}][pA, pB];

AvgDistDefinition[d_, fanoForm_, lMax_] := 1/(4 Pi) NIntegrate[Evaluate[Score[FromSphericalCoordinates[{1, \[Theta], \[Phi]}], fanoForm, d, lMax]* Sin[\[Theta]]], {\[Phi], 0, 2 Pi}, {\[Theta], 0, Pi}]; 
			 
PaulisAvg[d_, fanoForm_, lMax_] := Sum[Score[t, fanoForm, d, lMax], {t, pauliEigenvectors}] / Length[pauliEigenvectors];


AvgDistOfTelepFormula[chA_, chB_, Fid]:= With[{eigensystem = eigensystems[{chA, chB}]}, {pA, pB} |-> (2 * Max[eigensystem[pA, pB][[1]]] + 1)/3];
AvgDistOfTelepFormula[DC, DC, d_]:= With[{fanoForm = fanoForms[{DC, DC}], eigensystem = eigensystems[{DC, DC}]}, 
										{pA, pB} |-> PaulisAvg[d, fanoForm[pA, pB], GetLMax[eigensystem[pA, pB]]]];
										
AvgDistOfTelepFormula[chA_, chB_, d_]:= With[{fanoForm = fanoForms[{chA, chB}], eigensystem = eigensystems[{chA, chB}]}, 
										{pA, pB} |-> AvgDistDefinition[d, fanoForm[pA, pB], GetLMax[eigensystem[pA, pB]]]];


ForAllChannelsAndDistances[f_] := Do[
	f[chA, chB, d, AvgDistOfTelepFormula[chA, chB, d]]
, {chA, Channels}, {chB, Channels},{d, Distances}]


(*Map[ResourceFunction["ToPythonFunction"], eigensystems]*)


AvgDistanceOfTeleportation[chA_[pA_], chB_[pB_], d_]:= AvgDistOfTelepFormula[chA, chB, d][pA, pB];


Manipulate[ListPointPlot3D[{SpherePoints[100], (point |-> GetAffineDecompositionM[ch[p]] . point + GetAffineDecompositionC[ch[p]]) /@ SpherePoints[100]},BoxRatios->1 ],{p,0,1}, {ch, {ADC, MADC, DC, PDC}}]


MADCGate[\[Theta]_] = {{0, 0, 0, 1}, {Cos[\[Theta]/2], Sin[\[Theta]/2], 0, 0}, {-Sin[\[Theta]/2], Cos[\[Theta]/2], 0, 0}, {0, 0, 1, 0}};
BlochToDensityMatrix[t_] := 1/2*(IdentityMatrix[2] + Sum[t[[i]]*PauliMatrix[i], {i, 3}]);
DensityMatrixToBloch[\[Rho]_] := Table[Tr[\[Rho] . KroneckerProduct[PauliMatrix[i], IdentityMatrix[2]]], {i, 3}];
channel[p_] := With[{gate = MADCGate[2*ArcSin[Sqrt[p]]]},  t |-> DensityMatrixToBloch[gate . KroneckerProduct[BlochToDensityMatrix[t], DensityMatrix[ket[0]]] . ConjugateTranspose[gate]]];
Manipulate[ListPointPlot3D[{SpherePoints[100], channel[p] /@ SpherePoints[100]},BoxRatios->1 ],{p,0,1}]


Pt\[Rho][A1_,b1_,A2_,b2_,ww_]:=1/4 (IdentityMatrix[4]+KroneckerProduct[Sum[b1[[i]]*PauliMatrix[i],{i,1,3}],IdentityMatrix[2]]+KroneckerProduct[IdentityMatrix[2],Sum[b2[[i]]*Transpose[PauliMatrix[i]],{i,1,3}]]+Sum[(b1[[\[Alpha]]]*b2[[\[Beta]]]+Sum[ww[[i,i]]*A1[[\[Alpha],i]]*A2[[\[Beta],i]],{i,1,3}])*KroneckerProduct[PauliMatrix[\[Alpha]],Transpose[PauliMatrix[\[Beta]]]],{\[Alpha],1,3},{\[Beta],1,3}]);


(*AffectedBellPairDensityMatrixFromFano[{rA_, rB_, r_}] := KroneckerProduct[BlochToDensityMatrix[rA], BlochToDensityMatrix[rB]] + Sum[r[[i, j]]*KroneckerProduct[PauliMatrix[i], PauliMatrix[j]], {i, 1, 3}, {j, 1, 3}];
AffectedBellPairDensityMatrixFromFano[{rA_, rB_, r_}] := KroneckerProduct[BlochToDensityMatrix[rA], Transpose[BlochToDensityMatrix[rB]]] + Sum[r[[i, j]]*KroneckerProduct[PauliMatrix[i], Transpose[PauliMatrix[j]]], {i, 1, 3}, {j, 1, 3}];
AffectedBellPairDensityMatrixFromChannels[chA_[pA_], chB_[pB_]] := AffectedBellPairDensityMatrixFromFano[FanoFormForChannel[chA[pA], chB[pB]]]*)


PartialTranspose = ArrayFlatten @ Map[Transpose, Partition[#, {2, 2}], {2}] &;


TraceNorm[X_] := Tr[Sqrt[ConjugateTranspose[X] . X]]


ResourceNegativityWithChannels[chA_, chB_]:= Module[{eig},
eig=Eigenvalues[Pt\[Rho][GetAffineDecompositionM[chA], GetAffineDecompositionC[chA], GetAffineDecompositionM[chB],GetAffineDecompositionC[chB], w[[1]]]];
Abs[Total[Select[eig,#<0&]]]
]
ChannelsMakeResourceSeparable[chA_, chB_] := ResourceNegativityWithChannels[chA, chB] == 0


teleport[t_, chA_[pA_], chB_[pB_]] := Table[tBob[i, t, FanoFormForChannel[chA[pA], chB[pB]], GetLMax[eigensystems[{chA, chB}][pA, pB]]], {i, 1, 4}]
Manipulate[
 Graphics3D[{{Opacity[0.5], Sphere[]}, 
   {EdgeForm[{Dashed, Red}], 
    FaceForm[None], Cylinder[{{0, 0, -.001}, {0, 0, .001}}]}, 
   {Dashed, Line[{{0, 0, 0}, {1, 0, 0}}]}, 
   {Dashed, Line[{{0, 0, 0}, {0, 1, 0}}]}, 
   {Dashed, Line[{{0, 0, 0}, {0, 0, 1}}]}, 
   Arrow[{{1, 0, 0}, {1.3, 0, 0}}], 
   Arrow[{{0, 1, 0}, {0, 1.3, 0}}], 
   Arrow[{{0, 0, 1}, {0, 0, 1.3}}], 
   Text[Style["x", 18], {1.4, 0, 0}], 
   Text[Style["y", 18], {0, 1.4, 0}], 
   Text[Style["z", 18], {0, 0, 1.4}], 
   Text[Style["t", 18, Red], .5 {Cos[t[[1]]] Sin[t[[2]]], Sin[t[[1]]] Sin[t[[2]]], Cos[t[[2]]]}], 
   {Red, Arrow[{{0, 0, 0}, {Cos[t[[1]]] Sin[t[[2]]], Sin[t[[1]]] Sin[t[[2]]], Cos[t[[2]]]}}]},
  Sequence[({Green, Arrow[{{0, 0, 0}, #}]})& /@ teleport[{Cos[t[[1]]] Sin[t[[2]]], Sin[t[[1]]] Sin[t[[2]]], Cos[t[[2]]]}, chA[pA], chB[pB]]]
  }, 
  Boxed -> False], 
  {t,{0, 0},{2*\[Pi], \[Pi]}},
  {chA, {ADC, MADC, DC, PDC}},
  {pA, 0, 1},
  {chB, {ADC, MADC, DC, PDC}},
  {pB, 0, 1}
  ]


End[];


EndPackage[];
