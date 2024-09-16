(* ::Package:: *)

BeginPackage["NoisyQuantumTeleportationBenchmarking`"];


ForAllChannelsAndDistances::usage =
        "ForAllChannelsAndDistances[f] takes a function and evaluates it for all possible pairs of noisy channels and distance functions. 
		The last two arguments of f are a function to calculate the average distance of teleportation and the Pauli average of teleportation, respectively"
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

AvgDistanceOfTeleportation::"Average Distance Of Teleportation Formula"


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

Fid[r_,s_]:=(1+r . s)/2;
GetDistanceLabel[Fid] ^= "Fidelity";
GetDistanceThreshold[Fid] ^= N[2/3];
GetDistanceOutputFantasyName[Fid] ^= "fidelity";

Aff[r_,s_]:=(r . s+(1+Sqrt[1-r . r]) (1+Sqrt[1-s . s]))/((Sqrt[1+Sqrt[r . r]]+Sqrt[1-Sqrt[r . r]]) (Sqrt[1+Sqrt[s . s]]+Sqrt[1-Sqrt[s . s]]));
GetDistanceLabel[Aff] ^= "Affinity";
GetDistanceThreshold[Aff] ^= N[\[Sqrt]5/3];
GetDistanceOutputFantasyName[Aff] ^= "affinity";

WootersDist[r_,s_] := ArcCos[Sqrt[Fid[r, s]]];
GetDistanceLabel[WootersDist] ^= "Wooters Distance";
GetDistanceThreshold[WootersDist] ^= 0.589;
GetDistanceOutputFantasyName[WootersDist] ^= "wooters_distance";

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

GetLMax[chA_[pA_], chB_[pB_]] := FirstPosition[bellStates, GetMaxEigenvector[eigensystems[{chA, chB}][pA, pB]]][[1]];


RotOptFid[i_, chA_, chB_, lMax_]:= w[[i]] . w[[lMax]];
p[i_, t_, fano_]:=(1 + t . (w[[i]] . GetRA[fano])) / 4;
tBob[i_, t_, fano_, chA_, chB_, lMax_]:=1 / (4 * p[i, t, fano]) * RotOptFid[i, chA, chB, lMax] . (GetRB[fano] + Transpose[w[[i]] . GetR[fano]] . t);


Score[chA_, chB_, t_, fanoForm_, d_, lMax_] := Sum[p[i, t, fanoForm] * d[t, tBob[i, t, fanoForm, chA, chB, lMax]], {i, 1, 4}];

FanoFormForChannel[chA_[pA_], chB_[pB_]] := fanoForms[{chA, chB}][pA, pB];

AvgDistDefinition[chA_, chB_, d_] := With[{ 
			t = {Cos[\[Phi]] * Sin[\[Theta]], Sin[\[Phi]] * Sin[\[Theta]], Cos[\[Theta]]}, 
			lMax = GetLMax[chA, chB], fanoForm = FanoFormForChannel[chA, chB]}, 
			1/(4 Pi) NIntegrate[Evaluate[Score[chA, chB, t, fanoForm, d, lMax]* Sin[\[Theta]]], {\[Phi], 0, 2 Pi}, {\[Theta], 0, Pi}]]; 
			 
PaulisAvg[chA_, chB_, d_] := With[{lMax = GetLMax[chA, chB], fanoForm = FanoFormForChannel[chA, chB]}, 
	Sum[Score[chA, chB, t, fanoForm, d, lMax], {t, pauliEigenvectors}] / Length[pauliEigenvectors]];

OptimizedFidelityFormula[chA_[pA_], chB_[pB_]] := (2 * Max[eigensystems[{chA, chB}][pA, pB][[1]]] + 1)/3;


AvgDistOfTelepFormula[chA_, chB_, Fid]:= {pA, pB} |-> OptimizedFidelityFormula[chA[pA], chB[pB]];
AvgDistOfTelepFormula[DC, DC, d_]:= {pA, pB} |-> PaulisAvg[DC[pA], DC[pB], d];
AvgDistOfTelepFormula[chA_, chB_, d_]:= {pA, pB} |-> AvgDistDefinition[chA[pA], chB[pB], d];

PaulisAvgFormula[chA_, chB_, Fid]:= {pA, pB} |-> OptimizedFidelityFormula[chA[pA], chB[pB]];
PaulisAvgFormula[chA_, chB_, d_]:= {pA, pB} |-> PaulisAvg[chA[pA], chB[pB], d];


ForAllChannelsAndDistances[f_] := Do[
	f[chA, chB, d, AvgDistOfTelepFormula[chA, chB, d], PaulisAvgFormula[chA, chB, d]]
, {d, Distances}, {chA, Channels}, {chB, Channels}]


(*Map[ResourceFunction["ToPythonFunction"], eigensystems]*)


AvgDistanceOfTeleportation[chA_[pA_], chB_[pB_], d_]:= AvgDistOfTelepFormula[chA, chB, d][pA, pB];


End[];


EndPackage[];
