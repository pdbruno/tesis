(* ::Package:: *)

BeginPackage["NoisyQuantumTeleportationBenchmarking`"];


ForAllChannelsAndDistances::usage =
        "ForAllChannelsAndDistances[f] takes a function and evaluates it for all possible pairs of noisy channels and distance functions. 
		The last two arguments of f are a function to calculate the average distance of teleportation and the Pauli average of teleportation, respectively"
GetChannelLabel::usage = "GetChannelLabel[ch] returns a string representing the channel's description"

GetAffineDecompositionC::usage = "asd"
GetAffineDecompositionM::usage = "asd"
Channels::usage = "asd"


Begin["`Private`"]


(*Definici\[OAcute]n de matrices w*)
w={DiagonalMatrix[{1,-1,1}],DiagonalMatrix[{1,1,-1}], DiagonalMatrix[{-1,1,1}],DiagonalMatrix[{-1,-1,-1}]};
pauliEigenvectors = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};


(*Definici\[OAcute]n de las funciones de distancia*)
TrDist[r_,s_]:=Sqrt[(r-s) . (r-s)]/2;
Fid[r_,s_]:=(1+r . s)/2;
Aff[r_,s_]:=(r . s+(1+Sqrt[1-r . r]) (1+Sqrt[1-s . s]))/((Sqrt[1+Sqrt[r . r]]+Sqrt[1-Sqrt[r . r]]) (Sqrt[1+Sqrt[s . s]]+Sqrt[1-Sqrt[s . s]]));
WootersDist[r_,s_] := ArcCos[Sqrt[Fid[r, s]]];

Distances = {
{TrDist, "trace_distance",N[8(11\[Minus]2\[Sqrt]10)/81]}, 
{Fid, "fidelity",N[2/3]}, 
{Aff, "affinity",N[\[Sqrt]5/3] }, 
{WootersDist, "wooters_distance",0.589}
};


(*Clear[ChannelParameterQ];
ChannelParameterQ[p_] := RealQ[p] && 0<=p && p<=1;

ADC[p_] /; !ChannelParameterQ[p] := Throw[$Failed, ADC];
MADC[p_] /; !ChannelParameterQ[p] := Throw[$Failed, MADC];
DC[p_] /; !ChannelParameterQ[p] := Throw[$Failed, DC];
PDC[p_] /; !ChannelParameterQ[p] := Throw[$Failed, PDC];*)

GetAffineDecompositionM[ADC[p_]] := DiagonalMatrix[{Sqrt[1-p],Sqrt[1-p],1-p}];
GetAffineDecompositionC[ADC[p_]] := {0,0,p};
GetChannelLabel[ADC[_]] := "Amplitude Damping Channel";

GetAffineDecompositionM[MADC[p_]] := DiagonalMatrix[{Sqrt[1-p],Sqrt[1-p],1-p}];
GetAffineDecompositionC[MADC[p_]] := {0,0,-p};
GetChannelLabel[MADC[_]] := "Mirrored Amplitude Damping Channel";

GetAffineDecompositionM[DC[p_]] := DiagonalMatrix[{1-p,1-p,1-p}];
GetAffineDecompositionC[DC[p_]] := {0,0,0};
GetChannelLabel[DC[_]] := "Depolarizing Channel";

GetAffineDecompositionM[PDC[p_]] := DiagonalMatrix[{Sqrt[1-p],Sqrt[1-p],1}];
GetAffineDecompositionC[PDC[p_]] := {0,0,0};
GetChannelLabel[PDC[_]] := "Phase Damping Channel";

GetLMax[ADC[pA_], MADC[pB_]] := If[(1+Sqrt[1+2 pB-3 pB^2]<2 pA+pB), 3, 1];
GetLMax[MADC[pA_], ADC[pB_]] := GetLMax[ADC[pB], MADC[pA]];
GetLMax[_, _] := 1

Channels = { ADC, MADC, DC, PDC };


ForAllChannelsAndDistances[f_] := Do[
	
	rA = GetAffineDecompositionC[chA];
	rB = GetAffineDecompositionC[chB];
	r = GetAffineDecompositionM[chA] . w[[1]] . GetAffineDecompositionM[chB] + Outer[Times, rA, rB];
	lmax = GetLMax[chA, chB];
	RotOptFid[i_]:= w[[i]] . w[[lmax]];
	p[i_, t_]:=(1+t . (w[[i]] . rA))/4;
	tBob[i_, t_]:=1/(4*p[i, t])*RotOptFid[i] . (rB+Transpose[w[[i]] . r] . t);

	Do[
		(*Definici\[OAcute]n de AvgDist*)
		intVar = {Cos[\[Phi]]*Sin[\[Theta]],Sin[\[Phi]]*Sin[\[Theta]],Cos[\[Theta]]};
		AvgDist[pA_,pB_]:= 1/(4 Pi) NIntegrate[Evaluate[Sum[p[i, intVar]*d[[1]][intVar,tBob[i, intVar]],{i,1,4}]*Sin[\[Theta]]],{\[Phi],0,2 Pi},{\[Theta],0,Pi}];
		PaulisAvg[pA_,pB_]:= Sum[Sum[p[i, t]*d[[1]][t,tBob[i, t]],{i,1,4}], {t, pauliEigenvectors}]/6;
		
		f[chA, chB, d, AvgDist, PaulisAvg];
		
	,{d, Distances}]
	
,{chA, (#[pA])&/@ Channels},{chB, (#[pB])&/@ Channels}]


ForAllChannels[f_] := Do[
	
	rA = GetAffineDecompositionC[chA];
	rB = GetAffineDecompositionC[chB];
	r = GetAffineDecompositionM[chA] . w[[1]] . GetAffineDecompositionM[chB] + Outer[Times, rA, rB];
	lmax = GetLMax[chA, chB];

	f[chA, chB];

,{chA, (#[pA])&/@ Channels},{chB, (#[pB])&/@ Channels}]


End[];


EndPackage[];
