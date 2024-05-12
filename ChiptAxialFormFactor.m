(* ::Package:: *)

BeginPackage["ChiptAxialFormFactor`"];


PoleIntegrals::usage = "PoleIntegrals[expr,var] returns the contour integral of expression expr w.r.t. variable var according to Cauchy's integral formula";
Vap::usage = "Vap[p] pion- axial vector interaction";
Vann::usage = "Vann[\[Mu]] two nucleon - axial vector interaction";
Vpnn::usage = "Vann[p] two nucleon - pion interaction";
Vapnn::usage = "Vapnn[] two nucleon - pion - axial vector interaction";
Vppnn::usage = "Vppnn[p1,p2] two nucleon - two pion interaction";
Vappnn::usage = "Vappnn[] axial vector - two pion - 2 nucleon interaction";
Vpppnn::usage = "Vppnn[p1,p2] two nucleon - three pion interaction";
Vpn::usage = "Vpn[] nucleon - source pion interaction";
Vppn::usage = "Vppn[] nucleon - two source pion interaction";
VapnnNLO::usage = "VapnnNLO[p1,p2,k1,q] two nucleon - pion - axial vector interaction";
VppnnNLO::usage = "VppnnNLO[p1,p2,k1,k2] two nucleon - two pion interaction at NLO";
Pfermion::usage = "Pfermion[p,E,m] fermion propagator";
Ppion::usage = "Ppion[p,E] scalar pion propagator";
F\[Pi]::usage = "Pion decay constant";
gA::usage = "Axial coupling constant";
\[Tau]1::usage = "pauli matrix 1";
\[Tau]2::usage = "pauli matrix 2";
\[Tau]3::usage = "pauli matrix 3";
PreFactor::usage = "projection operator applied before taking the trace"


Begin["`Private`"];


\[Gamma][\[Mu]_] := Which[
\[Mu]==0,{{1,0,0,0},{0,1,0,0},{0,0,-1,0},{0,0,0,-1}}
,\[Mu]==1 ,{{0,0,0,1},{0,0,1,0},{0,-1,0,0},{-1,0,0,0}}
,\[Mu]==2 ,{{0,0,0,-I},{0,0,I,0},{0,I,0,0},{-I,0,0,0}}
,\[Mu]==3 ,{{0,0,1,0},{0,0,0,-1},{-1,0,0,0},{0,1,0,0}}]
\[Gamma]5 = I \[Gamma][0].\[Gamma][1].\[Gamma][2].\[Gamma][3];
\[Tau][i_] := Which[i==1,\[Tau]1,i==2,\[Tau]2,i==3,\[Tau]3]
slash[p_] := p[[1]]\[Gamma][0] - p[[2]]\[Gamma][1] - p[[3]]\[Gamma][2] - p[[4]]\[Gamma][3];
LC = LeviCivitaTensor[3];
KD[i_,j_] := KroneckerDelta[i,j]
Vap[p_,\[Mu]_] := I F\[Pi] p[[\[Mu]+1]]
Vann[\[Mu]_,i_] := -((I gA)/2)\[Gamma][\[Mu]].\[Gamma]5 \[Tau][i]
Vpnn[p_,i_] := gA/(2 F\[Pi]) slash[p].\[Gamma]5 \[Tau][i]
Vapnn[\[Mu]_,i_,j_] := 1/(2F\[Pi]) Sum[LC[[i,j,k]] \[Tau][k],{k,1,3}]\[Gamma][\[Mu]]
Vppnn[k1_,k2_,i_,j_] := I/(2F\[Pi]) Sum[LC[[i,j,k]] \[Tau][k],{k,1,3}](slash[k2]-slash[k1])
Vappnn[\[Mu]_,i_,j_,k_] := gA/(4F\[Pi]^2) (\[Tau][j].Sum[LC[[i,k,l]]\[Tau][l],{l,1,3}] + \[Tau][k].Sum[LC[[i,j,l]]\[Tau][l],{l,1,3}]) \[Gamma][\[Mu]].\[Gamma]5
Vpppnn[k1_,k2_,k3_,i_,j_,k_] := -gA/(8F\[Pi]^3) ((2 \[Tau][i]KD[j,k] - \[Tau][j].\[Tau][i].\[Tau][k] - \[Tau][k].\[Tau][i].\[Tau][j])slash[k1] + (2 \[Tau][j]KD[i,k] - \[Tau][i].\[Tau][j].\[Tau][k] - \[Tau][k].\[Tau][j].\[Tau][i])slash[k2] + (2 \[Tau][k]KD[i,j] - \[Tau][j].\[Tau][k].\[Tau][i] - \[Tau][i].\[Tau][k].\[Tau][j])slash[k3])
Vpn[i_]:= -(I/(2F))\[Tau][i] \[Gamma]5
Vppn[]:= -(1/(2F\[Pi]^2))
VppnnNLO[p1_,p2_,k1_,k2_,i_,j_] := -(c2/(m^2 F\[Pi]^2))KD[i,k](p1.k1 k2.p1 + k1.p2 k2.p2) - (2c3)/F^2 KD[i,j] (k1.k2) + (I c4)/(2F^2) Sum[LC[[i,j,k]]\[Tau][k],{k,1,3}](slash[k1].slash[k2] - slash[k2].slash[k1])
VapnnNLO[p1_,p2_,k_,q_,\[Mu]_,i_,j_] := -((I c2)/(m^2 F\[Pi]))KD[i,j](p1.k p1 + p2.k p2) - (2I c3)/F KD[i,j]k - c4/(2F) Sum[LC[[i,j,k]]\[Tau][k],{k,1,3}](slash[k].\[Gamma][\[Mu]] - \[Gamma][\[Mu]].slash[k]) - c6/(8m F) Sum[LC[[i,j,k]]\[Tau][k],{k,1,3}](slash[q].\[Gamma][\[Mu]] - \[Gamma][\[Mu]].slash[q])
Pfermion[p_,e_,m_] := (slash[p] + m IdentityMatrix[4])/(p[[1]]^2 - e^2)
Ppion[p_,e_] := 1/(p[[1]]^2 - e^2)
PreFactor[] := ((IdentityMatrix[4]+\[Gamma][0])/2).(IdentityMatrix[4] + \[Gamma][3].\[Gamma]5)


End[];


EndPackage[];
