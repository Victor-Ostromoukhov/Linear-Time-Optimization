(* Linear-Time-Optimization
	vo + Raanan Fattal May 2017
*)


(*myShiftPtsFunction := shiftPtsBitXor;*)			(* shiftPtsMixed or shiftPtsBitXor *)
myShiftPtsFunction := shiftPtsCranleyPatterson ;
myShiftPtsFunction := shiftPtsMixed;			(* shiftPtsMixed or shiftPtsBitXor *)
myShiftPtsFunction := shiftPtsBitXor;

(**************** System-dependent setup ******************)
SetDirectory[ToFileName[$HomeDirectory,"Linear-Time-Optimization"]];
SetOptions[Graphics, ImageSize -> { 1/4 1024,Automatic},AspectRatio->Automatic, PlotRange->All];
SetOptions[ListPlot, ImageSize -> { 600,Automatic},AspectRatio->.2, PlotRange->All];
SetOptions[ListLogLogPlot, ImageSize -> { 600,Automatic},AspectRatio->.61, PlotRange->All];
SetOptions[ListDensityPlot, ImageSize -> {512,Automatic},AspectRatio->Automatic, PlotRange->All];
pid := "_pid"<>ToString[$ProcessID]<>"_kid"<>ToString[$KernelID]
systemID = StringSplit[$System][[1]]; (* "Mac" or "Linux" *)
execPrefix = "~/bin/";

second[x_]:= If[Length[x] > 1, x[[2]], First[x] ] (* like First *)
third[x_]:= If[Length[x] > 2, x[[3]], First[x] ] (* like First *)
fourth[x_]:= If[Length[x] > 3, x[[4]], First[x] ] (* like First *)

epsilon = 10^-10.;
eps = 10^-6.;

Print["LinearTimeOptimization loaded."];


(**************** useful common constants, macros and routines ******************)
mf := MatrixForm
if := InputForm
T :=  Transpose
gr = GoldenRatio // N
PI = Pi//N;
known := ValueQ

n2PaddedString[n_,len_:5] := ToString[NumberForm[n, len-1, NumberPadding -> "0"]]

euclidlen[z_] := Sqrt[Total[z^2]]
euclidlen2[z_] := Total[z^2]
euclidlenN[z_] := Sqrt[Total[z^2]]//N

arr3D[{from_,to_},th_:.03] :={Cylinder[{from,.1from + .9 to},th],Cone[{.1from + .9 to,to},2 th]}

(*Complement*)
colSolidTable = {Red,Blue,Orange,Gray,Magenta,Yellow,Cyan,Green};
getColorSolid[ind_] := colSolidTable[[Mod[ind,Length[colSolidTable],1] ]]

order2permut[s_] := (Flatten[Drop[#,1]& /@ Sort[Table[{s[[i]],i},{i,Length[s]}]]]); (* 1..n *)
permut2order := order2permut

order2permut1toN := order2permut
order2permut0toNminus1[s_] := ((s+1)//order2permut)-1

vdc1dBase2[nbits_:4]    := (FromDigits[#,2]& /@ (Reverse /@ (IntegerDigits[#,2,nbits]& @ Range[0,2^nbits-1])))/2^nbits
vdc1dBase2Int[nbits_:4] := (FromDigits[#,2]& /@ (Reverse /@ (IntegerDigits[#,2,nbits]& @ Range[0,2^nbits-1]))) 

vdc1dBaseN[nbits_:4,base_:3] := (FromDigits[#,base]& /@ (Reverse /@ (IntegerDigits[#,base,nbits]& @  Range[0,base^nbits-1])))/base^nbits

reorg2D[lst2D_] :=
    With[ {sz = Floor[Length[lst2D]/2]},
        If[ sz == 0, Print["reorg3D: void data"]; Break[] ];
        (RotateLeft[#,sz]& /@ (RotateLeft[#,sz]& @ lst2D))
    ] (* reorg2D *)
Reorg := reorg2D

reorg3D[lst3D_] :=
    With[ {sz = Floor[Length[lst3D]/2]},
        If[ sz == 0, Print["reorg3D: void data"]; Break[] ];
        (RotateLeft[#,sz]& /@ #)& /@ (RotateLeft[#,sz]& /@ (RotateLeft[#,sz]& @ lst3D))
    ] (* reorg3D *)

reorg4D[lst4D_] :=
    With[ {sz = Floor[Length[lst4D]/2]},
        If[ sz == 0, Print["reorg4D: void data"]; Break[] ];
        (RotateLeft[#,sz]& /@ #)& /@   ((RotateLeft[#,sz]& /@ #)& /@ (RotateLeft[#,sz]& /@ (RotateLeft[#,sz]& @ lst4D)))
    ] (* reorg4D *)

getFourier1D[xvals_,fouriertabsz_] :=
    Module[ {tab,pos},
        tab = Table[0,{fouriertabsz}];
        If[ Min[xvals] <= 0,
            Print["getFourier1D: min=",Min[xvals]];
            Abort[]
        ];
        If[ Max[xvals] > fouriertabsz,
            Print["getFourier1D: max=",Max[xvals]];
            Abort[]
        ];
        Do[
            pos = xvals[[i]];
            If[ (pos <= 0) || (pos > fouriertabsz),
                Continue[]
            ];
            tab[[ pos ]] += 1
        ,{i,Length[xvals]}];
        (*tab[[xvals+1]] = 1;*)
        Return[Fourier[tab]// Abs]
    ] (* getFourier1D *)

getFourier2D[pts_,fouriertabsz_:512] :=
    Module[ {tab,sel},
        tab = Table[0,{fouriertabsz},{fouriertabsz}];
        sel = Select[Ceiling[pts], Min[#] > 0 && Max[#] <= fouriertabsz &];
        (tab[[#[[2]], #[[1]]]] = 1) & /@ sel;
        Return[Fourier[tab]// Abs]
    ] (* getFourier2D *)

getFourier2DPlusMorror[pts_,fouriertabsz_:512] :=
    Module[ {tab,sel},
        tab = Table[0,{fouriertabsz},{fouriertabsz}];
        sel = Select[Ceiling[pts], Min[#] > 0 && Max[#] <= fouriertabsz &];
        (tab[[#[[2]], #[[1]]]] = 1) & /@ sel;
        Return[
        	(
        	(Abs@Fourier[tab]) + (Abs@Fourier[Reverse/@tab]) + (Abs@Fourier[Reverse@tab]) + (Abs@Fourier[Reverse@(Reverse/@tab)]) +
        	(Abs@Fourier[T@tab]) + (Abs@Fourier[T@(Reverse/@tab)]) + (Abs@Fourier[T@(Reverse@tab)]) + (Abs@Fourier[T@(Reverse@(Reverse/@tab))])
        	)/8
        	]
    ] (* getFourier2DPlusMorror *)

getImageAndFourier2D[pts_,fouriertabsz_:512] :=
    Module[ {tab,sel},
        tab = Table[0,{fouriertabsz},{fouriertabsz}];
        sel = Select[Ceiling[pts], Min[#] > 0 && Max[#] <= fouriertabsz &];
        (tab[[#[[2]], #[[1]]]] = 1) & /@ sel;
        Return[{tab,Fourier[tab]// Abs}]
    ] (* getImageAndFourier2D *)

getFourier3D[pts_,fouriertabsz_:256] :=
    Module[ {tab,sel},
        tab = Table[0,{fouriertabsz},{fouriertabsz},{fouriertabsz}];
        sel = Select[Ceiling[pts], Min[#] > 0 && Max[#] <= fouriertabsz &];
        (tab[[#[[3]], #[[2]], #[[1]]]] = 1) & /@ sel;
        Return[Fourier[tab]// Abs]
    ] (* getFourier3D *)

getFourier4D[pts_,fouriertabsz_:64] :=
    Module[ {tab,sel},
        tab = Table[0,{fouriertabsz},{fouriertabsz},{fouriertabsz},{fouriertabsz}];
        sel = Select[Ceiling[pts], Min[#] > 0 && Max[#] <= fouriertabsz &];
        (tab[[#[[4]], #[[3]], #[[2]], #[[1]]]] = 1) & /@ sel;
        Return[Fourier[tab]// Abs]
    ] (* getFourier4D *)

niceRaster[img_,OptionsPattern[]] :=
    Block[ {sx,sy,z,lbl},
        z = OptionValue[zoom];
        lbl = OptionValue[PlotLabel];
        {sy,sx} = Take[#, 2] & @ Dimensions[img];
        Return[Graphics[Raster[img],PlotRange->{{0,sx},{0,sy}},ImageSize->{z sx,z sy},PlotLabel->lbl]];
    ];
Options[niceRaster] = {zoom->1,PlotLabel->None};

get2DfourierAndRadial4K[pts4K_,lbl_:""] := (* pts4K must be between 0 and 1 *)
    Module[ {fouriertabsz=4096,fsum,centralSz=128,rfsum,fs,fs1,maxval,maxpositions,rgbimage,iy,ix,radialSpectrum,radialSpectrumMax,
    	maxradialSpectrum,maxradialSpectrumMax,radialSpectrumScrambledLDS2dPureRandom,p1,p2,p3},
        radialSpectrumScrambledLDS2dPureRandom = {{1,4.29518*10^-8},{2,7.79978*10^-8},{3,9.82786*10^-8},{4,2.14975*10^-7},{5,3.89991*10^-7},{6,5.65343*10^-7},{7,7.58657*10^-7},{8,1.19737*10^-6},{9,3.06498*10^-6},{10,4.0001*10^-6},{11,6.12319*10^-6},{12,7.05342*10^-6},{13,7.81769*10^-6},{14,0.0000111474},{15,0.0000151127},{16,0.0000206135},{17,0.0000224409},{18,0.0000224069},{19,0.0000325448},{20,0.0000429995},{21,0.0000449287},{22,0.0000448171},{23,0.0000583691},{24,0.0000538246},{25,0.0000515352},{26,0.000053307},{27,0.0000571127},{28,0.0000724713},{29,0.0000915756},{30,0.0000918046},{31,0.000100565},{32,0.000124037},{33,0.000108819},{34,0.000106498},{35,0.000130497},{36,0.000138639},{37,0.000145227},{38,0.000175364},{39,0.000166251},{40,0.000176758},{41,0.000168615},{42,0.000165452},{43,0.000157801},{44,0.00015592},{45,0.000180432},{46,0.000170291},{47,0.000183831},{48,0.000197654},{49,0.000190239},{50,0.000205299},{51,0.000196384},{52,0.000199393},{53,0.000184867},{54,0.000199862},{55,0.00019525},{56,0.000183077},{57,0.000227733},{58,0.000232363},{59,0.000223094},{60,0.0001988},{61,0.000178908},{62,0.000201434},{63,0.00021009},{64,0.000208837},{65,0.00019644},{66,0.000182541},{67,0.000187736},{68,0.000190628},{69,0.000188631},{70,0.000197614},{71,0.000211417},{72,0.000204345},{73,0.000177879},{74,0.000169891},{75,0.000182841},{76,0.000200353},{77,0.000200096},{78,0.000192625},{79,0.000195345},{80,0.000200256},{81,0.000201079},{82,0.000216014},{83,0.000201631},{84,0.000198567},{85,0.000203169},{86,0.000205037},{87,0.000203988},{88,0.000199763},{89,0.000192393},{90,0.000186826},{91,0.000200103},{92,0.000231857},{93,0.000224188},{94,0.000218426},{95,0.0002249},{96,0.000220904},{97,0.000220412},{98,0.00021262},{99,0.000192152},{100,0.000188746},{101,0.000220298},{102,0.000218039},{103,0.000196443},{104,0.000189949},{105,0.000210372},{106,0.000226739},{107,0.000209216},{108,0.000201318},{109,0.000197378},{110,0.000182891},{111,0.000184064},{112,0.000199954},{113,0.000204013},{114,0.000228216},{115,0.000225159},{116,0.000207721},{117,0.000212594},{118,0.00019653},{119,0.000189115},{120,0.000194733}};
		fFourier := getFourier2D; (* getFourier2D or getFourier2DPlusMorror *)
        fsum = (Chop @ fFourier[Floor[1+ fouriertabsz pts4K], fouriertabsz]);
        fsum[[1,1]] = 0;
        rfsum = reorg2D[fsum];
        fs = rfsum[[fouriertabsz/2+1 - centralSz ;; fouriertabsz/2+1 + centralSz, fouriertabsz/2+1 - centralSz ;; fouriertabsz/2+1 + centralSz ]];
        fs1 = T[T[fs[[32;;225]]][[32;;225]]]; (* area [-1.5,11.5] *)
        maxval = Max[fs1];
        maxpositions = Position[fs1, Max[fs1]];
        fs1[[98,98]] = 1;
        rgbimage = Partition[#,Length[fs1]]& @ ({Flatten@fs1,Flatten@fs1,Flatten@fs1}//T);
        ({ix,iy} = #; rgbimage[[iy,ix]] = {1,0,0}) & /@ maxpositions;
        {radialSpectrum,radialSpectrumMax} = radialPowerSpectrum[fs^2];
        maxradialSpectrum = Max[radialSpectrum];
        maxradialSpectrumMax = Max[radialSpectrumMax];
        radialSpectrum = {Range[Length[radialSpectrum]],radialSpectrum}//T;
        radialSpectrumMax = {Range[Length[radialSpectrumMax]],radialSpectrumMax}//T;
		{p1,p2,p3} = {
		Graphics[{PointSize[.005],Point/@ pts4K}, ImageSize -> {256,256}, PlotLabel->lbl],
		niceRaster[30^2 rgbimage^2,zoom->2],
             ListPlot[{radialSpectrum[[;;120]],radialSpectrumScrambledLDS2dPureRandom[[;;120]],{{67,0},{67,.0004}},{{0,.0002},{120,.0002}}},
             PlotStyle->{Red,Blue,Cyan,Cyan},Joined->True,AspectRatio->256/512,Ticks->None, ImageSize -> {300, 300},PlotLabel->"Red: ours Blue:ScrambledLDS2dPureRandom"]};
        {p1,p2,p3}
    ] (* get2DfourierAndRadial *)

get2DfourierAndRadial16K[pts4K_,lbl_:""] := (* pts4K must be between 0 and 1 *)
    Module[ {fouriertabsz=2 4096,fsum,centralSz=2 128,rfsum,fs,fs1,maxval,maxpositions,rgbimage,iy,ix,radialSpectrum,radialSpectrumMax,
    	maxradialSpectrum,maxradialSpectrumMax,p1,p2,p3},
		fFourier := getFourier2D; (* getFourier2D or getFourier2DPlusMorror *)
        fsum = (Chop @ fFourier[Floor[1+ fouriertabsz pts4K], fouriertabsz]);
        fsum[[1,1]] = 0;
        rfsum = reorg2D[fsum];
        fs = rfsum[[fouriertabsz/2+1 - centralSz ;; fouriertabsz/2+1 + centralSz, fouriertabsz/2+1 - centralSz ;; fouriertabsz/2+1 + centralSz ]];
        (*fs1 = T[T[fs[[32;;225]]][[32;;225]]]; (* area [-1.5,11.5] *)*)
        fs1 = fs;
        maxval = Max[fs1];
        maxpositions = Position[fs1, Max[fs1]];
        (*fs1[[98,98]] = 1;*)
        rgbimage = Partition[#,Length[fs1]]& @ ({Flatten@fs1,Flatten@fs1,Flatten@fs1}//T);
        ({ix,iy} = #; rgbimage[[iy,ix]] = {1,0,0}) & /@ maxpositions;
        {radialSpectrum,radialSpectrumMax} = radialPowerSpectrum[fs^2];
        maxradialSpectrum = Max[radialSpectrum];
        maxradialSpectrumMax = Max[radialSpectrumMax];
        radialSpectrum = {Range[Length[radialSpectrum]],radialSpectrum}//T;
        radialSpectrumMax = {Range[Length[radialSpectrumMax]],radialSpectrumMax}//T;
		{p1,p2,p3} = {
		Graphics[{PointSize[.005],Point/@ pts4K}, ImageSize -> {256,256}, PlotLabel->lbl],
		niceRaster[30^2 rgbimage^2,zoom->2],
             ListPlot[{radialSpectrum[[;;2 120]],(*{{2 67,0},{2  67,.0004}},*){{0,.0002},{2 120,.0002}}},
             PlotStyle->{Red,Cyan,Cyan},Joined->True,AspectRatio->256/512,Ticks->None, ImageSize -> {300, 300}]};
        {radialSpectrum[[;;2 120]],p2,p3}
    ] (* get2DfourierAndRadial16K *)

radialPowerSpectrum[fpower_,dbg_:False] :=
    Module[ {fouriertabsz = Length[fpower], binsScalingFactor = 1, bins, r, radialSpectrum,radialSpectrumMax},
        bins = Table[{{},0},{ fouriertabsz}];
        Do[
            Do[
                r = Round[binsScalingFactor (euclidlen[{i,j}-{fouriertabsz/2+1,fouriertabsz/2+1}])];
                If[ r > 0,
                    AppendTo[bins[[r,1]],fpower[[j,i]] ];
                    bins[[r,2]]++;
                ];
            ,{i,fouriertabsz}]
        ,{j,fouriertabsz}];
        radialSpectrum = Table[If[ bins[[i, 2]] == 0, 0, Total[bins[[i, 1]]]/bins[[i, 2]] // N ], {i, Length[bins]}];
        radialSpectrumMax = Table[If[ bins[[i, 2]] == 0, 0, Max[bins[[i, 1]]]  // N ], {i, Length[bins]}];
        {radialSpectrum,radialSpectrumMax}
     ] (* radialPowerSpectrum *)

getDiscrepancy2Dexact[pts_] :=
    Module[ {execString},
        Export["tmp/tmp"<>pid<>".dat",N[pts]];
        execString =  "discrepancy -i tmp/tmp"<>pid<>".dat -o tmp/res"<>pid<>".dat -d 2 > /dev/null"; (* Linux: OpenMP version; Mac: non-parallel version *)
        (*Print[execString];*)
        Run[execPrefix<>execString];
        Last @ (Flatten@Import["tmp/res"<>pid<>".dat"])
    ] (* getDiscrepancy2Dexact *)

getGL2DDiscreapancy[dtab_,imagesize_:{1100,700},ourLabel_:"Ours"] :=
    Module[ {(*refpow05,refpow1LogSMinus1Halved,dtrefpow1LogSab,coef,npts,val,pDiscrepancy,discrepancyTabSobol,discrepancyHalton,discrepancyStratified,discrepancyPoissonDisk,discrepancyWhiteNoise,s*)},
        discrepancyWhiteNoise =  Get["data_discrepancy/discrepancyWhiteNoise.dat"];
        discrepancyStratified = Get["data_discrepancy/discrepancyStratified.dat"];
        discrepancyPoissonDisk = Get["data_discrepancy/discrepancyPoissonDisk.dat"];
        (*discrepancyTabSobol = Get["data_discrepancy/discrepancyTabSobol.dat"];
        discrepancyHalton = Get["data_discrepancy/discrepancyHalton.dat"];*)

        discrepancyTabSobol = Get["data_discrepancy/discrepancyTabSobol_short_step1over8.dat"];
        discrepancyHalton = Get["data_discrepancy/discrepancyHalton_short_step1over8.dat"];
        s=2;
        coef = 2.2;
        refpow05 = Table[
            npts = Round[2^(i/8)];
            val = coef npts^-.5;
            {npts,val}
        ,{i,3*8, 20*8}];
        coef = 1/2^(4 s) 1/((s - 1) Log[10,2])^((s - 1)/2) // N; (* Roth's consta *)
        coef = .25;
        refpow1LogS = Table[
            npts = Round[2^(i/8)];
            val = coef (Log[npts])^s npts^-1.0;
            {npts,val}
        ,{i,15*8, 20*8}];
        coef = .5;
        refpow1LogSMinus1Halved = Table[
            npts = Round[2^(i/8)];
            val = coef (Log[npts])^((s-1)/2) npts^-1.0;
            {npts,val}
        ,{i,15*8, 20*8}];
        k = .;
        pDiscrepancy = ListPlot[{
		            Log[10,#]&  /@ refpow05,
		            Log[10,#]&  /@ discrepancyWhiteNoise,
		            Log[10,#]&  /@ discrepancyPoissonDisk, 
		            Log[10,#]&  /@ discrepancyStratified, 
		            Log[10,#]& /@ discrepancyHalton, 
		            Log[10,#]& /@ discrepancyTabSobol,
		            Log[10,#]&  /@ dtab,
		            Log[10,#]&  /@ refpow1LogS,
		            Log[10,#]&  /@ refpow1LogSMinus1Halved
		        }
            , PlotLegends ->  Placed[#,{.25,.4}]& @ {
                Style[#,30]& @ (Subscript[k,1]  / N^Style["1/2",Italic]  ),
                Style[#,30]& @"WhiteNoise",
                Style[#,30]& @"PoissonDisk",
                Style[#,30]& @"Jitter",
                Style[#,30]& @"Halton",
                Style[#,30]& @"Sobol2dGrayCode_1_2",
                Style[#,30]& @ ourLabel,
                (*Style[#,30]& @ (Subscript[k,2] Style["(Log N)",Italic]  / N^Style["1",Italic] )*)
                Style[#,30]& @ (Subscript[k,2] Style["(Log N)"^Style["s",Italic],Italic] / N^Style["1",Italic] ),
                Style[#,30]& @ (Subscript[k,3] Style["(Log N)"^Style["(s-1)/2",Italic],Italic] / N^Style["1",Italic])
                }
            ,PlotRange->{{.5,6.03},Automatic}
            (*,FrameTicks->{{N[# Log[4]/Log[10]],ToString[#] }& /@ Range[12] ,Automatic}*)
            ,FrameTicks->{Automatic,Automatic}
            ,AspectRatio->.61
            ,FrameLabel-> {Style[ HoldForm@(Subscript[Log, 10] "(NSamples)"), 36],Style[ HoldForm@(Subscript[Log, 10] "(Discrepancy)"), 36] }
            ,FrameStyle->Directive[Black,24]
            ,RotateLabel -> True
            ,PlotMarkers->{{\[FilledCircle],1},{\[FilledCircle],10},{\[FilledCircle],10},{\[FilledCircle],10},{\[FilledCircle],10},{\[FilledCircle],10},{\[FilledCircle],10},{\[FilledCircle],1}}
            ,Frame->True
            ,ImageSize -> imagesize
            ,PlotStyle->  {{Black,AbsoluteDashing[{10,5}]},Black,Orange,Green,Blue,Gray,{Thickness[.003],Red},{Red,Dotted},{Blue,Dotted}}
            , Joined->True
            ,PlotRange->All
            ,PlotLabel-> Style["log-log 2D Discrepancy",Black,42]
            ];
        (*pDiscrepancy//Print;*)
        (*ListPlot[({Log[#[[1]]],Log[ #[[1]]/Log[#[[1]]] #[[2]]]} & /@ #) & /@ {refpow05,discrepancyWhiteNoise,discrepancyPoissonDisk,discrepancyStratified,discrepancyHalton,discrepancyTabSobol,dtab,refpow10}
        	, Joined -> True, ImageSize -> {1100, 700}, AspectRatio->.61, PlotStyle->{{Black,AbsoluteDashing[{10,5}]},Black,Orange,Green,Blue,Gray,Red,{Black,Dotted}}]//Print;*)
        pDiscrepancy
    ] (* getGL2DDiscreapancy *)

(*-------------------------------------------------------------------*)

exploreDataDiscreapancy[] :=
    Module[ {},
		(* D* *)
		names = {"pts_1092_23.dat","pts_17476_15.dat","pts_279620_24.dat","pts_4369_17.dat","pts_69905_15.dat"};
		names = {"pts_17476.dat","pts_279620.dat","pts_69905.dat"};
		counts = {1092,4369,17476,69905,279620};
		dtab = Sort @ Table[
			name = names[[iname]];
			fname = "data/sets_20170508/"<>name;
			fname = "data/sets_20170503/"<>name;
			pts = Import[fname];
			npts = Length[pts];
			d = getDiscrepancy2Dexact[pts];
			Print["Processing "fname-> {npts,d} ];
			{npts,d}
		,{iname,Length[names]}];
		Print[dtab];
		
		dtab = {{1092, 0.013725}, {4369, 0.00316811}, {17476, 0.00169935}, {69905, 0.000531731}, {279620, 0.000210969}};
		dtab = {{17476, 0.00158537},{69905, 0.000450239},{279620, 0.000191513}};
		getGL2DDiscreapancy[dtab,{1100,700},"LinearTimeOptim May3"]//Print;
   ]
   

exploreDataFourier[] :=
    Module[ {},
		(* D* *)
		names = {"pts_1092_23.dat","pts_17476_15.dat","pts_279620_24.dat","pts_4369_17.dat","pts_69905_15.dat"};
		names = {"pts_17476.dat","pts_279620.dat","pts_69905.dat"};
		counts = {1092,4369,17476,69905,279620};
		(* Fourier *)
    	fouriertabsz=2 4096;
    	centralSz=2 128;
    	fname = "data/sets_20170508/pts_69905_15.dat";
    	fname = "data/sets_20170508/pts_4369_17.dat";
    	
    	fname = "data/sets_20170509/pts_17476_A.dat";
        pts = Import[fname];
        Graphics[{Point/@pts},ImageSize -> 1/2{ 1024,1024}]//Print;
            {rPS1A,p2,p3} = get2DfourierAndRadial16K[pts,"LinearTimeOptim"];
            {p2,p3}//Print;

    	fname = "data/sets_20170509/pts_17476_B.dat";
        pts = Import[fname];
        Graphics[{Point/@pts},ImageSize -> 1/2{ 1024,1024}]//Print;
            {rPS1B,p2,p3} = get2DfourierAndRadial16K[pts,"LinearTimeOptim"];
            {p2,p3}//Print;

    	fname = "data/sets_20170503/pts_17476.dat";
        pts = Import[fname];
        Graphics[{Point/@pts},ImageSize -> 1/2{ 1024,1024}]//Print;
            {rPS2,p2,p3} = get2DfourierAndRadial16K[pts,"LinearTimeOptim"];
            {p2,p3}//Print;

		
        ListPlot[{
		            Log[10,#]&  /@ rPS1A,
		            Log[10,#]&  /@ rPS1B,
		            Log[10,#]&  /@ rPS2
		        },AspectRatio->1,PlotStyle->{Red,Blue,Black},PlotLabel->"Red:May9A Blue: May9B, Black:May3"]//Print;

   ]