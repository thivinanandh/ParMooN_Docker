cl__1 = 1e+22;
Point(1) = {0, 0, 18.3193878724927, cl__1};
Point(2) = {0, 58.0822473540764, 18.3193878724927, cl__1};
Point(3) = {0, 58.0822473540764, 0, cl__1};
Point(4) = {0, 0, 0, cl__1};
Point(5) = {152.90077690686, 58.0822473540764, 18.3193878724927, cl__1};
Point(6) = {152.90077690686, 58.0822473540764, 0, cl__1};
Point(7) = {152.90077690686, 0, 18.3193878724927, cl__1};
Point(8) = {152.90077690686, 0, 0, cl__1};
Point(9) = {119.666720578679, 29.0411236770382, 18.3193878724926, cl__1};
Point(10) = {25.9506160449542, 24.6145593465784, 18.3193878724926, cl__1};
Point(11) = {25.9506160449542, 33.467688007498, 18.3193878724927, cl__1};
Point(12) = {25.9506160449542, 29.0411236770382, 15.9940340390435, cl__1};

//End of outer Boundary lines for rectangular Component
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {4,3};
Line(4) = {1,4};
Line(5) = {2,5};
Line(6) = {5,6};
Line(7) = {3,6};
Line(8) = {5,7};
Line(9) = {7,8};
Line(10) = {6,8};
Line(11) = {7,1};
Line(12) = {8,4};

// End of outer Boundary lines for rectangular Component

 /////////////////////////////// Inner Ship Cavity Component ///////////////////////////////////////

//// Forward Ship water Line  ///////
p13 = newp;
Point(p13 + 1) = {117.9133965981383, 27.46790253563098, 18.31938787249264};
Point(p13 + 2) = {115.7660399263147, 26.34414698572264, 18.31938787249268};
Point(p13 + 3) = {113.5836966862716, 25.30599925190943, 18.31938787249265};
Point(p13 + 4) = {111.3897643920874, 24.28430363445318, 18.3193878724927};
Point(p13 + 5) = {109.1823274724141, 23.28893126380683, 18.31938787249264};
Point(p13 + 6) = {106.9426068891197, 22.35513489706123, 18.3193878724927};
Point(p13 + 7) = {104.6588378768967, 21.54322511655446, 18.31938787249269};
Point(p13 + 8) = {102.3392518183817, 20.89930521261095, 18.3193878724927};
Point(p13 + 9) = {99.96896737708012, 20.39297030138071, 18.31938787249263};
Point(p13 + 10) = {97.57751179490276, 20.00310998208104, 18.31938787249268};
Point(p13 + 11) = {95.17510434524844, 19.71697176106764, 18.31938787249263};
Point(p13 + 12) = {92.76833027564723, 19.52479872322295, 18.31938787249268};
Point(p13 + 13) = {90.34775963009206, 19.41110647836308, 18.3193878724927};
Point(p13 + 14) = {87.92877944909321, 19.3501893593509, 18.3193878724927};
Point(p13 + 15) = {85.50963101946979, 19.31865186295023, 18.31938787249267};
Point(p13 + 16) = {83.09033653876381, 19.30278277192399, 18.31938787249267};
Point(p13 + 17) = {80.67100393777014, 19.2952951983593, 18.31938787249266};
Point(p13 + 18) = {78.25166465116145, 19.29230878359422, 18.31938787249262};
Point(p13 + 19) = {75.83232451099001, 19.29153309833901, 18.3193878724927};
Point(p13 + 20) = {73.41222265412384, 19.29131181486977, 18.31938787249261};
Point(p13 + 21) = {70.99205293385624, 19.29118397834453, 18.3193878724927};
Point(p13 + 22) = {68.57188321131399, 19.29112997123665, 18.3193878724927};
Point(p13 + 23) = {66.15171348873022, 19.29112375200618, 18.31938787249268};
Point(p13 + 24) = {63.73154376808678, 19.29117243807257, 18.3193878724927};
Point(p13 + 25) = {61.31137405633939, 19.29137906759086, 18.31938787249278};
Point(p13 + 26) = {58.89120437811186, 19.29181791293745, 18.3193878724927};
Point(p13 + 27) = {56.47103479985163, 19.29262353945816, 18.3193878724927};
Point(p13 + 28) = {54.05086553399599, 19.29407951260235, 18.3193878724927};
Point(p13 + 29) = {51.6306974004311, 19.29680357044294, 18.31938787249262};
Point(p13 + 30) = {49.21582599975046, 19.30222134139602, 18.31938787249267};
Point(p13 + 31) = {46.80228954926865, 19.31391541058027, 18.31938787249267};
Point(p13 + 32) = {44.38891462006118, 19.3426980251507, 18.3193878724927};
Point(p13 + 33) = {41.95345087735404, 19.43229026710718, 18.31938787249262};
Point(p13 + 34) = {39.55462567842417, 19.74693439376802, 18.3193878724926};
Point(p13 + 35) = {37.22539261577334, 20.39850175940919, 18.3193878724927};
Point(p13 + 36) = {34.98581520948403, 21.3113918902262, 18.31938787249262};
Point(p13 + 37) = {32.79878862175321, 22.34798038072643, 18.31938787249262};
Point(p13 + 38) = {30.59351063532331, 23.34100599439493, 18.31938787249269};
Point(p13 + 39) = {28.32126022151572, 24.17755723187089, 18.31938787249268};

Spline(13) = {9, p13 + 1, p13 + 2, p13 + 3, p13 + 4, p13 + 5, p13 + 6, p13 + 7, p13 + 8, p13 + 9, p13 + 10, p13 + 11, p13 + 12, p13 + 13, p13 + 14, p13 + 15, p13 + 16, p13 + 17, p13 + 18, p13 + 19, p13 + 20, p13 + 21, p13 + 22, p13 + 23, p13 + 24, p13 + 25, p13 + 26, p13 + 27, p13 + 28, p13 + 29, p13 + 30, p13 + 31, p13 + 32, p13 + 33, p13 + 34, p13 + 35, p13 + 36, p13 + 37, p13 + 38, p13 + 39, 10};

//// Forward Ship water Line  - End  ///////

/// Back Water Line - Begin ///
p15 = newp;
Point(p15 + 1) = {117.91339659814, 30.61434481844445, 18.31938787249266};
Point(p15 + 2) = {115.7660399263209, 31.73810036835106, 18.31938787249268};
Point(p15 + 3) = {113.5836966862713, 32.77624810216712, 18.31938787249265};
Point(p15 + 4) = {111.389764392087, 33.7979437196234, 18.3193878724927};
Point(p15 + 5) = {109.1823274724139, 34.79331609026963, 18.31938787249265};
Point(p15 + 6) = {106.9426068891184, 35.72711245701564, 18.3193878724927};
Point(p15 + 7) = {104.6588378768959, 36.53902223752218, 18.31938787249269};
Point(p15 + 8) = {102.3392518183815, 37.18294214146557, 18.31938787249268};
Point(p15 + 9) = {99.9689673770792, 37.68927705269581, 18.3193878724926};
Point(p15 + 10) = {97.57751179490133, 38.07913737199559, 18.31938787249268};
Point(p15 + 11) = {95.17510434524638, 38.365275593009, 18.3193878724927};
Point(p15 + 12) = {92.7683302756469, 38.55744863085354, 18.3193878724927};
Point(p15 + 13) = {90.34775963009318, 38.67114087571332, 18.31938787249261};
Point(p15 + 14) = {87.92877944909424, 38.73205799472549, 18.3193878724927};
Point(p15 + 15) = {85.50963101947063, 38.76359549112617, 18.3193878724927};
Point(p15 + 16) = {83.09033653876443, 38.77946458215241, 18.31938787249263};
Point(p15 + 17) = {80.67100393777048, 38.7869521557171, 18.3193878724927};
Point(p15 + 18) = {78.25166465116166, 38.78993857048219, 18.3193878724927};
Point(p15 + 19) = {75.83232451099001, 38.79071425573742, 18.3193878724927};
Point(p15 + 20) = {73.41222265412378, 38.79093553920663, 18.3193878724926};
Point(p15 + 21) = {70.99205293385631, 38.79106337573187, 18.3193878724927};
Point(p15 + 22) = {68.571883211314, 38.79111738283974, 18.3193878724927};
Point(p15 + 23) = {66.15171348873018, 38.79112360207022, 18.31938787249268};
Point(p15 + 24) = {63.73154376808674, 38.79107491600383, 18.3193878724926};
Point(p15 + 25) = {61.31137405633949, 38.79086828648555, 18.31938787249275};
Point(p15 + 26) = {58.89120437811195, 38.79042944113895, 18.3193878724926};
Point(p15 + 27) = {56.47103479985167, 38.78962381461824, 18.3193878724927};
Point(p15 + 28) = {54.05086553399602, 38.78816784147406, 18.3193878724927};
Point(p15 + 29) = {51.63069740043104, 38.78544378363346, 18.31938787249262};
Point(p15 + 30) = {49.2158259997533, 38.78002601268038, 18.31938787249267};
Point(p15 + 31) = {46.80228954927511, 38.76833194349625, 18.31938787249268};
Point(p15 + 32) = {44.38891462007138, 38.73954932892589, 18.3193878724927};
Point(p15 + 33) = {41.95345087735408, 38.64995708696922, 18.31938787249262};
Point(p15 + 34) = {39.55462567842594, 38.33531296030879, 18.3193878724926};
Point(p15 + 35) = {37.22539261577451, 37.68374559466761, 18.31938787249265};
Point(p15 + 36) = {34.98581520948091, 36.7708554638488, 18.3193878724927};
Point(p15 + 37) = {32.79878862175263, 35.7342669733497, 18.3193878724927};
Point(p15 + 38) = {30.59351063532714, 34.74124135968306, 18.3193878724927};
Point(p15 + 39) = {28.32126022152166, 33.90469012220738, 18.31938787249268};
Spline(15) = {9, p15 + 1, p15 + 2, p15 + 3, p15 + 4, p15 + 5, p15 + 6, p15 + 7, p15 + 8, p15 + 9, p15 + 10, p15 + 11, p15 + 12, p15 + 13, p15 + 14, p15 + 15, p15 + 16, p15 + 17, p15 + 18, p15 + 19, p15 + 20, p15 + 21, p15 + 22, p15 + 23, p15 + 24, p15 + 25, p15 + 26, p15 + 27, p15 + 28, p15 + 29, p15 + 30, p15 + 31, p15 + 32, p15 + 33, p15 + 34, p15 + 35, p15 + 36, p15 + 37, p15 + 38, p15 + 39, 11};

/// Back Water Line - end ///

/// Back triangle Top part - begin ////

p14 = newp;
Point(p14 + 1) = {25.9506160449542, 24.83588756310139, 18.31938787249261};
Point(p14 + 2) = {25.9506160449542, 25.05721577962438, 18.31938787249261};
Point(p14 + 3) = {25.9506160449542, 25.27854399614737, 18.31938787249261};
Point(p14 + 4) = {25.9506160449542, 25.49987221267036, 18.31938787249261};
Point(p14 + 5) = {25.9506160449542, 25.72120042919335, 18.31938787249262};
Point(p14 + 6) = {25.9506160449542, 25.94252864571634, 18.31938787249262};
Point(p14 + 7) = {25.9506160449542, 26.16385686223933, 18.31938787249262};
Point(p14 + 8) = {25.9506160449542, 26.38518507876232, 18.31938787249262};
Point(p14 + 9) = {25.9506160449542, 26.60651329528531, 18.31938787249262};
Point(p14 + 10) = {25.9506160449542, 26.8278415118083, 18.31938787249263};
Point(p14 + 11) = {25.9506160449542, 27.04916972833129, 18.31938787249263};
Point(p14 + 12) = {25.9506160449542, 27.27049794485428, 18.31938787249263};
Point(p14 + 13) = {25.9506160449542, 27.49182616137727, 18.31938787249263};
Point(p14 + 14) = {25.9506160449542, 27.71315437790026, 18.31938787249264};
Point(p14 + 15) = {25.9506160449542, 27.93448259442325, 18.31938787249264};
Point(p14 + 16) = {25.9506160449542, 28.15581081094624, 18.31938787249264};
Point(p14 + 17) = {25.9506160449542, 28.37713902746923, 18.31938787249264};
Point(p14 + 18) = {25.9506160449542, 28.59846724399222, 18.31938787249265};
Point(p14 + 19) = {25.9506160449542, 28.81979546051521, 18.31938787249265};
Point(p14 + 20) = {25.9506160449542, 29.0411236770382, 18.31938787249265};
Point(p14 + 21) = {25.9506160449542, 29.26245189356119, 18.31938787249265};
Point(p14 + 22) = {25.9506160449542, 29.48378011008418, 18.31938787249265};
Point(p14 + 23) = {25.9506160449542, 29.70510832660717, 18.31938787249266};
Point(p14 + 24) = {25.9506160449542, 29.92643654313016, 18.31938787249266};
Point(p14 + 25) = {25.9506160449542, 30.14776475965315, 18.31938787249266};
Point(p14 + 26) = {25.9506160449542, 30.36909297617614, 18.31938787249267};
Point(p14 + 27) = {25.9506160449542, 30.59042119269913, 18.31938787249267};
Point(p14 + 28) = {25.9506160449542, 30.81174940922212, 18.31938787249267};
Point(p14 + 29) = {25.9506160449542, 31.03307762574511, 18.31938787249267};
Point(p14 + 30) = {25.9506160449542, 31.2544058422681, 18.31938787249268};
Point(p14 + 31) = {25.9506160449542, 31.47573405879109, 18.31938787249268};
Point(p14 + 32) = {25.9506160449542, 31.69706227531408, 18.31938787249268};
Point(p14 + 33) = {25.9506160449542, 31.91839049183707, 18.31938787249268};
Point(p14 + 34) = {25.9506160449542, 32.13971870836006, 18.31938787249269};
Point(p14 + 35) = {25.9506160449542, 32.36104692488305, 18.31938787249269};
Point(p14 + 36) = {25.9506160449542, 32.58237514140604, 18.31938787249269};
Point(p14 + 37) = {25.9506160449542, 32.80370335792902, 18.31938787249269};
Point(p14 + 38) = {25.9506160449542, 33.02503157445202, 18.31938787249269};
Point(p14 + 39) = {25.9506160449542, 33.24635979097501, 18.3193878724927};
Spline(14) = {10, p14 + 1, p14 + 2, p14 + 3, p14 + 4, p14 + 5, p14 + 6, p14 + 7, p14 + 8, p14 + 9, p14 + 10, p14 + 11, p14 + 12, p14 + 13, p14 + 14, p14 + 15, p14 + 16, p14 + 17, p14 + 18, p14 + 19, p14 + 20, p14 + 21, p14 + 22, p14 + 23, p14 + 24, p14 + 25, p14 + 26, p14 + 27, p14 + 28, p14 + 29, p14 + 30, p14 + 31, p14 + 32, p14 + 33, p14 + 34, p14 + 35, p14 + 36, p14 + 37, p14 + 38, p14 + 39, 11};

/// Back triangle Top part - end ////

// Back Triange - posterior(back) part - Begin ///

p16 = newp;
Point(p16 + 1) = {25.9506160449542, 29.1521732854042, 16.0513614010009};
Point(p16 + 2) = {25.9506160449542, 29.26332546450065, 16.10871234207317};
Point(p16 + 3) = {25.9506160449542, 29.37456985276203, 16.16608812388875};
Point(p16 + 4) = {25.9506160449542, 29.48589608862277, 16.22349000807613};
Point(p16 + 5) = {25.9506160449542, 29.59729381051735, 16.28091925626375};
Point(p16 + 6) = {25.9506160449542, 29.70875265688021, 16.3383771300801};
Point(p16 + 7) = {25.9506160449542, 29.82026226614583, 16.39586489115363};
Point(p16 + 8) = {25.9506160449542, 29.93181227674864, 16.45338380111281};
Point(p16 + 9) = {25.9506160449542, 30.04339232712314, 16.5109351215861};
Point(p16 + 10) = {25.9506160449542, 30.15499205570374, 16.56852011420196};
Point(p16 + 11) = {25.9506160449542, 30.26660110092492, 16.62614004058887};
Point(p16 + 12) = {25.9506160449542, 30.37820910122115, 16.68379616237527};
Point(p16 + 13) = {25.9506160449542, 30.48980569502687, 16.74148974118965};
Point(p16 + 14) = {25.9506160449542, 30.60138052077654, 16.79922203866045};
Point(p16 + 15) = {25.9506160449542, 30.71292321690463, 16.85699431641617};
Point(p16 + 16) = {25.9506160449542, 30.82442342184558, 16.91480783608523};
Point(p16 + 17) = {25.9506160449542, 30.93587077403387, 16.97266385929613};
Point(p16 + 18) = {25.9506160449542, 31.04725491190393, 17.03056364767732};
Point(p16 + 19) = {25.9506160449542, 31.15856547389025, 17.08850846285726};
Point(p16 + 20) = {25.9506160449542, 31.26979209842726, 17.14649956646442};
Point(p16 + 21) = {25.9506160449542, 31.38092442394943, 17.20453822012728};
Point(p16 + 22) = {25.9506160449542, 31.49195208889122, 17.26262568547427};
Point(p16 + 23) = {25.9506160449542, 31.60286473168708, 17.32076322413388};
Point(p16 + 24) = {25.9506160449542, 31.71365199077148, 17.37895209773457};
Point(p16 + 25) = {25.9506160449542, 31.82430350457887, 17.4371935679048};
Point(p16 + 26) = {25.9506160449542, 31.93480891154372, 17.49548889627304};
Point(p16 + 27) = {25.9506160449542, 32.04515785010046, 17.55383934446775};
Point(p16 + 28) = {25.9506160449542, 32.15533995868358, 17.6122461741174};
Point(p16 + 29) = {25.9506160449542, 32.265344887433, 17.67071064921402};
Point(p16 + 30) = {25.9506160449542, 32.37516955290094, 17.72923550098276};
Point(p16 + 31) = {25.9506160449542, 32.48483051792439, 17.78782742762431};
Point(p16 + 32) = {25.9506160449542, 32.59434763295404, 17.84649379117393};
Point(p16 + 33) = {25.9506160449542, 32.70374074844059, 17.90524195366688};
Point(p16 + 34) = {25.9506160449542, 32.81302971483474, 17.96407927713843};
Point(p16 + 35) = {25.9506160449542, 32.92223438258716, 18.02301312362382};
Point(p16 + 36) = {25.95061604495421, 33.03137460214859, 18.08205085515833};
Point(p16 + 37) = {25.9506160449542, 33.14047022396968, 18.14119983377721};
Point(p16 + 38) = {25.9506160449542, 33.24954109850115, 18.20046742151572};
Point(p16 + 39) = {25.9506160449542, 33.35860707619369, 18.25986098040913};
Spline(16) = {12, p16 + 1, p16 + 2, p16 + 3, p16 + 4, p16 + 5, p16 + 6, p16 + 7, p16 + 8, p16 + 9, p16 + 10, p16 + 11, p16 + 12, p16 + 13, p16 + 14, p16 + 15, p16 + 16, p16 + 17, p16 + 18, p16 + 19, p16 + 20, p16 + 21, p16 + 22, p16 + 23, p16 + 24, p16 + 25, p16 + 26, p16 + 27, p16 + 28, p16 + 29, p16 + 30, p16 + 31, p16 + 32, p16 + 33, p16 + 34, p16 + 35, p16 + 36, p16 + 37, p16 + 38, p16 + 39, 11};
p17 = newp;

// Back Triange - posterior(back) part - End///




// Back Triangle - anterior(front) part - Begin ///

p18 = newp;
Point(p18 + 1) = {25.9506160449542, 28.9300740686722, 16.0513614010009};
Point(p18 + 2) = {25.9506160449542, 28.81892188957575, 16.10871234207317};
Point(p18 + 3) = {25.9506160449542, 28.70767750131438, 16.16608812388875};
Point(p18 + 4) = {25.9506160449542, 28.59635126545363, 16.22349000807613};
Point(p18 + 5) = {25.9506160449542, 28.48495354355905, 16.28091925626375};
Point(p18 + 6) = {25.9506160449542, 28.37349469719619, 16.3383771300801};
Point(p18 + 7) = {25.9506160449542, 28.26198508793057, 16.39586489115363};
Point(p18 + 8) = {25.9506160449542, 28.15043507732776, 16.45338380111281};
Point(p18 + 9) = {25.9506160449542, 28.03885502695327, 16.5109351215861};
Point(p18 + 10) = {25.9506160449542, 27.92725529837266, 16.56852011420196};
Point(p18 + 11) = {25.9506160449542, 27.81564625315148, 16.62614004058887};
Point(p18 + 12) = {25.9506160449542, 27.70403825285525, 16.68379616237527};
Point(p18 + 13) = {25.9506160449542, 27.59244165904953, 16.74148974118965};
Point(p18 + 14) = {25.9506160449542, 27.48086683329986, 16.79922203866045};
Point(p18 + 15) = {25.9506160449542, 27.36932413717178, 16.85699431641617};
Point(p18 + 16) = {25.9506160449542, 27.25782393223082, 16.91480783608523};
Point(p18 + 17) = {25.9506160449542, 27.14637658004253, 16.97266385929613};
Point(p18 + 18) = {25.9506160449542, 27.03499244217247, 17.03056364767732};
Point(p18 + 19) = {25.9506160449542, 26.92368188018616, 17.08850846285726};
Point(p18 + 20) = {25.9506160449542, 26.81245525564914, 17.14649956646442};
Point(p18 + 21) = {25.9506160449542, 26.70132293012697, 17.20453822012728};
Point(p18 + 22) = {25.9506160449542, 26.59029526518518, 17.26262568547427};
Point(p18 + 23) = {25.9506160449542, 26.47938262238931, 17.32076322413388};
Point(p18 + 24) = {25.9506160449542, 26.36859536330492, 17.37895209773457};
Point(p18 + 25) = {25.9506160449542, 26.25794384949752, 17.4371935679048};
Point(p18 + 26) = {25.9506160449542, 26.14743844253268, 17.49548889627304};
Point(p18 + 27) = {25.9506160449542, 26.03708950397593, 17.55383934446775};
Point(p18 + 28) = {25.9506160449542, 25.92690739539282, 17.6122461741174};
Point(p18 + 29) = {25.9506160449542, 25.81690246664341, 17.67071064921402};
Point(p18 + 30) = {25.9506160449542, 25.70707780117547, 17.72923550098276};
Point(p18 + 31) = {25.9506160449542, 25.59741683615201, 17.78782742762431};
Point(p18 + 32) = {25.9506160449542, 25.48789972112236, 17.84649379117393};
Point(p18 + 33) = {25.9506160449542, 25.37850660563581, 17.90524195366688};
Point(p18 + 34) = {25.9506160449542, 25.26921763924166, 17.96407927713842};
Point(p18 + 35) = {25.9506160449542, 25.16001297148924, 18.02301312362381};
Point(p18 + 36) = {25.95061604495421, 25.05087275192782, 18.08205085515831};
Point(p18 + 37) = {25.9506160449542, 24.94177713010672, 18.14119983377717};
Point(p18 + 38) = {25.9506160449542, 24.83270625557526, 18.20046742151567};
Point(p18 + 39) = {25.9506160449542, 24.72364027788271, 18.25986098040906};
Spline(18) = {12, p18 + 1, p18 + 2, p18 + 3, p18 + 4, p18 + 5, p18 + 6, p18 + 7, p18 + 8, p18 + 9, p18 + 10, p18 + 11, p18 + 12, p18 + 13, p18 + 14, p18 + 15, p18 + 16, p18 + 17, p18 + 18, p18 + 19, p18 + 20, p18 + 21, p18 + 22, p18 + 23, p18 + 24, p18 + 25, p18 + 26, p18 + 27, p18 + 28, p18 + 29, p18 + 30, p18 + 31, p18 + 32, p18 + 33, p18 + 34, p18 + 35, p18 + 36, p18 + 37, p18 + 38, p18 + 39, 10};

// Back Triangle - anterior(front) part - End ///





// Bottom Spline of ship - Curved - Begin //
p17 = newp;
Point(p17 + 1) = {119.5204361157179, 29.04112367703819, 17.6300455559415};
Point(p17 + 2) = {119.5735418764575, 29.0411236770382, 17.19782524252674};
Point(p17 + 3) = {119.7920993035465, 29.04112367703819, 16.96965345646397};
Point(p17 + 4) = {120.2028881771678, 29.0411236770382, 16.91827040721037};
Point(p17 + 5) = {120.8364847124391, 29.04112367703821, 16.94555776879794};
Point(p17 + 6) = {121.6965663790702, 29.04112367703821, 16.88835893883565};
Point(p17 + 7) = {122.6266938789654, 29.0411236770382, 16.45518759453469};
Point(p17 + 8) = {123.2371668771979, 29.0411236770382, 15.4360352268405};
Point(p17 + 9) = {123.2568976772717, 29.0411236770382, 14.08799248638371};
Point(p17 + 10) = {122.6151254409542, 29.04112367703819, 12.82006695030139};
Point(p17 + 11) = {121.3695469198386, 29.0411236770382, 11.93854105478652};
Point(p17 + 12) = {119.7340056938318, 29.0411236770382, 11.4711481375563};
Point(p17 + 13) = {118.0494437095297, 29.0411236770382, 11.2883808023031};
Point(p17 + 14) = {116.6609997678795, 29.0411236770382, 11.25646086823386};
Point(p17 + 15) = {116.1088636614788, 29.0411236770382, 11.2564145534933};
Point(p17 + 16) = {115.3627637157957, 29.0411236770382, 11.2564145534933};
Point(p17 + 17) = {113.9587482171652, 29.0411236770382, 11.2564145534933};
Point(p17 + 18) = {111.5135746077952, 29.0411236770382, 11.2564145534933};
Point(p17 + 19) = {107.8832996101724, 29.0411236770382, 11.2564145534933};
Point(p17 + 20) = {103.4946810925131, 29.0411236770382, 11.2564145534933};
Point(p17 + 21) = {98.85774075274348, 29.0411236770382, 11.2564145534933};
Point(p17 + 22) = {94.47112410721384, 29.0411236770382, 11.2564145534933};
Point(p17 + 23) = {90.54591901517306, 29.0411236770382, 11.2564145534933};
Point(p17 + 24) = {86.9904820602807, 29.0411236770382, 11.2564145534933};
Point(p17 + 25) = {83.69875640516403, 29.0411236770382, 11.2564145534933};
Point(p17 + 26) = {80.58077514101225, 29.0411236770382, 11.2564145534933};
Point(p17 + 27) = {77.67169633255861, 29.0411236770382, 11.2564145534933};
Point(p17 + 28) = {75.06528521775081, 29.0411236770382, 11.2564145534933};
Point(p17 + 29) = {72.85559809126725, 29.0411236770382, 11.2564145534933};
Point(p17 + 30) = {71.0540759406878, 29.0411236770382, 11.2564145534933};
Point(p17 + 31) = {69.41801465284388, 29.0411236770382, 11.2564145534933};
Point(p17 + 32) = {67.65599923284221, 29.0411236770382, 11.2564145534933};
Point(p17 + 33) = {65.4760536471051, 29.0411236770382, 11.25782880950083};
Point(p17 + 34) = {62.56334016157224, 29.0411236770382, 11.31970093206809};
Point(p17 + 35) = {58.57269994752332, 29.0411236770382, 11.5775073900848};
Point(p17 + 36) = {53.15688386765454, 29.0411236770382, 12.17199386391333};
Point(p17 + 37) = {46.07319617755667, 29.0411236770382, 13.21226284462873};
Point(p17 + 38) = {38.17332672593353, 29.0411236770382, 14.4761999252455};
Point(p17 + 39) = {30.96149261574722, 29.0411236770382, 15.54420267041281};
Spline(17) = {9, p17 + 1, p17 + 2, p17 + 3, p17 + 4, p17 + 5, p17 + 6, p17 + 7, p17 + 8, p17 + 9, p17 + 10, p17 + 11, p17 + 12, p17 + 13, p17 + 14, p17 + 15, p17 + 16, p17 + 17, p17 + 18, p17 + 19, p17 + 20, p17 + 21, p17 + 22, p17 + 23, p17 + 24, p17 + 25, p17 + 26, p17 + 27, p17 + 28, p17 + 29, p17 + 30, p17 + 31, p17 + 32, p17 + 33, p17 + 34, p17 + 35, p17 + 36, p17 + 37, p17 + 38, p17 + 39, 12};


//Bottom Spline of ship - Curved - End //


/////////////////////////////////////////////////////////////////////////////////////
////////////                 LINE DETAILS                ///////////////////////
////    1,  lEFT - TOP 
///     2,  LEFT - bACK 
///     3,  left Bottom
///     4,  left Front
///     5,  back - top
///     6,  back - right
///     7,  back - bottom
///     8,  right top
///     9,  right front
///     10, right botttom
///     11, front - top
///     12, front - Bottom
/////////  Ship Parameters       
///     13, Ship front
///     14, triangle - top
///     15, ship - back
///     16, triangle - back
///     17, Ship - Bottom curved spline
///     18, triangle - right
////////////////////////////////////////////////////////////////////////////////



////// SURFACES ////////////////////////

// left wall
Line Loop(1) = {1, 2, -3, -4};
Plane Surface(1) = {1};

// back Wall
Line Loop(2) = {5, 6, -7, -2};
Plane Surface(2) = {2};

// Right Wall
Line Loop(3) = {8, 9, -10, -6};
Plane Surface(3) = {3};

// Front Wall
Line Loop(4) = {11, 4, -12, -9};
Plane Surface(4) = {4};

// Top Layer with the Cavity
// Top Layer has to be subracted from the Ships Top Layer 
// To Create a annulus like surface, Create outer line Loop , create inner line loop and then 
// create a annular surface using Plane surface ( outer loop , Inner Loop)
Line Loop(5) = {11, 1, 5, 8};
Line Loop(6) = {13,14,-15};
Plane Surface(6) = {5,6};

// Bottom Surface Water
Line Loop(7) = {3, 7, 10, 12};
Plane Surface(7) = {7};

// Ship Back Surface
Line Loop(8) = {15, -16, -17};
Surface(8) = {8};

//Ship Front Surface
Line Loop(9) = {17, 18, -13};
Surface(9) = {9};

// Triangular left surface
Line Loop(10) = {18, 14, -16};
Surface(10) = {10};

///// Surface Details 

// 1 --> Left Wall
// 2 --> back Wall
// 3 --> Right Wall
// 4 --> Front Wall
// 6 --> Top surface with cavity
// 7 --> Bottom Wall
// 8 --> -SHIP- Back
// 9 --> -SHIP- Front
//10 --> -SHIP- Triangular


////////////////   End of Surfaces /////////////////////



//////////////  VOLUME /////////////////////

Surface Loop(1) = {1, 2, 3, 4, 6,9,10, 8, 7};
Volume(1) = {1};

/////////  End of Volume ////////////////////

/////////// PHYSICAL TAGS //////////////////

Physical Surface(1001) = {1};        // Left Wall - water
Physical Surface(1002) = {2};        // back Wall - water
Physical Surface(1003) = {3};        // Right Wall - water
Physical Surface(1004) = {4};        // Front Wall- water
Physical Surface(1006) = {6};       // Top surface with cavity 
Physical Surface(1007) = {7};        // Bottom Wall -water
Physical Surface(1008) = {8,9,10};   // Fluid Structure interface


