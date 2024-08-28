//==================================================================================
// BSD 2-Clause License
//
// Copyright (c) 2014-2022, NJIT, Duality Technologies Inc. and other contributors
//
// All rights reserved.
//
// Author TPOC: contact@openfhe.org
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//==================================================================================

/*

Example for CKKS bootstrapping

*/

#define PROFILE

#include "openfhe.h"
#include "math/dftransform.h"
#include "math/chebyshev.h"  // To compute coefficients
#define Pi 3.14159265358979323846

using namespace std::literals;
using namespace lbcrypto;

using DggType = typename DCRTPoly::DggType;
using DugType = typename DCRTPoly::DugType;

uint32_t m_dim1;
uint32_t m_gs;
uint32_t m_levelEnc;
uint32_t m_levelDec;
std::vector<int32_t> m_paramsEnc;
std::vector<int32_t> m_paramsDec;
const uint32_t K_SPARSE  = 12;   // 28;   // upper bound for the number of overflows in the sparse secret case
const uint32_t K_UNIFORM = 512;  // upper bound for the number of overflows in the uniform secret case
static const uint32_t R_UNIFORM =
    6;  // number of double-angle iterations in CKKS bootstrapping. Must be static because it is used in a static function.
static const uint32_t R_SPARSE = 0;
// 1; // 3;  // number of double-angle iterations in CKKS bootstrapping. Must be static because it is used in a static function.

std::vector<ConstPlaintext> m_U0Pre;
std::vector<ConstPlaintext> m_U0hatTPre;
std::vector<std::vector<ConstPlaintext>> m_U0hatTPreFFT;
std::vector<std::vector<ConstPlaintext>> m_U0PreFFT;

lbcrypto::EvalKey<lbcrypto::DCRTPoly> m_conjKey;

/*// Coefficients for the function std::cos(2*Pi/2.0 * (x - 0.25)) in [-12, 12] of degree 58
static constexpr std::initializer_list<double> coeff_sin_12_mod4{
    0.129512345213235,      -0.1286677896351604,    0.1363383834589498,    -0.1142018378613941,   0.1545141704647699,
    -0.08141291251844558,   0.1761096162294667,     -0.02535548062057226,  0.1855256664050149,    0.05338405772283146,
    0.1600366563983409,     0.1382861408618142,     0.07933727250772064,   0.1887938172259454,    -0.0508684275179663,
    0.1510126626875632,     -0.1710404861989395,    0.005828988834750617,  -0.1762975230545142,   -0.1625227446591496,
    -0.01247763459786586,   -0.1757619261547813,    0.1833360208863378,    0.03821618963152264,   0.1367050888270508,
    0.212274514692657,      -0.1448327303876925,    0.01250052501594705,   -0.1627384135154802,   -0.2292386224419867,
    0.1899446989403046,     0.07306775506252569,    0.06977739015240741,   0.1915255316898047,    -0.2655271959116844,
    -0.2874207468631774,    0.268157851210104,      0.2247230237246914,    -0.1729534361021124,   -0.1239446371601144,
    0.0834897856568161,     0.05322619063942059,    -0.03228343602012518,  -0.01870676727689579,  0.01039083154643041,
    0.005548265040359028,   -0.00285467555494661,   -0.001418216083808576, 0.000681543347003348,  0.0003173197971419828,
    -0.0001433392217935311, -6.289928883143473e-05, 2.684326388310419e-05, 1.115274392486207e-05, -4.51562572159708e-06,
    -1.782620231175819e-06, 6.887004808313971e-07,  2.546110032749553e-07, -1.073273745396853e-07};

// Coefficients for the function std::cos(2*Pi/2.0 * x) in [-12, 12] of degree 58
static constexpr std::initializer_list<double> coeff_cos_12_mod4{
    0.1831581150953023,    9.371035021412339e-16,  0.19281159095967,       1.806464582440933e-16,
    0.2185160354501053,    -1.50538715203411e-16,  0.2490566077360331,     1.994637976445197e-16,
    0.2623729135982791,    4.930142922911712e-16,  0.2263260099553762,     -2.107542012847755e-16,
    0.1121998467821088,    6.661338147750939e-16,  -0.07193882009249968,   5.946279250534737e-16,
    -0.2418877752974279,   -1.430117794432405e-16, -0.2493223481164777,    -5.494663104924504e-16,
    -0.01764604007463705,  4.102179989292951e-16,  0.2592762872089749,     9.483939057814897e-16,
    0.1933301906646326,    1.686033610278204e-15,  -0.2048244115898029,    -7.338762366166289e-16,
    -0.2301468715126729,   1.686033610278204e-15,  0.2686223693422525,     -5.268855032119387e-17,
    0.09868013150053286,   1.63334505995701e-15,   -0.3755121616372034,    -9.785016488221718e-17,
    0.3792324700381559,    7.451666402568847e-16,  -0.244593094994635,     1.76882990364008e-16,
    0.1180723871954916,    1.159148107066265e-15,  -0.04565567305966564,   -8.787697499999121e-16,
    0.01469485489729492,   -3.801102558886129e-16, -0.004037120885979591,  -3.725833201284424e-16,
    0.000963847844676895,  7.583387778371831e-16,  -0.0002027122714816887, 1.378370111081232e-15,
    3.796210784230914e-05, 9.248722315309567e-16,  -6.386059138080906e-06, 1.310157255754687e-15,
    9.739695610267261e-07, 5.344124389721092e-16,  -1.517838290204281e-07};
*/

// Coefficients for the function std::sin(2*Pi*x) in [-12, 12] of degree 118
static constexpr std::initializer_list<double> coeff_sin_12_mod4{
    1.646383776275803e-16,  -0.1293058604986413,   -3.312592963705627e-16, -0.122241376393473,
    6.501568657862655e-18,  -0.107080271437615,    -6.157539464307591e-17, -0.08207830153574547,
    3.658225683162855e-16,  -0.04550823640876773,  2.814916833706903e-17,  0.003086187736459739,
    -4.095842479371458e-16, 0.06111285853253125,   -4.38010388930043e-16,  0.1209846067978416,
    1.437852521453549e-15,  0.1689789560194281,    3.502829445735017e-16,  0.186590159745329,
    -5.171079484313048e-16, 0.1562685554148288,    -9.249428953842724e-16, 0.07211623835734007,
    -7.71500044055433e-16,  -0.04769608275750077,  3.246702626634911e-16,  -0.1556788381731591,
    -4.995712060649174e-17, -0.1891567475932377,   7.711501838585974e-17,  -0.1092338568321133,
    1.597694898882841e-17,  0.0522614074707441,    -3.970913234084725e-17, 0.1825917885791433,
    -2.81302175764071e-16,  0.1587100214495787,    8.984118304575682e-17,  -0.02350855388312963,
    -5.526187584101241e-16, -0.1895135067432571,   4.380249664382445e-16,  -0.1341979856139577,
    4.039427522631708e-16,  0.1024023775695568,    -5.367730069951091e-16, 0.2006090574833604,
    -6.443695950302756e-16, -0.01535424569395262,  8.549708560171406e-16,  -0.213847336838837,
    -1.161827403658416e-16, -0.02124134925557792,  -7.018049773441295e-16, 0.2215475090706527,
    4.547162133289457e-16,  -0.006797591722846644, 6.29617156730372e-16,   -0.2274855512858265,
    8.479736520804275e-16,  0.1108403985739115,    8.603791115598918e-16,  0.1654877912554374,
    2.168841670217053e-16,  -0.2475897258182746,   1.883807652353402e-15,  0.07378033951098242,
    -1.032029270632387e-15, 0.1683720193493045,    -5.576334212314353e-16, -0.3064633632335527,
    -5.282524534513407e-16, 0.3074532190240889,    -1.333754535386739e-15, -0.2301917728599855,
    -1.322121683841954e-15, 0.1408463622942101,    9.441414736856291e-16,  -0.07355929357959146,
    1.44148960974982e-15,   0.03364639228060209,   -5.211459182031164e-16, -0.01371216164268495,
    3.168567182674948e-16,  0.005041043773675936,  9.815036272060371e-17,  -0.001687639359331398,
    -2.878401881924373e-16, 0.0005183571656007004, -3.225929677447794e-16, -0.0001469662620341661,
    -2.007614429508621e-16, 3.865965756989653e-05, -9.772907273358077e-16, -9.47618420422844e-06,
    3.592845558879191e-16,  2.172559835540091e-06, -1.033519820845989e-15, -4.674093240545777e-07,
    -1.638960180224197e-15, 9.463912104206526e-08, -1.699467772391514e-15, -1.808077164435976e-08,
    1.109661790559399e-15,  3.266996880143432e-09, -3.43085299899018e-16,  -5.594776938360508e-10,
    -4.628358853971733e-18, 9.09800934476106e-11,  -3.299400318783282e-16, -1.407375455889839e-11,
    2.390243042592694e-15,  2.075058411892811e-12, 4.847203695846538e-17,  -2.916671485330324e-13,
    3.168904287552107e-16,  3.866324532648082e-14, -1.304953479104785e-16};

// Coefficients for the function std::cos(2*Pi*x) in [-12, 12] of degree 118
static constexpr std::initializer_list<double> coeff_cos_12_mod4{
    0.1297324437190043,    3.078769732153795e-16,  0.1331623881971898,     -2.630948680204152e-16,
    0.1428900478488792,    -1.679328944811161e-16, 0.1570920099379769,     7.034522357708975e-16,
    0.1723323719163123,    -2.574971048710447e-16, 0.1831966630800875,     -5.653740780864242e-16,
    0.1822961626852133,    -2.089831575765001e-16, 0.1612222667210367,     8.95642103899286e-17,
    0.1130840211962715,    -3.955752625555179e-16, 0.03688481883077169,    1.287485524355224e-15,
    -0.05715487764584995,  5.187260518416697e-16,  -0.1442030733001987,    -1.492736839832143e-16,
    -0.1882007539119627,   3.414635521116027e-16,  -0.1565713066745595,    -2.761563153689465e-16,
    -0.04507455183928333,  -7.314410515177502e-16, 0.1004340665875224,     1.97787631277759e-16,
    0.1902571259712754,    -6.082902622315984e-16, 0.1445099886457774,     -4.179663151530001e-16,
    -0.02500892805632881,  -6.269494727295001e-16, -0.1807757487198436,    -8.396644724055806e-17,
    -0.1564559821094688,   -1.703585918458433e-15, 0.04965109565979883,    9.852063142892146e-16,
    0.2027187171310926,    1.0449157878825e-15,    0.08048488282190923,    -2.313742101739822e-16,
    -0.16961718167671,     -8.266030250570492e-16, -0.1496602648629826,    5.299215781404108e-16,
    0.1396360512989934,    3.955752625555179e-16,  0.1694985985980304,     1.640144602765567e-15,
    -0.1537215623863027,   -1.156871050869911e-16, -0.1434437943100919,    -2.724244732693661e-16,
    0.2125764637196288,    -2.780222364187367e-16, 0.03322835755860232,    -1.156871050869911e-16,
    -0.2433227424378535,   -2.015194733773393e-16, 0.1835661520585558,     6.885248673725761e-16,
    0.05244150467131627,   -2.257764470246117e-16, -0.2557267456434503,    -2.181261707204719e-15,
    0.3214459707895934,    -1.66253565536305e-15,  -0.2739018740663259,    -1.380781576844733e-15,
    0.1840501073055895,    7.202455252190091e-16,  -0.1036269589817286,    3.937093415057278e-16,
    0.05051948926821818,   -3.293350652879666e-16, -0.02177287098878604,   -5.709718412357947e-16,
    0.008416418909615591,  1.55990999762459e-15,   -0.002949597418465354,  -3.358657889622322e-16,
    0.0009450467013496873, -9.432230906689355e-16, -0.0002786900256678019, -1.296815129604174e-16,
    7.60644815668297e-05,  1.678395984286266e-15,  -1.93050902413091e-05,  -2.883780982450722e-15,
    4.574451087223549e-06, -8.835136170756497e-16, -1.015556045562637e-06, -1.296815129604174e-16,
    2.118872722469988e-07, 1.472211708284451e-15,  -4.166117652820384e-08, -3.731842099580358e-17,
    7.738385042465596e-09, -3.232708218761485e-16, -1.360895028321035e-09, 3.745836507453784e-16,
    2.270538332512849e-10, 9.506867748680962e-16,  -3.600744772375386e-11, 2.360390127984577e-16,
    5.434937047562565e-12, -7.169801633818763e-16, -7.825395327063854e-13, -2.145809207258706e-17,
    1.069068969671378e-13, 1.335299751256097e-16,  -1.658483108083036e-14};

// Coefficients for the function 0.5*(1-std::cos(2*Pi * x)) in [-12, 12] of degree 118
static constexpr std::initializer_list<double> coeff_cos_12_mod2{
    0.9351337781404983,     -1.01755818962945e-16,  -0.06658119409859502,  1.248819538001654e-16,
    -0.0714450239244396,    1.200019336278481e-17,  -0.07854600496898836,  -3.399500252395661e-16,
    -0.08616618595815619,   2.001245311179182e-16,  -0.09159833154004367,  2.278598918295959e-16,
    -0.09114808134260663,   -2.586572695276739e-16, -0.08061113336051826,  1.510484386633538e-16,
    -0.05654201059813585,   2.348744919723855e-16,  -0.01844240941538587,  -6.330514980815201e-16,
    0.02857743882292518,    -3.824355180437951e-16, 0.07210153665009927,   3.942896800793037e-16,
    0.0941003769559814,     -1.760433701916306e-16, 0.07828565333727965,   2.578714791978832e-16,
    0.02253727591964178,    5.135742693338384e-16,  -0.05021703329376143,  4.671789578601706e-16,
    -0.09512856298563767,   2.101202886313031e-16,  -0.07225499432288866,  2.951358039718454e-16,
    0.01250446402816436,    5.932541034839841e-16,  0.09038787435992175,   -1.569119566204444e-16,
    0.07822799105473424,    9.658273393238644e-16,  -0.02482554782989942,  -3.757936693411888e-16,
    -0.1013593585655461,    -2.450743915439825e-17, -0.04024244141095451,  1.245121064729441e-16,
    0.08480859083835478,    2.657160862138949e-16,  0.07483013243149089,   -2.705082427332634e-16,
    -0.06981802564949659,   -1.45224838023643e-16,  -0.08474929929901531,  -6.248281320975312e-16,
    0.07686078119315128,    4.781510098193152e-16,  0.07172189715504591,   4.27767781814887e-16,
    -0.1062882318598141,    5.369070375026453e-16,  -0.01661417877930132,  -7.895753427819478e-16,
    0.121661371218927,      1.248872922431102e-16,  -0.09178307602927802,  -4.124464505630812e-16,
    -0.02622075233565813,   1.105150222601214e-16,  0.1278633728217251,    1.068217877799184e-15,
    -0.1607229853947966,    9.820451873300615e-16,  0.1369509370331628,    5.450714670765227e-16,
    -0.09202505365279466,   1.046270044088757e-16,  0.05181347949086465,   4.46543384604835e-17,
    -0.02525974463410923,   -7.860605261828866e-17, 0.01088643549439487,   2.0254893871483e-16,
    -0.004208209454807728,  -3.5447463303043e-16,   0.00147479870923288,   -9.250280257544333e-17,
    -0.0004725233506745082, 2.259176096919221e-16,  0.0001393450128335047, 1.28673913885229e-16,
    -3.803224078333074e-05, 1.946196284040929e-17,  9.65254512078861e-06,  1.43653279429626e-15,
    -2.287225544178077e-06, 3.173508018389293e-16,  5.077780214504138e-07, -2.027571522255282e-16,
    -1.059436359867743e-07, 1.859362665863535e-16,  2.083058804010859e-08, -3.733388183836435e-16,
    -3.869192393106865e-09, 1.198285972325974e-15,  6.804476625517498e-10, 2.225686336024247e-16,
    -1.135266175226565e-10, -7.636819513060476e-16, 1.800467190561717e-11, 1.128496799542438e-16,
    -2.717162903015563e-12, 7.514197044559846e-16,  3.908569650776809e-13, 5.911306274748373e-17,
    -5.37864587260538e-14,  -7.71191681355089e-17,  8.018794618871971e-15};

// Coefficients for standard bootstrapping K = 12
static const inline std::vector<double> g_coefficientsSparse{
    0.051667950339505692,    -0.051331021411262792,  0.054391145603358268,     -0.045559941622459216,
    0.061642235519550802,    -0.032479052974230690,  0.070257571899204785,     -0.010115373259445478,
    0.074014032428612139,    0.021297157725027743,   0.063845388651374568,     0.055168188383325954,
    0.031650992415060121,    0.075317835969809566,   -0.020293566474452410,    0.060245336022068822,
    -0.068235281605174836,   0.0023254300981702058,  -0.070332535876492297,    -0.064837194371420742,
    -0.0049778560004883976,  -0.070118863627936762,  0.073140490252120063,     0.015246053839852879,
    0.054537439879144070,    0.084685278962595631,   -0.057779899737632208,    0.0049869879560770916,
    -0.064923233796776794,   -0.091452978793089529,  0.075776971345408659,     0.029149816828457239,
    0.027837151147861859,    0.076407632367427245,   -0.10593002504560581,     -0.11466428818827884,
    0.10697950466930695,     0.089651515543433494,   -0.068998438201839998,    -0.049446756192184235,
    0.033307605480156918,    0.021234177870771417,   -0.012879227585062743,    -0.0074629203963833763,
    0.0041453420323998430,   0.0022134375074726550,  -0.0011388507756961365,   -0.00056578635857651818,
    0.00027189645704615731,  0.00012659228348801623, -0.000057184076013433917, -0.000025093185722304338,
    0.000010708912907131761, 4.4493010938744704e-6,  -1.8014740230072435e-6,   -7.1116258059986816e-7,
    2.7475174030438802e-7,   1.0157509384824620e-7,  -4.2817427341919936e-8};

// Chebyshev series coefficients for the OPTIMIZED/uniform case
static const inline std::vector<double> g_coefficientsUniform{
    0.15421426400235561,    -0.0037671538417132409,  0.16032011744533031,      -0.0034539657223742453,
    0.17711481926851286,    -0.0027619720033372291,  0.19949802549604084,      -0.0015928034845171929,
    0.21756948616367638,    0.00010729951647566607,  0.21600427371240055,      0.0022171399198851363,
    0.17647500259573556,    0.0042856217194480991,   0.086174491919472254,     0.0054640252312780444,
    -0.046667988130649173,  0.0047346914623733714,   -0.17712686172280406,     0.0016205080004247200,
    -0.22703114241338604,   -0.0028145845916205865,  -0.13123089730288540,     -0.0056345646688793190,
    0.078818395388692147,   -0.0037868875028868542,  0.23226434602675575,      0.0021116338645426574,
    0.13985510526186795,    0.0059365649669377071,   -0.13918475289368595,     0.0018580676740836374,
    -0.23254376365752788,   -0.0054103844866927788,  0.056840618403875359,     -0.0035227192748552472,
    0.25667909012207590,    0.0055029673963982112,   -0.073334392714092062,    0.0027810273357488265,
    -0.24912792167850559,   -0.0069524866497120566,  0.21288810409948347,      0.0017810057298691725,
    0.088760951809475269,   0.0055957188940032095,   -0.31937177676259115,     -0.0087539416335935556,
    0.34748800245527145,    0.0075378299617709235,   -0.25116537379803394,     -0.0047285674679876204,
    0.13970502851683486,    0.0023672533925155220,   -0.063649401080083698,    -0.00098993213448982727,
    0.024597838934816905,   0.00035553235917057483,  -0.0082485030307578155,   -0.00011176184313622549,
    0.0024390574829093264,  0.000031180384864488629, -0.00064373524734389861,  -7.8036008952377965e-6,
    0.00015310015145922058, 1.7670804180220134e-6,   -0.000033066844379476900, -3.6460909134279425e-7,
    6.5276969021754105e-6,  6.8957843666189918e-8,   -1.1842811187642386e-6,   -1.2015133285307312e-8,
    1.9839339947648331e-7,  1.9372045971100854e-9,   -3.0815418032523593e-8,   -2.9013806338735810e-10,
    4.4540904298173700e-9,  4.0505136697916078e-11,  -6.0104912807134771e-10,  -5.2873323696828491e-12,
    7.5943206779351725e-11, 6.4679566322060472e-13,  -9.0081200925539902e-12,  -7.4396949275292252e-14,
    1.0057423059167244e-12, 8.1701187638005194e-15,  -1.0611736208855373e-13,  -8.9597492970451533e-16,
    1.1421575296031385e-14};

std::vector<std::complex<double>> DecryptWithoutDecode(const CryptoContextImpl<DCRTPoly>& cc,
                                                       ConstCiphertext<DCRTPoly> cTemp,
                                                       const PrivateKey<DCRTPoly> privateKey, uint32_t slots,
                                                       uint32_t ringDim, bool scale);

std::vector<Poly> EncryptBFVCoeff(std::vector<int64_t> input, BigInteger Q, BigInteger p,
                                  const PrivateKey<DCRTPoly> privateKey);

std::vector<int64_t> DecryptBFVCoeff(const std::vector<Poly>& input, BigInteger Q, BigInteger p,
                                     const PrivateKey<DCRTPoly> privateKey, uint32_t numSlots);

std::vector<Poly> EncryptBFVCoeff(std::vector<int64_t> input, BigInteger Q, BigInteger p,
                                  const PrivateKey<DCRTPoly> privateKey,
                                  std::shared_ptr<lbcrypto::M4DCRTParams>& elementParams);

std::vector<int64_t> DecryptBFVCoeff(const std::vector<Poly>& input, BigInteger Q, BigInteger p,
                                     const PrivateKey<DCRTPoly> privateKey,
                                     std::shared_ptr<lbcrypto::M4DCRTParams>& elementParams, uint32_t numSlots);

std::vector<double> DecryptCKKSCoeff(const std::vector<Poly>& input, BigInteger Q,
                                     const PrivateKey<DCRTPoly> privateKey, uint32_t numSlots);

std::vector<double> DecryptCKKSCoeff(const std::vector<DCRTPoly>& input, const PrivateKey<DCRTPoly> privateKey,
                                     uint32_t numSlots);

std::vector<DCRTPoly> ModSwitchUp(const std::vector<Poly>& input, BigInteger Qfrom, BigInteger Qto,
                                  const std::shared_ptr<lbcrypto::M4DCRTParams>& elementParams);
std::vector<DCRTPoly> ModSwitchDown(const std::vector<Poly>& input, BigInteger Qfrom, BigInteger Qto,
                                    const std::shared_ptr<lbcrypto::M4DCRTParams>& elementParams);

void FitToNativeVector(uint32_t ringDim, const std::vector<int64_t>& vec, int64_t bigBound, NativeVector* nativeVec);

Plaintext MakeAuxPlaintext(const CryptoContextImpl<DCRTPoly>& cc, const std::shared_ptr<DCRTPoly::Params> params,
                           const std::vector<std::complex<double>>& value, size_t noiseScaleDeg, uint32_t level,
                           usint slots);

std::vector<ConstPlaintext> EvalLinearTransformPrecompute(const CryptoContextImpl<DCRTPoly>& cc,
                                                          const std::vector<std::vector<std::complex<double>>>& A,
                                                          double scale, uint32_t L);
std::vector<ConstPlaintext> EvalLinearTransformPrecompute(const CryptoContextImpl<DCRTPoly>& cc,
                                                          const std::vector<std::vector<std::complex<double>>>& A,
                                                          const std::vector<std::vector<std::complex<double>>>& B,
                                                          uint32_t orientation, double scale, uint32_t L);
std::vector<std::vector<ConstPlaintext>> EvalCoeffsToSlotsPrecompute(const CryptoContextImpl<DCRTPoly>& cc,
                                                                     const std::vector<std::complex<double>>& A,
                                                                     const std::vector<uint32_t>& rotGroup, bool flag_i,
                                                                     double scale, uint32_t L);
std::vector<std::vector<ConstPlaintext>> EvalSlotsToCoeffsPrecompute(const CryptoContextImpl<DCRTPoly>& cc,
                                                                     const std::vector<std::complex<double>>& A,
                                                                     const std::vector<uint32_t>& rotGroup, bool flag_i,
                                                                     double scale, uint32_t L);

void EvalFuncBTSetup(const CryptoContextImpl<DCRTPoly>& cc, uint32_t numSlots, uint32_t digitSize,
                     std::vector<uint32_t> dim1, std::vector<uint32_t> levelBudget, double scaleMod);

void EvalFuncBTKeyGen(const PrivateKey<DCRTPoly> privateKey, uint32_t slots);

EvalKey<DCRTPoly> ConjugateKeyGen(const PrivateKey<DCRTPoly> privateKey);
Ciphertext<DCRTPoly> Conjugate(ConstCiphertext<DCRTPoly> ciphertext);

std::vector<int32_t> FindBootstrapRotationIndices(uint32_t slots, uint32_t M);
std::vector<int32_t> FindLinearTransformRotationIndices(uint32_t slots, uint32_t M);
std::vector<int32_t> FindCoeffsToSlotsRotationIndices(uint32_t slots, uint32_t M);
std::vector<int32_t> FindSlotsToCoeffsRotationIndices(uint32_t slots, uint32_t M);

Ciphertext<DCRTPoly> EvalLinearTransform(const std::vector<ConstPlaintext>& A, ConstCiphertext<DCRTPoly> ct);
Ciphertext<DCRTPoly> EvalCoeffsToSlots(const std::vector<std::vector<ConstPlaintext>>& A,
                                       ConstCiphertext<DCRTPoly> ctxt);
Ciphertext<DCRTPoly> EvalSlotsToCoeffs(const std::vector<std::vector<ConstPlaintext>>& A,
                                       ConstCiphertext<DCRTPoly> ctxt);

Ciphertext<DCRTPoly> EvalMultExt(ConstCiphertext<DCRTPoly> ciphertext, ConstPlaintext plaintext);

void EvalAddExtInPlace(Ciphertext<DCRTPoly>& ciphertext1, ConstCiphertext<DCRTPoly> ciphertext2);

Ciphertext<DCRTPoly> EvalAddExt(ConstCiphertext<DCRTPoly> ciphertext1, ConstCiphertext<DCRTPoly> ciphertext2);

void AdjustCiphertext(Ciphertext<DCRTPoly>& ciphertext, double correction);

void ApplyDoubleAngleIterations(Ciphertext<DCRTPoly>& ciphertext, uint32_t numIter);

Ciphertext<DCRTPoly> EvalFuncBT(ConstCiphertext<DCRTPoly> ciphertext, uint32_t digitBitSize, BigInteger initialScaling,
                                uint64_t postScaling, bool step = false);

void SimpleBootstrapExample();
void TestModApprox();
void Floor();
void Sign();

std::vector<DCRTPoly> ModSwitchUp(const std::vector<Poly>& input, BigInteger Qfrom, BigInteger Qto,
                                  const std::shared_ptr<lbcrypto::M4DCRTParams>& elementParams) {
    Poly bPoly = input[0];
    bPoly.SwitchModulus(Qto, 1, 0, 0);  // need to switch to modulus before because the new modulus is bigger
    bPoly = bPoly.MultiplyAndRound(Qto, Qfrom);

    Poly aPoly = input[1];
    aPoly.SwitchModulus(Qto, 1, 0, 0);  // need to switch to modulus before because the new modulus is bigger
    aPoly = aPoly.MultiplyAndRound(Qto, Qfrom);

    // Going to Double-CRT
    std::vector<DCRTPoly> output(2);
    output[0] = DCRTPoly(bPoly, elementParams);
    output[1] = DCRTPoly(aPoly, elementParams);

    // Switching to NTT representation
    output[0].SetFormat(Format::EVALUATION);
    output[1].SetFormat(Format::EVALUATION);

    return output;
}

std::vector<DCRTPoly> ModSwitchDown(const std::vector<Poly>& input, BigInteger Qfrom, BigInteger Qto,
                                    const std::shared_ptr<lbcrypto::M4DCRTParams>& elementParams) {
    Poly bPoly = input[0];
    bPoly      = bPoly.MultiplyAndRound(Qto, Qfrom);
    bPoly.SwitchModulus(Qto, 1, 0, 0);

    Poly aPoly = input[1];
    aPoly      = aPoly.MultiplyAndRound(Qto, Qfrom);
    aPoly.SwitchModulus(Qto, 1, 0, 0);

    // Going to Double-CRT
    std::vector<DCRTPoly> output(2);
    output[0] = DCRTPoly(bPoly, elementParams);
    output[1] = DCRTPoly(aPoly, elementParams);

    // Switching to NTT representation
    output[0].SetFormat(Format::EVALUATION);
    output[1].SetFormat(Format::EVALUATION);

    return output;
}

std::vector<Poly> EncryptBFVCoeff(std::vector<int64_t> input, BigInteger Q, BigInteger p,
                                  const PrivateKey<DCRTPoly> privateKey) {
    // Generate encryption of 0 using the existing CKKS cryptocontext

    const auto cryptoParams =
        std::dynamic_pointer_cast<CryptoParametersRLWE<DCRTPoly>>(privateKey->GetCryptoParameters());

    const DCRTPoly& s = privateKey->GetPrivateElement();

    auto elementParams = cryptoParams->GetElementParams();

    const DggType& dgg = cryptoParams->GetDiscreteGaussianGenerator();
    DugType dug;

    DCRTPoly a(dug, elementParams, Format::EVALUATION);
    DCRTPoly e(dgg, elementParams, Format::EVALUATION);

    DCRTPoly b = e - a * s;  // encryption of 0 using Q'

    a.SetFormat(Format::COEFFICIENT);
    b.SetFormat(Format::COEFFICIENT);

    auto aPoly           = a.CRTInterpolate();
    auto bPoly           = b.CRTInterpolate();
    BigInteger bigQPrime = b.GetModulus();

    // Do modulus switching from Q' to Q
    bPoly = bPoly.MultiplyAndRound(Q, bigQPrime);
    bPoly.SwitchModulus(Q, 1, 0, 0);

    aPoly = aPoly.MultiplyAndRound(Q, bigQPrime);
    aPoly.SwitchModulus(Q, 1, 0, 0);

    auto mPoly = bPoly;
    mPoly.SetValuesToZero();

    BigInteger delta = Q / p;
    uint32_t gap     = mPoly.GetLength() / (2.0 * input.size());

    for (size_t i = 0; i < input.size() && i < mPoly.GetLength(); i++) {
        BigInteger entry{input[i]};

        if (input[i] < 0) {
            entry = mPoly.GetModulus() - BigInteger(static_cast<uint64_t>(llabs(input[i])));
        }
        mPoly[i * gap] = delta * entry;
    }

    bPoly += mPoly;  // Adds the message

    return {bPoly, aPoly};
}

std::vector<Poly> EncryptBFVCoeff(std::vector<int64_t> input, BigInteger Q, BigInteger p,
                                  const PrivateKey<DCRTPoly> privateKey,
                                  std::shared_ptr<lbcrypto::M4DCRTParams>& elementParams) {
    // Generate encryption of 0 using the existing CKKS cryptocontext but using fewer limbs

    const auto cryptoParams =
        std::dynamic_pointer_cast<CryptoParametersRLWE<DCRTPoly>>(privateKey->GetCryptoParameters());

    const DCRTPoly& s = privateKey->GetPrivateElement();

    // auto elementParams = cryptoParams->GetElementParams();

    const DggType& dgg = cryptoParams->GetDiscreteGaussianGenerator();
    DugType dug;

    DCRTPoly a(dug, elementParams, Format::EVALUATION);
    DCRTPoly e(dgg, elementParams, Format::EVALUATION);

    size_t sizeQ  = s.GetParams()->GetParams().size();
    size_t sizeQl = elementParams->GetParams().size();
    size_t diffQl = sizeQ - sizeQl;

    auto scopy(s);
    scopy.DropLastElements(diffQl);

    DCRTPoly b = e - a * scopy;  // encryption of 0 using Q'

    a.SetFormat(Format::COEFFICIENT);
    b.SetFormat(Format::COEFFICIENT);

    auto aPoly           = a.CRTInterpolate();
    auto bPoly           = b.CRTInterpolate();
    BigInteger bigQPrime = b.GetModulus();

    // Do modulus switching from Q' to Q
    if (Q < bigQPrime) {
        bPoly = bPoly.MultiplyAndRound(Q, bigQPrime);
        bPoly.SwitchModulus(Q, 1, 0, 0);

        aPoly = aPoly.MultiplyAndRound(Q, bigQPrime);
        aPoly.SwitchModulus(Q, 1, 0, 0);
    }
    else {
        bPoly.SwitchModulus(Q, 1, 0, 0);
        bPoly = bPoly.MultiplyAndRound(Q, bigQPrime);

        aPoly.SwitchModulus(Q, 1, 0, 0);
        aPoly = aPoly.MultiplyAndRound(Q, bigQPrime);
    }

    auto mPoly = bPoly;
    mPoly.SetValuesToZero();

    BigInteger delta = Q / p;
    uint32_t gap     = mPoly.GetLength() / (2.0 * input.size());

    for (size_t i = 0; i < input.size() && i < mPoly.GetLength(); i++) {
        BigInteger entry{input[i]};

        if (input[i] < 0) {
            entry = mPoly.GetModulus() - BigInteger(static_cast<uint64_t>(llabs(input[i])));
        }
        mPoly[i * gap] = delta * entry;
    }

    bPoly += mPoly;  // Adds the message

    return {bPoly, aPoly};
}

std::vector<int64_t> DecryptBFVCoeff(const std::vector<Poly>& input, BigInteger Q, BigInteger p,
                                     const PrivateKey<DCRTPoly> privateKey, uint32_t numSlots) {
    const DCRTPoly& s = privateKey->GetPrivateElement();

    BigInteger bigQPrime = s.GetModulus();

    // Poly bPoly = input[0];
    // bPoly.SwitchModulus(bigQPrime, 1, 0, 0); //need to switch to modulus before because the new modulus is bigger
    // bPoly = bPoly.MultiplyAndRound(bigQPrime, Q);

    // Poly aPoly = input[1];
    // aPoly.SwitchModulus(bigQPrime, 1, 0, 0); //need to switch to modulus before because the new modulus is bigger
    // aPoly = aPoly.MultiplyAndRound(bigQPrime, Q);

    // // Going back to Double-CRT
    // DCRTPoly b = DCRTPoly(bPoly, s.GetParams());
    // DCRTPoly a = DCRTPoly(aPoly, s.GetParams());

    // // Switching to NTT representation
    // b.SetFormat(Format::EVALUATION);
    // a.SetFormat(Format::EVALUATION);

    auto ba = ModSwitchUp(input, Q, bigQPrime, s.GetParams());

    auto m = ba[0] + ba[1] * s;

    m.SetFormat(Format::COEFFICIENT);

    auto mPoly   = m.CRTInterpolate();
    uint32_t gap = mPoly.GetLength() / (2 * numSlots);

    std::cerr << "\nmPolys BFV:" << std::endl;
    for (size_t i = 0; i < numSlots; i++) {
        std::cerr << mPoly[i] << " ";
    }
    std::cerr << "\n";

    mPoly = mPoly.MultiplyAndRound(Q, bigQPrime);
    mPoly.SwitchModulus(Q, 1, 0, 0);

    mPoly = mPoly.MultiplyAndRound(p, Q);
    mPoly.SwitchModulus(p, 1, 0, 0);

    BigInteger half = p >> 1;

    std::vector<int64_t> output(numSlots);

    for (size_t i = 0, idx = 0; i < numSlots; ++i, idx += gap) {
        int64_t val;
        if (mPoly[idx] > half) {
            val = (-(p - mPoly[idx]).ConvertToInt());
        }
        else {
            val = mPoly[idx].ConvertToInt();
        }
        output[i] = val;
    }

    return output;
}

std::vector<int64_t> DecryptBFVCoeff(const std::vector<Poly>& input, BigInteger Q, BigInteger p,
                                     const PrivateKey<DCRTPoly> privateKey,
                                     std::shared_ptr<lbcrypto::M4DCRTParams>& elementParams, uint32_t numSlots) {
    const DCRTPoly& s = privateKey->GetPrivateElement();

    auto bigQPrime = elementParams->GetModulus();

    std::vector<lbcrypto::DCRTPoly> ba(2);
    if (Q < bigQPrime) {
        ba = ModSwitchUp(input, Q, bigQPrime, elementParams);
    }
    else {
        ba = ModSwitchDown(input, Q, bigQPrime, elementParams);
    }

    size_t sizeQ  = s.GetParams()->GetParams().size();
    size_t sizeQl = elementParams->GetParams().size();
    size_t diffQl = sizeQ - sizeQl;

    auto scopy(s);
    scopy.DropLastElements(diffQl);

    auto m = ba[0] + ba[1] * scopy;

    m.SetFormat(Format::COEFFICIENT);

    auto mPoly   = m.CRTInterpolate();
    uint32_t gap = mPoly.GetLength() / (2 * numSlots);

    std::cerr << "\nmPolys BFV:" << std::endl;
    for (size_t i = 0; i < mPoly.GetLength() / 2; i++) {
        std::cerr << mPoly[i] << " ";
    }
    std::cerr << "\n";

    if (Q < bigQPrime) {
        mPoly = mPoly.MultiplyAndRound(Q, bigQPrime);
        mPoly.SwitchModulus(Q, 1, 0, 0);
    }
    else {
        mPoly.SwitchModulus(Q, 1, 0, 0);
        mPoly = mPoly.MultiplyAndRound(Q, bigQPrime);
    }

    mPoly = mPoly.MultiplyAndRound(p, Q);
    mPoly.SwitchModulus(p, 1, 0, 0);

    BigInteger half = p >> 1;

    std::vector<int64_t> output(numSlots);
    for (size_t i = 0, idx = 0; i < numSlots; ++i, idx += gap) {
        int64_t val;
        if (mPoly[idx] > half) {
            val = (-(p - mPoly[idx]).ConvertToInt());
        }
        else {
            val = mPoly[idx].ConvertToInt();
        }
        output[i] = val;
    }

    return output;
}

std::vector<double> DecryptCKKSCoeff(const std::vector<Poly>& input, BigInteger Q,
                                     const PrivateKey<DCRTPoly> privateKey, uint32_t numSlots) {
    const DCRTPoly& s = privateKey->GetPrivateElement();

    BigInteger bigQPrime = s.GetModulus();

    auto ba = ModSwitchUp(input, Q, bigQPrime, s.GetParams());

    auto m = ba[0] + ba[1] * s;

    m.SetFormat(Format::COEFFICIENT);

    auto mPoly = m.CRTInterpolate();

    std::cerr << "\nmPolys CKKS:" << std::endl;
    for (size_t i = 0; i < mPoly.GetLength() / 2; i++) {
        std::cerr << mPoly[i] << " ";
    }
    std::cerr << "\n";

    BigInteger half = bigQPrime >> 1;

    std::vector<double> output(numSlots);
    uint32_t gap = mPoly.GetLength() / (2 * numSlots);

    // for (size_t i = 0; i < output.size(); i++) {
    //     if (mPoly[i] > half) {
    //         output[i] = (-(bigQPrime - mPoly[i]).ConvertToDouble());
    //     }
    //     else {
    //         output[i] = mPoly[i].ConvertToDouble();
    //     }
    // }

    for (size_t i = 0, idx = 0; i < numSlots; ++i, idx += gap) {
        if (mPoly[idx] > half) {
            output[i] = (-(bigQPrime - mPoly[idx]).ConvertToDouble());
        }
        else {
            output[i] = mPoly[idx].ConvertToDouble();
        }
    }

    return output;
}

std::vector<double> DecryptCKKSCoeff(const std::vector<DCRTPoly>& input, const PrivateKey<DCRTPoly> privateKey,
                                     uint32_t numSlots) {
    const DCRTPoly& s = privateKey->GetPrivateElement();

    // std::cout << "s modulus = " << s.GetModulus() << ", input modulus = " << input[0].GetModulus() << std::endl;
    // std::cout << "s towers = " << s.GetNumOfElements() << ", input towers " << input[0].GetNumOfElements() << std::endl;
    // BigInteger bigQPrime = s.GetModulus();
    BigInteger QPrime = input[0].GetModulus();

    size_t sizeQ  = s.GetParams()->GetParams().size();
    size_t sizeQl = input[0].GetParams()->GetParams().size();
    size_t diffQl = sizeQ - sizeQl;

    auto scopy(s);
    scopy.DropLastElements(diffQl);

    auto m = input[0] + input[1] * scopy;
    m.SetFormat(Format::COEFFICIENT);
    auto mPoly = m.CRTInterpolate();

    std::cerr << "\nmPolys CKKS:" << std::endl;
    for (size_t i = 0; i < mPoly.GetLength() / 2; i++) {
        std::cerr << mPoly[i] << " ";
    }
    std::cerr << "\n";

    BigInteger half = QPrime >> 1;

    std::vector<double> output(2 * numSlots);

    // for (size_t i = 0; i < output.size(); i++) {
    //     if (mPoly[i] > half) {
    //         output[i] = (-(QPrime - mPoly[i]).ConvertToDouble());
    //     }
    //     else {
    //         output[i] = mPoly[i].ConvertToDouble();
    //     }
    // }

    uint32_t gap = mPoly.GetLength() / (2 * numSlots);

    for (size_t i = 0, idx = 0; i < 2 * numSlots; ++i, idx += gap) {
        if (mPoly[idx] > half) {
            output[i] = (-(QPrime - mPoly[idx]).ConvertToDouble());
        }
        else {
            output[i] = mPoly[idx].ConvertToDouble();
        }
    }

    return output;
}

std::vector<std::complex<double>> DecryptWithoutDecode(const CryptoContextImpl<DCRTPoly>& cc,
                                                       ConstCiphertext<DCRTPoly> cTemp,
                                                       const PrivateKey<DCRTPoly> privateKey, uint32_t slots,
                                                       uint32_t ringDim, bool scale) {
    Plaintext decrypted = cc.GetPlaintextForDecrypt(cTemp->GetEncodingType(), cTemp->GetElements()[0].GetParams(),
                                                    cc.GetEncodingParams());
    bool isNativePoly   = true;
    DecryptResult result;

    if ((cTemp->GetEncodingType() == CKKS_PACKED_ENCODING) &&
        (cTemp->GetElements()[0].GetParams()->GetParams().size() > 1)) {
        result       = cc.GetScheme()->Decrypt(cTemp, privateKey, &decrypted->GetElement<Poly>());
        isNativePoly = false;
    }
    else {
        result       = cc.GetScheme()->Decrypt(cTemp, privateKey, &decrypted->GetElement<NativePoly>());
        isNativePoly = true;
    }

    auto elemModulus   = decrypted->GetElementModulus();
    auto noiseScaleDeg = cTemp->GetNoiseScaleDeg();
    auto scalingFactor = cTemp->GetScalingFactor();

    decrypted->SetScalingFactorInt(result.scalingFactorInt);

    double p     = cc.GetEncodingParams()->GetPlaintextModulus();
    double powP  = 0.0;
    uint32_t Nh  = ringDim / 2;
    uint32_t gap = Nh / slots;
    std::vector<std::complex<double>> curValues(slots);

    const auto cryptoParamsCKKS = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc.GetCryptoParameters());

    auto scalTech = cryptoParamsCKKS->GetScalingTechnique();

    if (isNativePoly) {
        if (scalTech == FLEXIBLEAUTO || scalTech == FLEXIBLEAUTOEXT) {
            powP = pow(scalingFactor, -1);
        }
        else {
            powP = pow(2, -p);
        }

        const NativeInteger& q = decrypted->GetElementModulus().ConvertToInt();
        NativeInteger qHalf    = q >> 1;

        for (size_t i = 0, idx = 0; i < slots; ++i, idx += gap) {
            std::complex<double> cur;

            if (decrypted->GetElement<NativePoly>()[idx] > qHalf)
                cur.real(-((q - decrypted->GetElement<NativePoly>()[idx])).ConvertToDouble());
            else
                cur.real((decrypted->GetElement<NativePoly>()[idx]).ConvertToDouble());

            if (decrypted->GetElement<NativePoly>()[idx + Nh] > qHalf)
                cur.imag(-((q - decrypted->GetElement<NativePoly>()[idx + Nh])).ConvertToDouble());
            else
                cur.imag((decrypted->GetElement<NativePoly>()[idx + Nh]).ConvertToDouble());

            curValues[i] = (scale) ? cur * powP : cur;
        }

        // for (size_t i = 0; i < 2 * slots; ++i) {
        //     std::cout << decrypted->GetElement<NativePoly>()[i] << " ";
        // }
        // std::cout << std::endl;
    }
    else {
        powP = pow(2, -p);

        // we will bring down the scaling factor to 2^p
        double scalingFactorPre = 0.0;
        if (scalTech == FLEXIBLEAUTO || scalTech == FLEXIBLEAUTOEXT)
            scalingFactorPre = pow(scalingFactor, -1) * pow(2, p);
        else
            scalingFactorPre = pow(2, -p * (noiseScaleDeg - 1));

        const BigInteger& q = decrypted->GetElementModulus();
        BigInteger qHalf    = q >> 1;

        for (size_t i = 0, idx = 0; i < slots; ++i, idx += gap) {
            std::complex<double> cur;

            if (decrypted->GetElement<Poly>()[idx] > qHalf)
                cur.real(-((q - decrypted->GetElement<Poly>()[idx])).ConvertToDouble() * scalingFactorPre);
            else
                cur.real((decrypted->GetElement<Poly>()[idx]).ConvertToDouble() * scalingFactorPre);

            if (decrypted->GetElement<Poly>()[idx + Nh] > qHalf)
                cur.imag(-((q - decrypted->GetElement<Poly>()[idx + Nh])).ConvertToDouble() * scalingFactorPre);
            else
                cur.imag((decrypted->GetElement<Poly>()[idx + Nh]).ConvertToDouble() * scalingFactorPre);

            curValues[i] = (scale) ? cur * powP : cur;
        }
    }
    return curValues;
}

void TestModApprox() {
    CCParams<CryptoContextCKKSRNS> parameters;
    SecretKeyDist secretKeyDist = SPARSE_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetRingDim(32);

    uint32_t dcrtBits = 41;
    uint32_t firstMod = 41;
    uint32_t numSlots = 16;

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(FIXEDMANUAL);
    parameters.SetFirstModSize(firstMod);
    parameters.SetNumLargeDigits(3);
    parameters.SetBatchSize(numSlots);

    usint depth = 7 + 1 + 5;
    parameters.SetMultiplicativeDepth(depth);

    CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);

    cryptoContext->Enable(PKE);
    cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);
    cryptoContext->Enable(ADVANCEDSHE);

    usint ringDim = cryptoContext->GetRingDimension();
    std::cout << "CKKS scheme is using ring dimension " << ringDim << std::endl << std::endl;

    auto keyPair = cryptoContext->KeyGen();
    cryptoContext->EvalMultKeyGen(keyPair.secretKey);

    std::vector<double> x = {-8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7};
    std::transform(x.begin(), x.end(), x.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, 1 / 8.0));
    std::cout << "Input = " << x << std::endl << std::endl;

    Plaintext ptxt            = cryptoContext->MakeCKKSPackedPlaintext(x);
    Ciphertext<DCRTPoly> ctxt = cryptoContext->Encrypt(keyPair.publicKey, ptxt);

    double a_cheby = -12;
    double b_cheby = 12;

    auto ctxt_cos = cryptoContext->EvalChebyshevSeries(ctxt, coeff_cos_12_mod2, a_cheby, b_cheby);

    // Mod 2
    auto result = ctxt_cos;  // Andreea: for benchmarking, the multiplication by 0.5 should be in the evaluation
    std::cerr << "result levels: " << result->GetLevel() << " and depth: " << result->GetNoiseScaleDeg() << std::endl;

    Plaintext ptxt_result;
    cryptoContext->Decrypt(keyPair.secretKey, result, &ptxt_result);
    ptxt_result->SetLength(numSlots);
    std::cout << "Output after polynomial approximation for mod 2:" << ptxt_result << std::endl;
    std::vector<double> exact2 = {0,   (2 - sqrt(2)) / 4, 0.5, (2 + sqrt(2)) / 4, 1,   (2 + sqrt(2)) / 4,
                                  0.5, (2 - sqrt(2)) / 4, 0,   (2 - sqrt(2)) / 4, 0.5, (2 + sqrt(2)) / 4,
                                  1,   (2 + sqrt(2)) / 4, 0.5, (2 - sqrt(2)) / 4};

    std::vector<double> error2(ptxt_result->GetRealPackedValue());
    std::transform(exact2.begin(), exact2.end(), error2.begin(), error2.begin(), std::minus<double>());
    std::transform(error2.begin(), error2.end(), error2.begin(), [&](const double& elem) { return std::abs(elem); });
    std::cout << "Absolute error mod 2: " << error2 << std::endl;
    double mean_error = std::accumulate(error2.begin(), error2.end(), 0.0) / numSlots;
    std::cout << "Average absolute error mod 2: " << mean_error << std::endl;
    auto max_error_it = std::max_element(error2.begin(), error2.end());
    std::cout << "Max absolute error mod 2: " << *max_error_it << std::endl << std::endl;

    // Mod 4
    ctxt_cos      = cryptoContext->EvalChebyshevSeries(ctxt, coeff_cos_12_mod4, a_cheby, b_cheby);
    auto ctxt_sin = cryptoContext->EvalChebyshevSeries(ctxt, coeff_sin_12_mod4, a_cheby, b_cheby);

    // cryptoContext->EvalSquareInPlace(ctxt_sin);
    // cryptoContext->EvalAddInPlace(ctxt_sin, ctxt_sin);
    // cryptoContext->EvalSubInPlace(ctxt_sin, 1.0);
    // cryptoContext->ModReduceInPlace(ctxt_sin);

    // cryptoContext->EvalSquareInPlace(ctxt_cos);
    // cryptoContext->EvalAddInPlace(ctxt_cos, ctxt_cos);
    // cryptoContext->EvalSubInPlace(ctxt_cos, 1.0);
    // cryptoContext->ModReduceInPlace(ctxt_cos);

    result = cryptoContext->EvalAdd(cryptoContext->EvalSub(ctxt_cos, ctxt_sin), 1.0);
    cryptoContext->EvalSquareInPlace(ctxt_cos);
    cryptoContext->ModReduceInPlace(ctxt_cos);
    result = cryptoContext->EvalMult(result, ctxt_cos);
    cryptoContext->ModReduceInPlace(result);
    ctxt_sin = cryptoContext->EvalSub(2.0, ctxt_sin);
    result   = cryptoContext->EvalSub(ctxt_sin, result);

    std::cerr << "result levels: " << result->GetLevel() << " and depth: " << result->GetNoiseScaleDeg() << std::endl;

    cryptoContext->Decrypt(keyPair.secretKey, result, &ptxt_result);
    ptxt_result->SetLength(numSlots);
    std::cout << "Output after polynomial approximation for mod 4:" << ptxt_result << std::endl;
    std::vector<double> exact4 = {0, 1.5 - sqrt(2) / 2, 1, 1.5, 2, 1.5 + sqrt(2) / 2, 3, 1.5,
                                  0, 1.5 - sqrt(2) / 2, 1, 1.5, 2, 1.5 + sqrt(2) / 2, 3, 1.5};

    std::vector<double> error4(ptxt_result->GetRealPackedValue());
    std::transform(exact4.begin(), exact4.end(), error4.begin(), error4.begin(), std::minus<double>());
    std::transform(error4.begin(), error4.end(), error4.begin(), [&](const double& elem) { return std::abs(elem); });
    std::cout << "Absolute error mod 4: " << error4 << std::endl;
    mean_error = std::accumulate(error4.begin(), error4.end(), 0.0) / numSlots;
    std::cout << "Average absolute error mod 4: " << mean_error << std::endl;
    max_error_it = std::max_element(error4.begin(), error4.end());
    std::cout << "Max absolute error mod 4: " << *max_error_it << std::endl << std::endl;

    // Mod 8
    ctxt_cos = cryptoContext->EvalChebyshevSeries(ctxt, coeff_cos_12_mod4, a_cheby, b_cheby);
    ctxt_sin = cryptoContext->EvalChebyshevSeries(ctxt, coeff_sin_12_mod4, a_cheby, b_cheby);

    // cryptoContext->EvalSquareInPlace(ctxt_sin);
    // cryptoContext->EvalAddInPlace(ctxt_sin, ctxt_sin);
    // cryptoContext->EvalSubInPlace(ctxt_sin, 1.0);
    // cryptoContext->ModReduceInPlace(ctxt_sin);

    // cryptoContext->EvalSquareInPlace(ctxt_cos);
    // cryptoContext->EvalAddInPlace(ctxt_cos, ctxt_cos);
    // cryptoContext->EvalSubInPlace(ctxt_cos, 1.0);
    // cryptoContext->ModReduceInPlace(ctxt_cos);

    /*
    cryptoContext->GetScheme()->MultByIntegerInPlace(ctxt_sin, 2); // 2sin
    auto ctxt_sin_sqrt = cryptoContext->EvalMult(ctxt_sin, 2.0 + 2*sqrt(2)); // (4+4sqrt(2))sin
    cryptoContext->ModReduceInPlace(ctxt_sin_sqrt);
    auto ctxt_cos2 = cryptoContext->EvalSquare(ctxt_cos); // cos^2
    cryptoContext->ModReduceInPlace(ctxt_cos2);
    ctxt_cos = cryptoContext->EvalMult(ctxt_cos2, ctxt_cos); // cos^3
    cryptoContext->ModReduceInPlace(ctxt_cos);
    auto coeff = cryptoContext->EvalAdd(-2, cryptoContext->EvalAdd(cryptoContext->GetScheme()->MultByInteger(ctxt_sin, 7), ctxt_sin_sqrt));
    result = cryptoContext->EvalMult(ctxt_cos2, coeff);
    cryptoContext->LevelReduceInPlace(result, nullptr, 1);
    coeff = cryptoContext->EvalSub(-2, cryptoContext->GetScheme()->MultByInteger(ctxt_sin, 4));
    cryptoContext->LevelReduceInPlace(coeff, nullptr, 2);
    cryptoContext->EvalAddInPlace(result, cryptoContext->EvalMult(ctxt_cos, coeff));
    auto ctxt_cos5 = cryptoContext->EvalMult(ctxt_cos2, ctxt_cos); // cos^5
    cryptoContext->ModReduceInPlace(ctxt_cos5);
    cryptoContext->EvalSquareInPlace(ctxt_cos2); // cos^4
    cryptoContext->ModReduceInPlace(ctxt_cos2);
    coeff = cryptoContext->EvalSub(8, cryptoContext->EvalAdd(cryptoContext->GetScheme()->MultByInteger(ctxt_sin, 6), ctxt_sin_sqrt));
    cryptoContext->LevelReduceInPlace(coeff, nullptr, 1);
    cryptoContext->EvalAddInPlace(result, cryptoContext->EvalMult(ctxt_cos2, coeff));
    coeff = cryptoContext->EvalAdd(8, cryptoContext->GetScheme()->MultByInteger(ctxt_sin, 4));
    cryptoContext->LevelReduceInPlace(coeff, nullptr, 2);
    cryptoContext->EvalAddInPlace(result, cryptoContext->EvalMult(ctxt_cos5, coeff));
    ctxt_cos5 = cryptoContext->EvalSquare(ctxt_cos); // cos^6
    cryptoContext->ModReduceInPlace(ctxt_cos5);
    ctxt_cos = cryptoContext->EvalMult(ctxt_cos, ctxt_cos2); // cos^7
    cryptoContext->ModReduceInPlace(ctxt_cos);
    coeff = cryptoContext->EvalAdd(-8, cryptoContext->GetScheme()->MultByInteger(ctxt_sin_sqrt, 2));
    cryptoContext->LevelReduceInPlace(coeff, nullptr, 2);
    cryptoContext->EvalAddInPlace(result, cryptoContext->EvalMult(ctxt_cos5, coeff));
    cryptoContext->ModReduceInPlace(result);
    ctxt_sin = cryptoContext->EvalAdd(-4, ctxt_sin);
    cryptoContext->LevelReduceInPlace(ctxt_sin, nullptr, 3);
    cryptoContext->LevelReduceInPlace(ctxt_cos, nullptr, 1);
    cryptoContext->EvalSubInPlace(result, cryptoContext->EvalAdd(ctxt_sin, cryptoContext->GetScheme()->MultByInteger(ctxt_cos, 8)));
*/
    auto scheme = cryptoContext->GetScheme();

    cryptoContext->GetScheme()->MultByIntegerInPlace(ctxt_sin, 2);              // 2sin
    auto ctxt_sin_sqrt = cryptoContext->EvalMult(ctxt_sin, 2.0 + 2 * sqrt(2));  // (4+4sqrt(2))sin
    cryptoContext->ModReduceInPlace(ctxt_sin_sqrt);
    auto ctxt_cos2 = cryptoContext->EvalSquare(ctxt_cos);  // cos^2
    cryptoContext->ModReduceInPlace(ctxt_cos2);
    auto ctxt_cos3 = cryptoContext->EvalMult(ctxt_cos2, ctxt_cos);  // cos^3
    cryptoContext->ModReduceInPlace(ctxt_cos3);

    result = cryptoContext->EvalSub(scheme->MultByInteger(ctxt_sin, 7), cryptoContext->EvalAdd(2, ctxt_sin_sqrt));
    cryptoContext->LevelReduceInPlace(result, nullptr, 1);

    auto term = cryptoContext->EvalMult(ctxt_cos, cryptoContext->EvalAdd(2, scheme->MultByInteger(ctxt_sin, 4)));
    cryptoContext->ModReduceInPlace(term);
    cryptoContext->LevelReduceInPlace(term, nullptr, 1);
    cryptoContext->EvalSubInPlace(result, term);
    term = cryptoContext->EvalMult(
        ctxt_cos2,
        cryptoContext->EvalSub(8, cryptoContext->EvalAdd(ctxt_sin_sqrt, scheme->MultByInteger(ctxt_sin, 6))));
    cryptoContext->ModReduceInPlace(term);
    cryptoContext->EvalAddInPlace(result, term);

    auto result2 = cryptoContext->EvalMult(ctxt_cos, cryptoContext->EvalAdd(8, scheme->MultByInteger(ctxt_sin, 4)));
    cryptoContext->ModReduceInPlace(result2);
    cryptoContext->LevelReduceInPlace(result, nullptr, 1);
    term = cryptoContext->EvalMult(ctxt_cos2, cryptoContext->EvalAdd(-8, scheme->MultByInteger(ctxt_sin_sqrt, 2)));
    cryptoContext->ModReduceInPlace(term);
    cryptoContext->EvalAddInPlace(result2, term);
    cryptoContext->EvalSubInPlace(result2, scheme->MultByInteger(ctxt_cos3, 8));

    cryptoContext->LevelReduceInPlace(ctxt_cos2, nullptr, 1);
    result2 = cryptoContext->EvalMult(ctxt_cos2, result2);
    cryptoContext->ModReduceInPlace(result2);
    cryptoContext->EvalAddInPlace(result2, result);
    cryptoContext->LevelReduceInPlace(ctxt_cos2, nullptr, 1);
    result2 = cryptoContext->EvalMult(ctxt_cos2, result2);
    cryptoContext->ModReduceInPlace(result2);

    result = cryptoContext->EvalSub(4, ctxt_sin);
    cryptoContext->LevelReduceInPlace(result, nullptr, 4);
    cryptoContext->EvalAddInPlace(result, result2);

    std::cerr << "result levels: " << result->GetLevel() << " and depth: " << result->GetNoiseScaleDeg() << std::endl;

    cryptoContext->Decrypt(keyPair.secretKey, result, &ptxt_result);
    ptxt_result->SetLength(numSlots);
    std::cout << "Output after polynomial approximation for mod 8:" << ptxt_result << std::endl;
    std::vector<double> exact8 = {0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7};

    std::vector<double> error8(ptxt_result->GetRealPackedValue());
    std::transform(exact8.begin(), exact8.end(), error8.begin(), error8.begin(), std::minus<double>());
    std::transform(error8.begin(), error8.end(), error8.begin(), [&](const double& elem) { return std::abs(elem); });
    std::cout << "Absolute error mod 8: " << error8 << std::endl;
    mean_error = std::accumulate(error8.begin(), error8.end(), 0.0) / numSlots;
    std::cout << "Average absolute error mod 8: " << mean_error << std::endl;
    max_error_it = std::max_element(error8.begin(), error8.end());
    std::cout << "Max absolute error mod 8: " << *max_error_it << std::endl << std::endl;

    // double a = -12;
    // double b = 12;
    // int degree = 118;
    // auto coefficients = EvalChebyshevCoefficients([](double x) -> double { return std::cos(2*Pi* x); }, a, b, degree);
    // std::cout.precision(16);
    // std::cout << "\n";
    // std::cout << "coefficients for cos approx of size " << coefficients.size() << ": " << std::endl;
    // for (uint32_t i = 0; i < coefficients.size(); i++) {
    //     std::cout << coefficients[i] << ", ";
    //     if ((i+1) % 4 == 0){
    //         std::cout << "\n";
    //     }
    // }
    // std::cout << std::endl << std::endl;

    // coefficients = EvalChebyshevCoefficients([](double x) -> double { return std::sin(2*Pi*x); }, a, b, degree);
    // std::cout.precision(16);
    // std::cout << "\n";
    // std::cout << "coefficients for sine approx of size " << coefficients.size() << ": " << std::endl;
    // for (uint32_t i = 0; i < coefficients.size(); i++) {
    //     std::cout << coefficients[i] << ", ";
    //     if ((i+1) % 4 == 0){
    //         std::cout << "\n";
    //     }
    // }
    // std::cout << std::endl << std::endl;
}

//------------------------------------------------------------------------------
// Precomputations for CoeffsToSlots and SlotsToCoeffs
//------------------------------------------------------------------------------
void FitToNativeVector(uint32_t ringDim, const std::vector<int64_t>& vec, int64_t bigBound, NativeVector* nativeVec) {
    if (nativeVec == nullptr)
        OPENFHE_THROW("The passed native vector is empty.");
    NativeInteger bigValueHf(bigBound >> 1);
    NativeInteger modulus(nativeVec->GetModulus());
    NativeInteger diff = bigBound - modulus;
    uint32_t dslots    = vec.size();
    uint32_t gap       = ringDim / dslots;
    for (usint i = 0; i < vec.size(); i++) {
        NativeInteger n(vec[i]);
        if (n > bigValueHf) {
            (*nativeVec)[gap * i] = n.ModSub(diff, modulus);
        }
        else {
            (*nativeVec)[gap * i] = n.Mod(modulus);
        }
    }
}

#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
void FitToNativeVector(uint32_t ringDim, const std::vector<int128_t>& vec, int128_t bigBound, NativeVector* nativeVec) {
    if (nativeVec == nullptr)
        OPENFHE_THROW("The passed native vector is empty.");
    NativeInteger bigValueHf((uint128_t)bigBound >> 1);
    NativeInteger modulus(nativeVec->GetModulus());
    NativeInteger diff = NativeInteger((uint128_t)bigBound) - modulus;
    uint32_t dslots    = vec.size();
    uint32_t gap       = ringDim / dslots;
    for (usint i = 0; i < vec.size(); i++) {
        NativeInteger n((uint128_t)vec[i]);
        if (n > bigValueHf) {
            (*nativeVec)[gap * i] = n.ModSub(diff, modulus);
        }
        else {
            (*nativeVec)[gap * i] = n.Mod(modulus);
        }
    }
}
#endif

#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
Plaintext MakeAuxPlaintext(const CryptoContextImpl<DCRTPoly>& cc, const std::shared_ptr<DCRTPoly::Params> params,
                           const std::vector<std::complex<double>>& value, size_t noiseScaleDeg, uint32_t level,
                           usint slots) const {
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc.GetCryptoParameters());

    double scFact = cryptoParams->GetScalingFactorReal(level);

    Plaintext p = Plaintext(std::make_shared<CKKSPackedEncoding>(params, cc.GetEncodingParams(), value, noiseScaleDeg,
                                                                 level, scFact, slots));

    DCRTPoly& plainElement = p->GetElement<DCRTPoly>();

    usint N = cc.GetRingDimension();

    std::vector<std::complex<double>> inverse = value;

    inverse.resize(slots);

    DiscreteFourierTransform::FFTSpecialInv(inverse, N * 2);
    uint64_t pBits = cc.GetEncodingParams()->GetPlaintextModulus();

    double powP      = std::pow(2.0, MAX_DOUBLE_PRECISION);
    int32_t pCurrent = pBits - MAX_DOUBLE_PRECISION;

    std::vector<int128_t> temp(2 * slots);
    for (size_t i = 0; i < slots; ++i) {
        // extract the mantissa of real part and multiply it by 2^52
        int32_t n1 = 0;
        double dre = std::frexp(inverse[i].real(), &n1) * powP;
        // extract the mantissa of imaginary part and multiply it by 2^52
        int32_t n2 = 0;
        double dim = std::frexp(inverse[i].imag(), &n2) * powP;

        // Check for possible overflow
        if (is128BitOverflow(dre) || is128BitOverflow(dim)) {
            DiscreteFourierTransform::FFTSpecial(inverse, N * 2);

            double invLen = static_cast<double>(inverse.size());
            double factor = 2 * M_PI * i;

            double realMax = -1, imagMax = -1;
            uint32_t realMaxIdx = -1, imagMaxIdx = -1;

            for (uint32_t idx = 0; idx < inverse.size(); idx++) {
                // exp( j*2*pi*n*k/N )
                std::complex<double> expFactor = {cos((factor * idx) / invLen), sin((factor * idx) / invLen)};

                // X[k] * exp( j*2*pi*n*k/N )
                std::complex<double> prodFactor = inverse[idx] * expFactor;

                double realVal = prodFactor.real();
                double imagVal = prodFactor.imag();

                if (realVal > realMax) {
                    realMax    = realVal;
                    realMaxIdx = idx;
                }
                if (imagVal > imagMax) {
                    imagMax    = imagVal;
                    imagMaxIdx = idx;
                }
            }

            auto scaledInputSize = ceil(log2(dre));

            std::stringstream buffer;
            buffer << std::endl
                   << "Overflow in data encoding - scaled input is too large to fit "
                      "into a NativeInteger (60 bits). Try decreasing scaling factor."
                   << std::endl;
            buffer << "Overflow at slot number " << i << std::endl;
            buffer << "- Max real part contribution from input[" << realMaxIdx << "]: " << realMax << std::endl;
            buffer << "- Max imaginary part contribution from input[" << imagMaxIdx << "]: " << imagMax << std::endl;
            buffer << "Scaling factor is " << ceil(log2(powP)) << " bits " << std::endl;
            buffer << "Scaled input is " << scaledInputSize << " bits " << std::endl;
            OPENFHE_THROW(buffer.str());
        }

        int64_t re64       = std::llround(dre);
        int32_t pRemaining = pCurrent + n1;
        int128_t re        = 0;
        if (pRemaining < 0) {
            re = re64 >> (-pRemaining);
        }
        else {
            int128_t pPowRemaining = ((int128_t)1) << pRemaining;
            re                     = pPowRemaining * re64;
        }

        int64_t im64 = std::llround(dim);
        pRemaining   = pCurrent + n2;
        int128_t im  = 0;
        if (pRemaining < 0) {
            im = im64 >> (-pRemaining);
        }
        else {
            int128_t pPowRemaining = ((int64_t)1) << pRemaining;
            im                     = pPowRemaining * im64;
        }

        temp[i]         = (re < 0) ? Max128BitValue() + re : re;
        temp[i + slots] = (im < 0) ? Max128BitValue() + im : im;

        if (is128BitOverflow(temp[i]) || is128BitOverflow(temp[i + slots])) {
            OPENFHE_THROW("Overflow, try to decrease scaling factor");
        }
    }

    const std::shared_ptr<ILDCRTParams<BigInteger>> bigParams        = plainElement.GetParams();
    const std::vector<std::shared_ptr<ILNativeParams>>& nativeParams = bigParams->GetParams();

    for (size_t i = 0; i < nativeParams.size(); i++) {
        NativeVector nativeVec(N, nativeParams[i]->GetModulus());
        FitToNativeVector(N, temp, Max128BitValue(), &nativeVec);
        NativePoly element = plainElement.GetElementAtIndex(i);
        element.SetValues(nativeVec, Format::COEFFICIENT);
        plainElement.SetElementAtIndex(i, element);
    }

    usint numTowers = nativeParams.size();
    std::vector<DCRTPoly::Integer> moduli(numTowers);
    for (usint i = 0; i < numTowers; i++) {
        moduli[i] = nativeParams[i]->GetModulus();
    }

    DCRTPoly::Integer intPowP = NativeInteger(1) << pBits;
    std::vector<DCRTPoly::Integer> crtPowP(numTowers, intPowP);

    auto currPowP = crtPowP;

    // We want to scale temp by 2^(pd), and the loop starts from j=2
    // because temp is already scaled by 2^p in the re/im loop above,
    // and currPowP already is 2^p.
    for (size_t i = 2; i < noiseScaleDeg; i++) {
        currPowP = CKKSPackedEncoding::CRTMult(currPowP, crtPowP, moduli);
    }

    if (noiseScaleDeg > 1) {
        plainElement = plainElement.Times(currPowP);
    }

    p->SetFormat(Format::EVALUATION);
    p->SetScalingFactor(pow(p->GetScalingFactor(), noiseScaleDeg));

    return p;
}
#else
Plaintext MakeAuxPlaintext(const CryptoContextImpl<DCRTPoly>& cc, const std::shared_ptr<DCRTPoly::Params> params,
                           const std::vector<std::complex<double>>& value, size_t noiseScaleDeg, uint32_t level,
                           usint slots) {
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc.GetCryptoParameters());

    double scFact = cryptoParams->GetScalingFactorReal(level);

    Plaintext p = Plaintext(std::make_shared<CKKSPackedEncoding>(params, cc.GetEncodingParams(), value, noiseScaleDeg,
                                                                 level, scFact, slots));

    DCRTPoly& plainElement = p->GetElement<DCRTPoly>();

    usint N = cc.GetRingDimension();

    std::vector<std::complex<double>> inverse = value;

    inverse.resize(slots);

    DiscreteFourierTransform::FFTSpecialInv(inverse, N * 2);
    double powP = scFact;

    // Compute approxFactor, a value to scale down by, in case the value exceeds a 64-bit integer.
    constexpr int32_t MAX_BITS_IN_WORD = 61;

    int32_t logc = 0;
    for (size_t i = 0; i < slots; ++i) {
        inverse[i] *= powP;
        if (inverse[i].real() != 0) {
            int32_t logci = static_cast<int32_t>(ceil(log2(std::abs(inverse[i].real()))));
            if (logc < logci)
                logc = logci;
        }
        if (inverse[i].imag() != 0) {
            int32_t logci = static_cast<int32_t>(ceil(log2(std::abs(inverse[i].imag()))));
            if (logc < logci)
                logc = logci;
        }
    }
    if (logc < 0) {
        OPENFHE_THROW("Too small scaling factor");
    }
    int32_t logValid    = (logc <= MAX_BITS_IN_WORD) ? logc : MAX_BITS_IN_WORD;
    int32_t logApprox   = logc - logValid;
    double approxFactor = pow(2, logApprox);

    std::vector<int64_t> temp(2 * slots);

    for (size_t i = 0; i < slots; ++i) {
        // Scale down by approxFactor in case the value exceeds a 64-bit integer.
        double dre = inverse[i].real() / approxFactor;
        double dim = inverse[i].imag() / approxFactor;

        // Check for possible overflow
        if (is64BitOverflow(dre) || is64BitOverflow(dim)) {
            DiscreteFourierTransform::FFTSpecial(inverse, N * 2);

            double invLen = static_cast<double>(inverse.size());
            double factor = 2 * M_PI * i;

            double realMax = -1, imagMax = -1;
            uint32_t realMaxIdx = -1, imagMaxIdx = -1;

            for (uint32_t idx = 0; idx < inverse.size(); idx++) {
                // exp( j*2*pi*n*k/N )
                std::complex<double> expFactor = {cos((factor * idx) / invLen), sin((factor * idx) / invLen)};

                // X[k] * exp( j*2*pi*n*k/N )
                std::complex<double> prodFactor = inverse[idx] * expFactor;

                double realVal = prodFactor.real();
                double imagVal = prodFactor.imag();

                if (realVal > realMax) {
                    realMax    = realVal;
                    realMaxIdx = idx;
                }
                if (imagVal > imagMax) {
                    imagMax    = imagVal;
                    imagMaxIdx = idx;
                }
            }

            auto scaledInputSize = ceil(log2(dre));

            std::stringstream buffer;
            buffer << std::endl
                   << "Overflow in data encoding - scaled input is too large to fit "
                      "into a NativeInteger (60 bits). Try decreasing scaling factor."
                   << std::endl;
            buffer << "Overflow at slot number " << i << std::endl;
            buffer << "- Max real part contribution from input[" << realMaxIdx << "]: " << realMax << std::endl;
            buffer << "- Max imaginary part contribution from input[" << imagMaxIdx << "]: " << imagMax << std::endl;
            buffer << "Scaling factor is " << ceil(log2(powP)) << " bits " << std::endl;
            buffer << "Scaled input is " << scaledInputSize << " bits " << std::endl;
            OPENFHE_THROW(buffer.str());
        }

        int64_t re = std::llround(dre);
        int64_t im = std::llround(dim);

        temp[i]         = (re < 0) ? Max64BitValue() + re : re;
        temp[i + slots] = (im < 0) ? Max64BitValue() + im : im;
    }

    const std::shared_ptr<ILDCRTParams<BigInteger>> bigParams        = plainElement.GetParams();
    const std::vector<std::shared_ptr<ILNativeParams>>& nativeParams = bigParams->GetParams();

    for (size_t i = 0; i < nativeParams.size(); i++) {
        NativeVector nativeVec(N, nativeParams[i]->GetModulus());
        FitToNativeVector(N, temp, Max64BitValue(), &nativeVec);
        NativePoly element = plainElement.GetElementAtIndex(i);
        element.SetValues(nativeVec, Format::COEFFICIENT);
        plainElement.SetElementAtIndex(i, element);
    }

    usint numTowers = nativeParams.size();
    std::vector<DCRTPoly::Integer> moduli(numTowers);
    for (usint i = 0; i < numTowers; i++) {
        moduli[i] = nativeParams[i]->GetModulus();
    }

    DCRTPoly::Integer intPowP{static_cast<uint64_t>(std::llround(powP))};
    std::vector<DCRTPoly::Integer> crtPowP(numTowers, intPowP);

    auto currPowP = crtPowP;

    // We want to scale temp by 2^(pd), and the loop starts from j=2
    // because temp is already scaled by 2^p in the re/im loop above,
    // and currPowP already is 2^p.
    for (size_t i = 2; i < noiseScaleDeg; i++) {
        currPowP = CKKSPackedEncoding::CRTMult(currPowP, crtPowP, moduli);
    }

    if (noiseScaleDeg > 1) {
        plainElement = plainElement.Times(currPowP);
    }

    // Scale back up by the approxFactor to get the correct encoding.
    if (logApprox > 0) {
        int32_t logStep = (logApprox <= MAX_LOG_STEP) ? logApprox : MAX_LOG_STEP;
        auto intStep    = DCRTPoly::Integer(uint64_t(1) << logStep);
        std::vector<DCRTPoly::Integer> crtApprox(numTowers, intStep);
        logApprox -= logStep;

        while (logApprox > 0) {
            logStep = (logApprox <= MAX_LOG_STEP) ? logApprox : MAX_LOG_STEP;
            intStep = DCRTPoly::Integer(uint64_t(1) << logStep);
            std::vector<DCRTPoly::Integer> crtSF(numTowers, intStep);
            crtApprox = CKKSPackedEncoding::CRTMult(crtApprox, crtSF, moduli);
            logApprox -= logStep;
        }
        plainElement = plainElement.Times(crtApprox);
    }

    p->SetFormat(Format::EVALUATION);
    p->SetScalingFactor(pow(p->GetScalingFactor(), noiseScaleDeg));

    return p;
}
#endif

std::vector<ConstPlaintext> EvalLinearTransformPrecompute(const CryptoContextImpl<DCRTPoly>& cc,
                                                          const std::vector<std::vector<std::complex<double>>>& A,
                                                          double scale, uint32_t L) {
    if (A[0].size() != A.size()) {
        OPENFHE_THROW("The matrix passed to EvalLTPrecompute is not square");
    }

    uint32_t slots = A.size();

    uint32_t M = cc.GetCyclotomicOrder();

    // Computing the baby-step bStep and the giant-step gStep.
    int bStep = (m_dim1 == 0) ? ceil(sqrt(slots)) : m_dim1;
    int gStep = ceil(static_cast<double>(slots) / bStep);

    // make sure the plaintext is created only with the necessary amount of moduli

    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc.GetCryptoParameters());

    ILDCRTParams<DCRTPoly::Integer> elementParams = *(cryptoParams->GetElementParams());

    uint32_t towersToDrop = 0;
    if (L != 0) {
        towersToDrop = elementParams.GetParams().size() - L - 1;
    }

    for (uint32_t i = 0; i < towersToDrop; i++) {
        elementParams.PopLastParam();
    }

    auto paramsQ = elementParams.GetParams();
    usint sizeQ  = paramsQ.size();
    auto paramsP = cryptoParams->GetParamsP()->GetParams();
    usint sizeP  = paramsP.size();

    std::vector<NativeInteger> moduli(sizeQ + sizeP);
    std::vector<NativeInteger> roots(sizeQ + sizeP);

    for (size_t i = 0; i < sizeQ; i++) {
        moduli[i] = paramsQ[i]->GetModulus();
        roots[i]  = paramsQ[i]->GetRootOfUnity();
    }

    for (size_t i = 0; i < sizeP; i++) {
        moduli[sizeQ + i] = paramsP[i]->GetModulus();
        roots[sizeQ + i]  = paramsP[i]->GetRootOfUnity();
    }

    auto elementParamsPtr = std::make_shared<ILDCRTParams<DCRTPoly::Integer>>(M, moduli, roots);
    //  auto elementParamsPtr2 = std::dynamic_pointer_cast<typename DCRTPoly::Params>(elementParamsPtr);

    std::vector<ConstPlaintext> result(slots);
// parallelizing the loop (below) with OMP causes a segfault on MinGW
// see https://github.com/openfheorg/openfhe-development/issues/176
#if !defined(__MINGW32__) && !defined(__MINGW64__)
    #pragma omp parallel for
#endif
    for (int j = 0; j < gStep; j++) {
        int offset = -bStep * j;
        for (int i = 0; i < bStep; i++) {
            if (bStep * j + i < static_cast<int>(slots)) {
                auto diag = ExtractShiftedDiagonal(A, bStep * j + i);
                for (uint32_t k = 0; k < diag.size(); k++)
                    diag[k] *= scale;

                result[bStep * j + i] =
                    MakeAuxPlaintext(cc, elementParamsPtr, Rotate(diag, offset), 1, towersToDrop, diag.size());
            }
        }
    }
    return result;
}

std::vector<ConstPlaintext> EvalLinearTransformPrecompute(const CryptoContextImpl<DCRTPoly>& cc,
                                                          const std::vector<std::vector<std::complex<double>>>& A,
                                                          const std::vector<std::vector<std::complex<double>>>& B,
                                                          uint32_t orientation, double scale, uint32_t L) {
    uint32_t slots = A.size();

    uint32_t M = cc.GetCyclotomicOrder();

    // Computing the baby-step bStep and the giant-step gStep.
    int bStep = (m_dim1 == 0) ? ceil(sqrt(slots)) : m_dim1;
    int gStep = ceil(static_cast<double>(slots) / bStep);

    // make sure the plaintext is created only with the necessary amount of moduli

    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc.GetCryptoParameters());

    auto elementParams = *(cryptoParams->GetElementParams());

    uint32_t towersToDrop = 0;
    if (L != 0) {
        towersToDrop = elementParams.GetParams().size() - L - 1;
    }

    for (uint32_t i = 0; i < towersToDrop; i++) {
        elementParams.PopLastParam();
    }

    auto paramsQ = elementParams.GetParams();
    usint sizeQ  = paramsQ.size();
    auto paramsP = cryptoParams->GetParamsP()->GetParams();
    usint sizeP  = paramsP.size();

    std::vector<NativeInteger> moduli(sizeQ + sizeP);
    std::vector<NativeInteger> roots(sizeQ + sizeP);
    for (size_t i = 0; i < sizeQ; i++) {
        moduli[i] = paramsQ[i]->GetModulus();
        roots[i]  = paramsQ[i]->GetRootOfUnity();
    }

    for (size_t i = 0; i < sizeP; i++) {
        moduli[sizeQ + i] = paramsP[i]->GetModulus();
        roots[sizeQ + i]  = paramsP[i]->GetRootOfUnity();
    }

    auto elementParamsPtr = std::make_shared<ILDCRTParams<DCRTPoly::Integer>>(M, moduli, roots);
    //  auto elementParamsPtr2 = std::dynamic_pointer_cast<typename DCRTPoly::Params>(elementParamsPtr);

    std::vector<ConstPlaintext> result(slots);

    if (orientation == 0) {
        // vertical concatenation - used during homomorphic encoding
        // #pragma omp parallel for
        for (int j = 0; j < gStep; j++) {
            int offset = -bStep * j;
            for (int i = 0; i < bStep; i++) {
                if (bStep * j + i < static_cast<int>(slots)) {
                    auto vecA = ExtractShiftedDiagonal(A, bStep * j + i);
                    auto vecB = ExtractShiftedDiagonal(B, bStep * j + i);

                    vecA.insert(vecA.end(), vecB.begin(), vecB.end());
                    for (uint32_t k = 0; k < vecA.size(); k++)
                        vecA[k] *= scale;

                    result[bStep * j + i] =
                        MakeAuxPlaintext(cc, elementParamsPtr, Rotate(vecA, offset), 1, towersToDrop, vecA.size());
                }
            }
        }
    }
    else {
        // horizontal concatenation - used during homomorphic decoding
        std::vector<std::vector<std::complex<double>>> newA(slots);

        //  A and B are concatenated horizontally
        for (uint32_t i = 0; i < slots; ++i) {
            newA[i].reserve(A[i].size() + B[i].size());
            newA[i].insert(newA[i].end(), A[i].begin(), A[i].end());
            newA[i].insert(newA[i].end(), B[i].begin(), B[i].end());
        }

#pragma omp parallel for
        for (int j = 0; j < gStep; j++) {
            int offset = -bStep * j;
            for (int i = 0; i < bStep; i++) {
                if (bStep * j + i < static_cast<int>(slots)) {
                    // shifted diagonal is computed for rectangular map newA of dimension
                    // slots x 2*slots
                    auto vec = ExtractShiftedDiagonal(newA, bStep * j + i);
                    for (uint32_t k = 0; k < vec.size(); k++)
                        vec[k] *= scale;

                    result[bStep * j + i] =
                        MakeAuxPlaintext(cc, elementParamsPtr, Rotate(vec, offset), 1, towersToDrop, vec.size());
                }
            }
        }
    }

    return result;
}

std::vector<std::vector<ConstPlaintext>> EvalCoeffsToSlotsPrecompute(const CryptoContextImpl<DCRTPoly>& cc,
                                                                     const std::vector<std::complex<double>>& A,
                                                                     const std::vector<uint32_t>& rotGroup, bool flag_i,
                                                                     double scale, uint32_t L) {
    uint32_t slots = rotGroup.size();

    uint32_t M = cc.GetCyclotomicOrder();

    int32_t levelBudget     = m_paramsEnc[CKKS_BOOT_PARAMS::LEVEL_BUDGET];  // Andreea: this is not currently filled
    int32_t layersCollapse  = m_paramsEnc[CKKS_BOOT_PARAMS::LAYERS_COLL];
    int32_t remCollapse     = m_paramsEnc[CKKS_BOOT_PARAMS::LAYERS_REM];
    int32_t numRotations    = m_paramsEnc[CKKS_BOOT_PARAMS::NUM_ROTATIONS];
    int32_t b               = m_paramsEnc[CKKS_BOOT_PARAMS::BABY_STEP];
    int32_t g               = m_paramsEnc[CKKS_BOOT_PARAMS::GIANT_STEP];
    int32_t numRotationsRem = m_paramsEnc[CKKS_BOOT_PARAMS::NUM_ROTATIONS_REM];
    int32_t bRem            = m_paramsEnc[CKKS_BOOT_PARAMS::BABY_STEP_REM];
    int32_t gRem            = m_paramsEnc[CKKS_BOOT_PARAMS::GIANT_STEP_REM];

    int32_t stop    = -1;
    int32_t flagRem = 0;

    if (remCollapse != 0) {
        stop    = 0;
        flagRem = 1;
    }

    // result is the rotated plaintext version of the coefficients
    std::vector<std::vector<ConstPlaintext>> result(levelBudget);
    for (uint32_t i = 0; i < uint32_t(levelBudget); i++) {
        if (flagRem == 1 && i == 0) {
            // remainder corresponds to index 0 in encoding and to last index in decoding
            result[i] = std::vector<ConstPlaintext>(numRotationsRem);
        }
        else {
            result[i] = std::vector<ConstPlaintext>(numRotations);
        }
    }

    // make sure the plaintext is created only with the necessary amount of moduli

    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc.GetCryptoParameters());

    auto elementParams = *(cryptoParams->GetElementParams());

    uint32_t towersToDrop = 0;

    if (L != 0) {
        towersToDrop = elementParams.GetParams().size() - L - levelBudget;
    }

    for (uint32_t i = 0; i < towersToDrop; i++) {
        elementParams.PopLastParam();
    }

    uint32_t level0 = towersToDrop + levelBudget - 1;

    auto paramsQ = elementParams.GetParams();
    usint sizeQ  = paramsQ.size();
    auto paramsP = cryptoParams->GetParamsP()->GetParams();
    usint sizeP  = paramsP.size();

    std::vector<NativeInteger> moduli(sizeQ + sizeP);
    std::vector<NativeInteger> roots(sizeQ + sizeP);
    for (size_t i = 0; i < sizeQ; i++) {
        moduli[i] = paramsQ[i]->GetModulus();
        roots[i]  = paramsQ[i]->GetRootOfUnity();
    }

    for (size_t i = 0; i < sizeP; i++) {
        moduli[sizeQ + i] = paramsP[i]->GetModulus();
        roots[sizeQ + i]  = paramsP[i]->GetRootOfUnity();
    }

    // we need to pre-compute the plaintexts in the extended basis P*Q
    std::vector<std::shared_ptr<ILDCRTParams<BigInteger>>> paramsVector(levelBudget - stop);
    for (int32_t s = levelBudget - 1; s >= stop; s--) {
        paramsVector[s - stop] = std::make_shared<ILDCRTParams<BigInteger>>(M, moduli, roots);
        moduli.erase(moduli.begin() + sizeQ - 1);
        roots.erase(roots.begin() + sizeQ - 1);
        sizeQ--;
    }

    if (slots == M / 4) {
        //------------------------------------------------------------------------------
        // fully-packed mode
        //------------------------------------------------------------------------------

        auto coeff = CoeffEncodingCollapse(A, rotGroup, levelBudget, flag_i);

        for (int32_t s = levelBudget - 1; s > stop; s--) {
            for (int32_t i = 0; i < b; i++) {
#if !defined(__MINGW32__) && !defined(__MINGW64__)
    #pragma omp parallel for
#endif
                for (int32_t j = 0; j < g; j++) {
                    if (g * i + j != int32_t(numRotations)) {
                        uint32_t rot =
                            ReduceRotation(-g * i * (1 << ((s - flagRem) * layersCollapse + remCollapse)), slots);
                        if ((flagRem == 0) && (s == stop + 1)) {
                            // do the scaling only at the last set of coefficients
                            for (uint32_t k = 0; k < slots; k++) {
                                coeff[s][g * i + j][k] *= scale;
                            }
                        }

                        auto rotateTemp = Rotate(coeff[s][g * i + j], rot);

                        result[s][g * i + j] =
                            MakeAuxPlaintext(cc, paramsVector[s - stop], rotateTemp, 1, level0 - s, rotateTemp.size());
                    }
                }
            }
        }

        if (flagRem) {
            for (int32_t i = 0; i < bRem; i++) {
#pragma omp parallel for
                for (int32_t j = 0; j < gRem; j++) {
                    if (gRem * i + j != int32_t(numRotationsRem)) {
                        uint32_t rot = ReduceRotation(-gRem * i, slots);
                        for (uint32_t k = 0; k < slots; k++) {
                            coeff[stop][gRem * i + j][k] *= scale;
                        }

                        auto rotateTemp = Rotate(coeff[stop][gRem * i + j], rot);
                        result[stop][gRem * i + j] =
                            MakeAuxPlaintext(cc, paramsVector[0], rotateTemp, 1, level0, rotateTemp.size());
                    }
                }
            }
        }
    }
    else {
        //------------------------------------------------------------------------------
        // sparsely-packed mode
        //------------------------------------------------------------------------------

        auto coeff  = CoeffEncodingCollapse(A, rotGroup, levelBudget, false);
        auto coeffi = CoeffEncodingCollapse(A, rotGroup, levelBudget, true);

        for (int32_t s = levelBudget - 1; s > stop; s--) {
            for (int32_t i = 0; i < b; i++) {
#if !defined(__MINGW32__) && !defined(__MINGW64__)
    #pragma omp parallel for
#endif
                for (int32_t j = 0; j < g; j++) {
                    if (g * i + j != int32_t(numRotations)) {
                        uint32_t rot =
                            ReduceRotation(-g * i * (1 << ((s - flagRem) * layersCollapse + remCollapse)), M / 4);
                        // concatenate the coefficients horizontally on their third dimension, which corresponds to the # of slots
                        auto clearTemp  = coeff[s][g * i + j];
                        auto clearTempi = coeffi[s][g * i + j];
                        clearTemp.insert(clearTemp.end(), clearTempi.begin(), clearTempi.end());
                        if ((flagRem == 0) && (s == stop + 1)) {
                            // do the scaling only at the last set of coefficients
                            for (uint32_t k = 0; k < clearTemp.size(); k++) {
                                clearTemp[k] *= scale;
                            }
                        }

                        auto rotateTemp = Rotate(clearTemp, rot);
                        result[s][g * i + j] =
                            MakeAuxPlaintext(cc, paramsVector[s - stop], rotateTemp, 1, level0 - s, rotateTemp.size());
                    }
                }
            }
        }

        if (flagRem) {
            for (int32_t i = 0; i < bRem; i++) {
#pragma omp parallel for
                for (int32_t j = 0; j < gRem; j++) {
                    if (gRem * i + j != int32_t(numRotationsRem)) {
                        uint32_t rot = ReduceRotation(-gRem * i, M / 4);
                        // concatenate the coefficients on their third dimension, which corresponds to the # of slots
                        auto clearTemp  = coeff[stop][gRem * i + j];
                        auto clearTempi = coeffi[stop][gRem * i + j];
                        clearTemp.insert(clearTemp.end(), clearTempi.begin(), clearTempi.end());
                        for (uint32_t k = 0; k < clearTemp.size(); k++) {
                            clearTemp[k] *= scale;
                        }

                        auto rotateTemp = Rotate(clearTemp, rot);
                        result[stop][gRem * i + j] =
                            MakeAuxPlaintext(cc, paramsVector[0], rotateTemp, 1, level0, rotateTemp.size());
                    }
                }
            }
        }
    }
    return result;
}

std::vector<std::vector<ConstPlaintext>> EvalSlotsToCoeffsPrecompute(const CryptoContextImpl<DCRTPoly>& cc,
                                                                     const std::vector<std::complex<double>>& A,
                                                                     const std::vector<uint32_t>& rotGroup, bool flag_i,
                                                                     double scale, uint32_t L) {
    uint32_t slots = rotGroup.size();

    uint32_t M = cc.GetCyclotomicOrder();

    int32_t levelBudget     = m_paramsDec[CKKS_BOOT_PARAMS::LEVEL_BUDGET];  // Andreea: this is not currently filled
    int32_t layersCollapse  = m_paramsDec[CKKS_BOOT_PARAMS::LAYERS_COLL];
    int32_t remCollapse     = m_paramsDec[CKKS_BOOT_PARAMS::LAYERS_REM];
    int32_t numRotations    = m_paramsDec[CKKS_BOOT_PARAMS::NUM_ROTATIONS];
    int32_t b               = m_paramsDec[CKKS_BOOT_PARAMS::BABY_STEP];
    int32_t g               = m_paramsDec[CKKS_BOOT_PARAMS::GIANT_STEP];
    int32_t numRotationsRem = m_paramsDec[CKKS_BOOT_PARAMS::NUM_ROTATIONS_REM];
    int32_t bRem            = m_paramsDec[CKKS_BOOT_PARAMS::BABY_STEP_REM];
    int32_t gRem            = m_paramsDec[CKKS_BOOT_PARAMS::GIANT_STEP_REM];

    int32_t flagRem = 0;

    if (remCollapse != 0) {
        flagRem = 1;
    }

    // result is the rotated plaintext version of coeff
    std::vector<std::vector<ConstPlaintext>> result(levelBudget);
    for (uint32_t i = 0; i < uint32_t(levelBudget); i++) {
        if (flagRem == 1 && i == uint32_t(levelBudget - 1)) {
            // remainder corresponds to index 0 in encoding and to last index in decoding
            result[i] = std::vector<ConstPlaintext>(numRotationsRem);
        }
        else {
            result[i] = std::vector<ConstPlaintext>(numRotations);
        }
    }

    // make sure the plaintext is created only with the necessary amount of moduli

    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc.GetCryptoParameters());

    auto elementParams = *(cryptoParams->GetElementParams());

    uint32_t towersToDrop = 0;

    if (L != 0) {
        towersToDrop = elementParams.GetParams().size() - L - levelBudget;
    }

    for (uint32_t i = 0; i < towersToDrop; i++) {
        elementParams.PopLastParam();
    }

    uint32_t level0 = towersToDrop;

    auto paramsQ = elementParams.GetParams();
    usint sizeQ  = paramsQ.size();
    auto paramsP = cryptoParams->GetParamsP()->GetParams();
    usint sizeP  = paramsP.size();

    std::vector<NativeInteger> moduli(sizeQ + sizeP);
    std::vector<NativeInteger> roots(sizeQ + sizeP);
    for (size_t i = 0; i < sizeQ; i++) {
        moduli[i] = paramsQ[i]->GetModulus();
        roots[i]  = paramsQ[i]->GetRootOfUnity();
    }

    for (size_t i = 0; i < sizeP; i++) {
        moduli[sizeQ + i] = paramsP[i]->GetModulus();
        roots[sizeQ + i]  = paramsP[i]->GetRootOfUnity();
    }

    // we need to pre-compute the plaintexts in the extended basis P*Q
    std::vector<std::shared_ptr<ILDCRTParams<BigInteger>>> paramsVector(levelBudget - flagRem + 1);
    for (int32_t s = 0; s < levelBudget - flagRem + 1; s++) {
        paramsVector[s] = std::make_shared<ILDCRTParams<BigInteger>>(M, moduli, roots);
        moduli.erase(moduli.begin() + sizeQ - 1);
        roots.erase(roots.begin() + sizeQ - 1);
        sizeQ--;
    }

    if (slots == M / 4) {
        // fully-packed
        auto coeff = CoeffDecodingCollapse(A, rotGroup, levelBudget, flag_i);

        for (int32_t s = 0; s < levelBudget - flagRem; s++) {
            for (int32_t i = 0; i < b; i++) {
#pragma omp parallel for
                for (int32_t j = 0; j < g; j++) {
                    if (g * i + j != int32_t(numRotations)) {
                        uint32_t rot = ReduceRotation(-g * i * (1 << (s * layersCollapse)), slots);
                        if ((flagRem == 0) && (s == levelBudget - flagRem - 1)) {
                            // do the scaling only at the last set of coefficients
                            for (uint32_t k = 0; k < slots; k++) {
                                coeff[s][g * i + j][k] *= scale;
                            }
                        }

                        auto rotateTemp = Rotate(coeff[s][g * i + j], rot);
                        result[s][g * i + j] =
                            MakeAuxPlaintext(cc, paramsVector[s], rotateTemp, 1, level0 + s, rotateTemp.size());
                    }
                }
            }
        }

        if (flagRem) {
            int32_t s = levelBudget - flagRem;
            for (int32_t i = 0; i < bRem; i++) {
#pragma omp parallel for
                for (int32_t j = 0; j < gRem; j++) {
                    if (gRem * i + j != int32_t(numRotationsRem)) {
                        uint32_t rot = ReduceRotation(-gRem * i * (1 << (s * layersCollapse)), slots);
                        for (uint32_t k = 0; k < slots; k++) {
                            coeff[s][gRem * i + j][k] *= scale;
                        }

                        auto rotateTemp = Rotate(coeff[s][gRem * i + j], rot);
                        result[s][gRem * i + j] =
                            MakeAuxPlaintext(cc, paramsVector[s], rotateTemp, 1, level0 + s, rotateTemp.size());
                    }
                }
            }
        }
    }
    else {
        //------------------------------------------------------------------------------
        // sparsely-packed mode
        //------------------------------------------------------------------------------

        auto coeff  = CoeffDecodingCollapse(A, rotGroup, levelBudget, false);
        auto coeffi = CoeffDecodingCollapse(A, rotGroup, levelBudget, true);

        for (int32_t s = 0; s < levelBudget - flagRem; s++) {
            for (int32_t i = 0; i < b; i++) {
#pragma omp parallel for
                for (int32_t j = 0; j < g; j++) {
                    if (g * i + j != int32_t(numRotations)) {
                        uint32_t rot = ReduceRotation(-g * i * (1 << (s * layersCollapse)), M / 4);
                        // concatenate the coefficients horizontally on their third dimension, which corresponds to the # of slots
                        auto clearTemp  = coeff[s][g * i + j];
                        auto clearTempi = coeffi[s][g * i + j];
                        clearTemp.insert(clearTemp.end(), clearTempi.begin(), clearTempi.end());
                        if ((flagRem == 0) && (s == levelBudget - flagRem - 1)) {
                            // do the scaling only at the last set of coefficients
                            for (uint32_t k = 0; k < clearTemp.size(); k++) {
                                clearTemp[k] *= scale;
                            }
                        }

                        auto rotateTemp = Rotate(clearTemp, rot);
                        result[s][g * i + j] =
                            MakeAuxPlaintext(cc, paramsVector[s], rotateTemp, 1, level0 + s, rotateTemp.size());
                    }
                }
            }
        }

        if (flagRem) {
            int32_t s = levelBudget - flagRem;
            for (int32_t i = 0; i < bRem; i++) {
#pragma omp parallel for
                for (int32_t j = 0; j < gRem; j++) {
                    if (gRem * i + j != int32_t(numRotationsRem)) {
                        uint32_t rot = ReduceRotation(-gRem * i * (1 << (s * layersCollapse)), M / 4);
                        // concatenate the coefficients horizontally on their third dimension, which corresponds to the # of slots
                        auto clearTemp  = coeff[s][gRem * i + j];
                        auto clearTempi = coeffi[s][gRem * i + j];
                        clearTemp.insert(clearTemp.end(), clearTempi.begin(), clearTempi.end());
                        for (uint32_t k = 0; k < clearTemp.size(); k++) {
                            clearTemp[k] *= scale;
                        }

                        auto rotateTemp = Rotate(clearTemp, rot);
                        result[s][gRem * i + j] =
                            MakeAuxPlaintext(cc, paramsVector[s], rotateTemp, 1, level0 + s, rotateTemp.size());
                    }
                }
            }
        }
    }
    return result;
}

void EvalFuncBTSetup(const CryptoContextImpl<DCRTPoly>& cc, uint32_t numSlots, uint32_t digitBitSize,
                     std::vector<uint32_t> dim1, std::vector<uint32_t> levelBudget, double scaleMod) {
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc.GetCryptoParameters());

    if (cryptoParams->GetKeySwitchTechnique() != HYBRID)
        OPENFHE_THROW("CKKS Bootstrapping is only supported for the Hybrid key switching method.");
#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
    if (cryptoParams->GetScalingTechnique() == FLEXIBLEAUTO || cryptoParams->GetScalingTechnique() == FLEXIBLEAUTOEXT)
        OPENFHE_THROW("128-bit CKKS Bootstrapping is supported for FIXEDMANUAL and FIXEDAUTO methods only.");
#endif

    uint32_t M     = cc.GetCyclotomicOrder();
    uint32_t slots = (numSlots == 0) ? M / 4 : numSlots;

    m_dim1     = dim1[0];
    m_gs       = dim1[1];
    m_levelEnc = levelBudget[0];
    m_levelDec = levelBudget[1];

    m_paramsEnc = GetCollapsedFFTParams(slots, levelBudget[0], dim1[0]);
    m_paramsDec = GetCollapsedFFTParams(slots, levelBudget[1], dim1[1]);

    uint32_t m    = 4 * slots;
    bool isSparse = (M != m) ? true : false;

    // computes indices for all primitive roots of unity
    std::vector<uint32_t> rotGroup(slots);
    uint32_t fivePows = 1;
    for (uint32_t i = 0; i < slots; ++i) {
        rotGroup[i] = fivePows;
        fivePows *= 5;
        fivePows %= m;
    }

    // computes all powers of a primitive root of unity exp(2 * M_PI/m)
    std::vector<std::complex<double>> ksiPows(m + 1);
    for (uint32_t j = 0; j < m; ++j) {
        double angle = 2.0 * M_PI * j / m;
        ksiPows[j].real(cos(angle));
        ksiPows[j].imag(sin(angle));
    }
    ksiPows[m] = ksiPows[0];

    // Extract the modulus prior to bootstrapping
    NativeInteger q  = cryptoParams->GetElementParams()->GetParams()[0]->GetModulus().ConvertToInt();
    double qDouble   = q.ConvertToDouble();
    uint128_t factor = ((uint128_t)1 << ((uint32_t)std::round(std::log2(qDouble))));
    double pre       = qDouble / factor;
    double k         = (cryptoParams->GetSecretKeyDist() == SPARSE_TERNARY) ? K_SPARSE : 1.0;
    std::cerr << "pre in setup = q / factor = " << pre << ", k = " << k << std::endl;
    double scaleEnc = pre / k;
    // double scaleDec  = 1.0; // 1.0 / pre;
    double scaleDec = scaleMod / pre;

    uint32_t approxModDepth = 8;  // Andreea: automate?
    uint32_t extraDepth(0);
    switch (digitBitSize) {
        case 1:
            extraDepth = 0;
            break;
        case 2:
            extraDepth = 2;
            break;
        case 3:
            extraDepth = 4;
            break;
        default:
            OPENFHE_THROW("Digit sizes of more than 3 bits are not currently allowed.");
            break;
    }
    uint32_t depthBT = approxModDepth + m_levelEnc + m_levelDec + extraDepth;

    // compute # of levels to remain when encoding the coefficients
    uint32_t L0 = cryptoParams->GetElementParams()->GetParams().size();
    // for FLEXIBLEAUTOEXT we do not need extra modulus in auxiliary plaintexts
    if (cryptoParams->GetScalingTechnique() == FLEXIBLEAUTOEXT)
        L0 -= 1;
    uint32_t lEnc = L0 - m_levelEnc - 1;
    uint32_t lDec = L0 - depthBT;

    bool isLTBootstrap = (m_levelEnc == 1) && (m_levelDec == 1);

    if (isLTBootstrap) {
        // allocate all vectors
        std::vector<std::vector<std::complex<double>>> U0(slots, std::vector<std::complex<double>>(slots));
        std::vector<std::vector<std::complex<double>>> U1(slots, std::vector<std::complex<double>>(slots));
        std::vector<std::vector<std::complex<double>>> U0hatT(slots, std::vector<std::complex<double>>(slots));
        std::vector<std::vector<std::complex<double>>> U1hatT(slots, std::vector<std::complex<double>>(slots));

        for (size_t i = 0; i < slots; i++) {
            for (size_t j = 0; j < slots; j++) {
                U0[i][j]     = ksiPows[(j * rotGroup[i]) % m];
                U0hatT[j][i] = std::conj(U0[i][j]);
                U1[i][j]     = std::complex<double>(0, 1) * U0[i][j];
                U1hatT[j][i] = std::conj(U1[i][j]);
            }
        }

        if (!isSparse) {
            m_U0hatTPre = EvalLinearTransformPrecompute(cc, U0hatT, scaleEnc, lEnc);
            m_U0Pre     = EvalLinearTransformPrecompute(cc, U0, scaleDec, lDec);
        }
        else {
            m_U0hatTPre = EvalLinearTransformPrecompute(cc, U0hatT, U1hatT, 0, scaleEnc, lEnc);
            m_U0Pre     = EvalLinearTransformPrecompute(cc, U0, U1, 1, scaleDec, lDec);
        }
    }
    else {
        m_U0hatTPreFFT = EvalCoeffsToSlotsPrecompute(cc, ksiPows, rotGroup, false, scaleEnc, lEnc);
        m_U0PreFFT     = EvalSlotsToCoeffsPrecompute(cc, ksiPows, rotGroup, false, scaleDec, lDec);
    }
}

void EvalFuncBTKeyGen(const PrivateKey<DCRTPoly> privateKey, uint32_t slots) {
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(privateKey->GetCryptoParameters());

    if (cryptoParams->GetKeySwitchTechnique() != HYBRID)
        OPENFHE_THROW("CKKS Bootstrapping is only supported for the Hybrid key switching method.");
#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
    if (cryptoParams->GetScalingTechnique() == FLEXIBLEAUTO || cryptoParams->GetScalingTechnique() == FLEXIBLEAUTOEXT)
        OPENFHE_THROW("128-bit CKKS Bootstrapping is supported for FIXEDMANUAL and FIXEDAUTO methods only.");
#endif
    auto cc    = privateKey->GetCryptoContext();
    uint32_t M = cc->GetCyclotomicOrder();

    if (slots == 0)
        slots = M / 4;
    // computing all indices for baby-step giant-step procedure
    auto algo = cc->GetScheme();
    // auto evalKeys = algo->EvalAtIndexKeyGen(nullptr, privateKey, FindBootstrapRotationIndices(slots, M));
    cc->EvalAtIndexKeyGen(privateKey, FindBootstrapRotationIndices(slots, M));

    m_conjKey = ConjugateKeyGen(privateKey);
    // (*evalKeys)[M - 1] = conjKey;

    // return evalKeys;
}

EvalKey<DCRTPoly> ConjugateKeyGen(const PrivateKey<DCRTPoly> privateKey) {
    const auto cc = privateKey->GetCryptoContext();
    auto algo     = cc->GetScheme();

    const DCRTPoly& s = privateKey->GetPrivateElement();
    usint N           = s.GetRingDimension();

    PrivateKey<DCRTPoly> privateKeyPermuted = std::make_shared<PrivateKeyImpl<DCRTPoly>>(cc);

    usint index = 2 * N - 1;
    std::vector<usint> vec(N);
    PrecomputeAutoMap(N, index, &vec);

    DCRTPoly sPermuted = s.AutomorphismTransform(index, vec);

    privateKeyPermuted->SetPrivateElement(sPermuted);
    privateKeyPermuted->SetKeyTag(privateKey->GetKeyTag());

    return algo->KeySwitchGen(privateKey, privateKeyPermuted);
}

//------------------------------------------------------------------------------
// Find Rotation Indices
//------------------------------------------------------------------------------

std::vector<int32_t> FindBootstrapRotationIndices(uint32_t slots, uint32_t M) {
    std::vector<int32_t> fullIndexList;

    bool isLTBootstrap = (m_levelEnc == 1) && (m_levelDec == 1);

    if (isLTBootstrap) {
        fullIndexList = FindLinearTransformRotationIndices(slots, M);
    }
    else {
        fullIndexList = FindCoeffsToSlotsRotationIndices(slots, M);

        std::vector<int32_t> indexListStC = FindSlotsToCoeffsRotationIndices(slots, M);
        fullIndexList.insert(fullIndexList.end(), indexListStC.begin(), indexListStC.end());
    }

    // Remove possible duplicates
    sort(fullIndexList.begin(), fullIndexList.end());
    fullIndexList.erase(unique(fullIndexList.begin(), fullIndexList.end()), fullIndexList.end());

    // remove automorphisms corresponding to 0
    fullIndexList.erase(std::remove(fullIndexList.begin(), fullIndexList.end(), 0), fullIndexList.end());
    fullIndexList.erase(std::remove(fullIndexList.begin(), fullIndexList.end(), M / 4), fullIndexList.end());

    return fullIndexList;
}

std::vector<int32_t> FindLinearTransformRotationIndices(uint32_t slots, uint32_t M) {
    std::vector<int32_t> indexList;

    // Computing the baby-step g and the giant-step h.
    int g = (m_dim1 == 0) ? ceil(sqrt(slots)) : m_dim1;
    int h = ceil(static_cast<double>(slots) / g);

    // computing all indices for baby-step giant-step procedure
    // ATTN: resize() is used as indexListEvalLT may be empty here
    indexList.reserve(g + h + M - 2);
    for (int i = 0; i < g; i++) {
        indexList.emplace_back(i + 1);
    }
    for (int i = 2; i < h; i++) {
        indexList.emplace_back(g * i);
    }

    uint32_t m = slots * 4;
    // additional automorphisms are needed for sparse bootstrapping
    if (m != M) {
        for (uint32_t j = 1; j < M / m; j <<= 1) {
            indexList.emplace_back(j * slots);
        }
    }
    // Remove possible duplicates
    sort(indexList.begin(), indexList.end());
    indexList.erase(unique(indexList.begin(), indexList.end()), indexList.end());

    // remove automorphisms corresponding to 0
    indexList.erase(std::remove(indexList.begin(), indexList.end(), 0), indexList.end());
    indexList.erase(std::remove(indexList.begin(), indexList.end(), M / 4), indexList.end());

    return indexList;
}

std::vector<int32_t> FindCoeffsToSlotsRotationIndices(uint32_t slots, uint32_t M) {
    std::vector<int32_t> indexList;

    int32_t levelBudget     = m_paramsEnc[CKKS_BOOT_PARAMS::LEVEL_BUDGET];  // Andreea: these haven't been filled yet
    int32_t layersCollapse  = m_paramsEnc[CKKS_BOOT_PARAMS::LAYERS_COLL];
    int32_t remCollapse     = m_paramsEnc[CKKS_BOOT_PARAMS::LAYERS_REM];
    int32_t numRotations    = m_paramsEnc[CKKS_BOOT_PARAMS::NUM_ROTATIONS];
    int32_t b               = m_paramsEnc[CKKS_BOOT_PARAMS::BABY_STEP];
    int32_t g               = m_paramsEnc[CKKS_BOOT_PARAMS::GIANT_STEP];
    int32_t numRotationsRem = m_paramsEnc[CKKS_BOOT_PARAMS::NUM_ROTATIONS_REM];
    int32_t bRem            = m_paramsEnc[CKKS_BOOT_PARAMS::BABY_STEP_REM];
    int32_t gRem            = m_paramsEnc[CKKS_BOOT_PARAMS::GIANT_STEP_REM];

    int32_t stop;
    int32_t flagRem;
    if (remCollapse == 0) {
        stop    = -1;
        flagRem = 0;
    }
    else {
        stop    = 0;
        flagRem = 1;
    }

    // Computing all indices for baby-step giant-step procedure for encoding and decoding
    indexList.reserve(b + g - 2 + bRem + gRem - 2 + 1 + M);

    for (int32_t s = int32_t(levelBudget) - 1; s > stop; s--) {
        for (int32_t j = 0; j < g; j++) {
            indexList.emplace_back(ReduceRotation(
                (j - int32_t((numRotations + 1) / 2) + 1) * (1 << ((s - flagRem) * layersCollapse + remCollapse)),
                slots));
        }

        for (int32_t i = 0; i < b; i++) {
            indexList.emplace_back(
                ReduceRotation((g * i) * (1 << ((s - flagRem) * layersCollapse + remCollapse)), M / 4));
        }
    }

    if (flagRem) {
        for (int32_t j = 0; j < gRem; j++) {
            indexList.emplace_back(ReduceRotation((j - int32_t((numRotationsRem + 1) / 2) + 1), slots));
        }
        for (int32_t i = 0; i < bRem; i++) {
            indexList.emplace_back(ReduceRotation(gRem * i, M / 4));
        }
    }

    uint32_t m = slots * 4;
    // additional automorphisms are needed for sparse bootstrapping
    if (m != M) {
        for (uint32_t j = 1; j < M / m; j <<= 1) {
            indexList.emplace_back(j * slots);
        }
    }

    // Remove possible duplicates
    sort(indexList.begin(), indexList.end());
    indexList.erase(unique(indexList.begin(), indexList.end()), indexList.end());

    // remove automorphisms corresponding to 0
    indexList.erase(std::remove(indexList.begin(), indexList.end(), 0), indexList.end());
    indexList.erase(std::remove(indexList.begin(), indexList.end(), M / 4), indexList.end());

    return indexList;
}

std::vector<int32_t> FindSlotsToCoeffsRotationIndices(uint32_t slots, uint32_t M) {
    std::vector<int32_t> indexList;

    int32_t levelBudget     = m_paramsDec[CKKS_BOOT_PARAMS::LEVEL_BUDGET];  // Andreea: these haven't been filled yet
    int32_t layersCollapse  = m_paramsDec[CKKS_BOOT_PARAMS::LAYERS_COLL];
    int32_t remCollapse     = m_paramsDec[CKKS_BOOT_PARAMS::LAYERS_REM];
    int32_t numRotations    = m_paramsDec[CKKS_BOOT_PARAMS::NUM_ROTATIONS];
    int32_t b               = m_paramsDec[CKKS_BOOT_PARAMS::BABY_STEP];
    int32_t g               = m_paramsDec[CKKS_BOOT_PARAMS::GIANT_STEP];
    int32_t numRotationsRem = m_paramsDec[CKKS_BOOT_PARAMS::NUM_ROTATIONS_REM];
    int32_t bRem            = m_paramsDec[CKKS_BOOT_PARAMS::BABY_STEP_REM];
    int32_t gRem            = m_paramsDec[CKKS_BOOT_PARAMS::GIANT_STEP_REM];

    int32_t flagRem;
    if (remCollapse == 0) {
        flagRem = 0;
    }
    else {
        flagRem = 1;
    }

    // Computing all indices for baby-step giant-step procedure for encoding and decoding
    indexList.reserve(b + g - 2 + bRem + gRem - 2 + 1 + M);

    for (int32_t s = 0; s < int32_t(levelBudget); s++) {
        for (int32_t j = 0; j < g; j++) {
            indexList.emplace_back(
                ReduceRotation((j - (numRotations + 1) / 2 + 1) * (1 << (s * layersCollapse)), M / 4));
        }
        for (int32_t i = 0; i < b; i++) {
            indexList.emplace_back(ReduceRotation((g * i) * (1 << (s * layersCollapse)), M / 4));
        }
    }

    if (flagRem) {
        int32_t s = int32_t(levelBudget) - flagRem;
        for (int32_t j = 0; j < gRem; j++) {
            indexList.emplace_back(
                ReduceRotation((j - (numRotationsRem + 1) / 2 + 1) * (1 << (s * layersCollapse)), M / 4));
        }
        for (int32_t i = 0; i < bRem; i++) {
            indexList.emplace_back(ReduceRotation((gRem * i) * (1 << (s * layersCollapse)), M / 4));
        }
    }

    uint32_t m = slots * 4;
    // additional automorphisms are needed for sparse bootstrapping
    if (m != M) {
        for (uint32_t j = 1; j < M / m; j <<= 1) {
            indexList.emplace_back(j * slots);
        }
    }

    // Remove possible duplicates
    sort(indexList.begin(), indexList.end());
    indexList.erase(unique(indexList.begin(), indexList.end()), indexList.end());

    // remove automorphisms corresponding to 0
    indexList.erase(std::remove(indexList.begin(), indexList.end(), 0), indexList.end());
    indexList.erase(std::remove(indexList.begin(), indexList.end(), M / 4), indexList.end());

    return indexList;
}

//------------------------------------------------------------------------------
// EVALUATION: CoeffsToSlots and SlotsToCoeffs
//------------------------------------------------------------------------------

Ciphertext<DCRTPoly> EvalLinearTransform(const std::vector<ConstPlaintext>& A, ConstCiphertext<DCRTPoly> ct) {
    uint32_t slots = A.size();

    auto cc = ct->GetCryptoContext();
    // Computing the baby-step bStep and the giant-step gStep.
    uint32_t bStep = (m_dim1 == 0) ? ceil(sqrt(slots)) : m_dim1;
    uint32_t gStep = ceil(static_cast<double>(slots) / bStep);

    uint32_t M = cc->GetCyclotomicOrder();
    uint32_t N = cc->GetRingDimension();

    // computes the NTTs for each CRT limb (for the hoisted automorphisms used
    // later on)
    auto digits = cc->EvalFastRotationPrecompute(ct);

    std::vector<Ciphertext<DCRTPoly>> fastRotation(bStep - 1);

    // hoisted automorphisms
#pragma omp parallel for
    for (uint32_t j = 1; j < bStep; j++) {
        fastRotation[j - 1] = cc->EvalFastRotationExt(ct, j, digits, true);
    }

    Ciphertext<DCRTPoly> result;
    DCRTPoly first;

    for (uint32_t j = 0; j < gStep; j++) {
        Ciphertext<DCRTPoly> inner = EvalMultExt(cc->KeySwitchExt(ct, true), A[bStep * j]);
        for (uint32_t i = 1; i < bStep; i++) {
            if (bStep * j + i < slots) {
                EvalAddExtInPlace(inner, EvalMultExt(fastRotation[i - 1], A[bStep * j + i]));
            }
        }

        if (j == 0) {
            first         = cc->KeySwitchDownFirstElement(inner);
            auto elements = inner->GetElements();
            elements[0].SetValuesToZero();
            inner->SetElements(elements);
            result = inner;
        }
        else {
            inner = cc->KeySwitchDown(inner);
            // Find the automorphism index that corresponds to rotation index index.
            usint autoIndex = FindAutomorphismIndex2nComplex(bStep * j, M);
            std::vector<usint> map(N);
            PrecomputeAutoMap(N, autoIndex, &map);
            DCRTPoly firstCurrent = inner->GetElements()[0].AutomorphismTransform(autoIndex, map);
            first += firstCurrent;

            auto innerDigits = cc->EvalFastRotationPrecompute(inner);
            EvalAddExtInPlace(result, cc->EvalFastRotationExt(inner, bStep * j, innerDigits, false));
        }
    }

    result        = cc->KeySwitchDown(result);
    auto elements = result->GetElements();
    elements[0] += first;
    result->SetElements(elements);

    return result;
}

Ciphertext<DCRTPoly> EvalCoeffsToSlots(const std::vector<std::vector<ConstPlaintext>>& A,
                                       ConstCiphertext<DCRTPoly> ctxt) {
    uint32_t slots = ctxt->GetSlots();

    auto cc    = ctxt->GetCryptoContext();
    uint32_t M = cc->GetCyclotomicOrder();
    uint32_t N = cc->GetRingDimension();

    int32_t levelBudget     = m_paramsEnc[CKKS_BOOT_PARAMS::LEVEL_BUDGET];
    int32_t layersCollapse  = m_paramsEnc[CKKS_BOOT_PARAMS::LAYERS_COLL];
    int32_t remCollapse     = m_paramsEnc[CKKS_BOOT_PARAMS::LAYERS_REM];
    int32_t numRotations    = m_paramsEnc[CKKS_BOOT_PARAMS::NUM_ROTATIONS];
    int32_t b               = m_paramsEnc[CKKS_BOOT_PARAMS::BABY_STEP];
    int32_t g               = m_paramsEnc[CKKS_BOOT_PARAMS::GIANT_STEP];
    int32_t numRotationsRem = m_paramsEnc[CKKS_BOOT_PARAMS::NUM_ROTATIONS_REM];
    int32_t bRem            = m_paramsEnc[CKKS_BOOT_PARAMS::BABY_STEP_REM];
    int32_t gRem            = m_paramsEnc[CKKS_BOOT_PARAMS::GIANT_STEP_REM];

    int32_t stop    = -1;
    int32_t flagRem = 0;

    auto algo = cc->GetScheme();

    if (remCollapse != 0) {
        stop    = 0;
        flagRem = 1;
    }

    // precompute the inner and outer rotations
    std::vector<std::vector<int32_t>> rot_in(levelBudget);
    for (uint32_t i = 0; i < uint32_t(levelBudget); i++) {
        if (flagRem == 1 && i == 0) {
            // remainder corresponds to index 0 in encoding and to last index in decoding
            rot_in[i] = std::vector<int32_t>(numRotationsRem + 1);
        }
        else {
            rot_in[i] = std::vector<int32_t>(numRotations + 1);
        }
    }

    std::vector<std::vector<int32_t>> rot_out(levelBudget);
    for (uint32_t i = 0; i < uint32_t(levelBudget); i++) {
        rot_out[i] = std::vector<int32_t>(b + bRem);
    }

    for (int32_t s = levelBudget - 1; s > stop; s--) {
        for (int32_t j = 0; j < g; j++) {
            rot_in[s][j] = ReduceRotation(
                (j - int32_t((numRotations + 1) / 2) + 1) * (1 << ((s - flagRem) * layersCollapse + remCollapse)),
                slots);
        }

        for (int32_t i = 0; i < b; i++) {
            rot_out[s][i] = ReduceRotation((g * i) * (1 << ((s - flagRem) * layersCollapse + remCollapse)), M / 4);
        }
    }

    if (flagRem) {
        for (int32_t j = 0; j < gRem; j++) {
            rot_in[stop][j] = ReduceRotation((j - int32_t((numRotationsRem + 1) / 2) + 1), slots);
        }

        for (int32_t i = 0; i < bRem; i++) {
            rot_out[stop][i] = ReduceRotation((gRem * i), M / 4);
        }
    }

    Ciphertext<DCRTPoly> result = ctxt->Clone();

    // hoisted automorphisms
    for (int32_t s = levelBudget - 1; s > stop; s--) {
        if (s != levelBudget - 1) {
            algo->ModReduceInternalInPlace(result, BASE_NUM_LEVELS_TO_DROP);
        }

        // computes the NTTs for each CRT limb (for the hoisted automorphisms used later on)
        auto digits = cc->EvalFastRotationPrecompute(result);

        std::vector<Ciphertext<DCRTPoly>> fastRotation(g);
#pragma omp parallel for
        for (int32_t j = 0; j < g; j++) {
            if (rot_in[s][j] != 0) {
                fastRotation[j] = cc->EvalFastRotationExt(result, rot_in[s][j], digits, true);
            }
            else {
                fastRotation[j] = cc->KeySwitchExt(result, true);
            }
        }

        Ciphertext<DCRTPoly> outer;
        DCRTPoly first;
        for (int32_t i = 0; i < b; i++) {
            // for the first iteration with j=0:
            int32_t G                  = g * i;
            Ciphertext<DCRTPoly> inner = EvalMultExt(fastRotation[0], A[s][G]);
            // continue the loop
            for (int32_t j = 1; j < g; j++) {
                if ((G + j) != int32_t(numRotations)) {
                    EvalAddExtInPlace(inner, EvalMultExt(fastRotation[j], A[s][G + j]));
                }
            }

            if (i == 0) {
                first         = cc->KeySwitchDownFirstElement(inner);
                auto elements = inner->GetElements();
                elements[0].SetValuesToZero();
                inner->SetElements(elements);
                outer = inner;
            }
            else {
                if (rot_out[s][i] != 0) {
                    inner = cc->KeySwitchDown(inner);
                    // Find the automorphism index that corresponds to rotation index index.
                    usint autoIndex = FindAutomorphismIndex2nComplex(rot_out[s][i], M);
                    std::vector<usint> map(N);
                    PrecomputeAutoMap(N, autoIndex, &map);
                    first += inner->GetElements()[0].AutomorphismTransform(autoIndex, map);
                    auto innerDigits = cc->EvalFastRotationPrecompute(inner);
                    EvalAddExtInPlace(outer, cc->EvalFastRotationExt(inner, rot_out[s][i], innerDigits, false));
                }
                else {
                    first += cc->KeySwitchDownFirstElement(inner);
                    auto elements = inner->GetElements();
                    elements[0].SetValuesToZero();
                    inner->SetElements(elements);
                    EvalAddExtInPlace(outer, inner);
                }
            }
        }
        result                          = cc->KeySwitchDown(outer);
        std::vector<DCRTPoly>& elements = result->GetElements();
        elements[0] += first;
    }

    if (flagRem) {
        algo->ModReduceInternalInPlace(result, BASE_NUM_LEVELS_TO_DROP);

        // computes the NTTs for each CRT limb (for the hoisted automorphisms used later on)
        auto digits = cc->EvalFastRotationPrecompute(result);
        std::vector<Ciphertext<DCRTPoly>> fastRotation(gRem);

#pragma omp parallel for
        for (int32_t j = 0; j < gRem; j++) {
            if (rot_in[stop][j] != 0) {
                fastRotation[j] = cc->EvalFastRotationExt(result, rot_in[stop][j], digits, true);
            }
            else {
                fastRotation[j] = cc->KeySwitchExt(result, true);
            }
        }

        Ciphertext<DCRTPoly> outer;
        DCRTPoly first;
        for (int32_t i = 0; i < bRem; i++) {
            Ciphertext<DCRTPoly> inner;
            // for the first iteration with j=0:
            int32_t GRem = gRem * i;
            inner        = EvalMultExt(fastRotation[0], A[stop][GRem]);
            // continue the loop
            for (int32_t j = 1; j < gRem; j++) {
                if ((GRem + j) != int32_t(numRotationsRem)) {
                    EvalAddExtInPlace(inner, EvalMultExt(fastRotation[j], A[stop][GRem + j]));
                }
            }

            if (i == 0) {
                first         = cc->KeySwitchDownFirstElement(inner);
                auto elements = inner->GetElements();
                elements[0].SetValuesToZero();
                inner->SetElements(elements);
                outer = inner;
            }
            else {
                if (rot_out[stop][i] != 0) {
                    inner = cc->KeySwitchDown(inner);
                    // Find the automorphism index that corresponds to rotation index index.
                    usint autoIndex = FindAutomorphismIndex2nComplex(rot_out[stop][i], M);
                    std::vector<usint> map(N);
                    PrecomputeAutoMap(N, autoIndex, &map);
                    first += inner->GetElements()[0].AutomorphismTransform(autoIndex, map);
                    auto innerDigits = cc->EvalFastRotationPrecompute(inner);
                    EvalAddExtInPlace(outer, cc->EvalFastRotationExt(inner, rot_out[stop][i], innerDigits, false));
                }
                else {
                    first += cc->KeySwitchDownFirstElement(inner);
                    auto elements = inner->GetElements();
                    elements[0].SetValuesToZero();
                    inner->SetElements(elements);
                    EvalAddExtInPlace(outer, inner);
                }
            }
        }

        result                          = cc->KeySwitchDown(outer);
        std::vector<DCRTPoly>& elements = result->GetElements();
        elements[0] += first;
    }

    return result;
}

Ciphertext<DCRTPoly> EvalSlotsToCoeffs(const std::vector<std::vector<ConstPlaintext>>& A,
                                       ConstCiphertext<DCRTPoly> ctxt) {
    auto cc = ctxt->GetCryptoContext();

    uint32_t M = cc->GetCyclotomicOrder();
    uint32_t N = cc->GetRingDimension();

    int32_t levelBudget     = m_paramsDec[CKKS_BOOT_PARAMS::LEVEL_BUDGET];
    int32_t layersCollapse  = m_paramsDec[CKKS_BOOT_PARAMS::LAYERS_COLL];
    int32_t remCollapse     = m_paramsDec[CKKS_BOOT_PARAMS::LAYERS_REM];
    int32_t numRotations    = m_paramsDec[CKKS_BOOT_PARAMS::NUM_ROTATIONS];
    int32_t b               = m_paramsDec[CKKS_BOOT_PARAMS::BABY_STEP];
    int32_t g               = m_paramsDec[CKKS_BOOT_PARAMS::GIANT_STEP];
    int32_t numRotationsRem = m_paramsDec[CKKS_BOOT_PARAMS::NUM_ROTATIONS_REM];
    int32_t bRem            = m_paramsDec[CKKS_BOOT_PARAMS::BABY_STEP_REM];
    int32_t gRem            = m_paramsDec[CKKS_BOOT_PARAMS::GIANT_STEP_REM];

    auto algo = cc->GetScheme();

    int32_t flagRem = 0;

    if (remCollapse != 0) {
        flagRem = 1;
    }

    // precompute the inner and outer rotations

    std::vector<std::vector<int32_t>> rot_in(levelBudget);
    for (uint32_t i = 0; i < uint32_t(levelBudget); i++) {
        if (flagRem == 1 && i == uint32_t(levelBudget - 1)) {
            // remainder corresponds to index 0 in encoding and to last index in decoding
            rot_in[i] = std::vector<int32_t>(numRotationsRem + 1);
        }
        else {
            rot_in[i] = std::vector<int32_t>(numRotations + 1);
        }
    }

    std::vector<std::vector<int32_t>> rot_out(levelBudget);
    for (uint32_t i = 0; i < uint32_t(levelBudget); i++) {
        rot_out[i] = std::vector<int32_t>(b + bRem);
    }

    for (int32_t s = 0; s < levelBudget - flagRem; s++) {
        for (int32_t j = 0; j < g; j++) {
            rot_in[s][j] =
                ReduceRotation((j - int32_t((numRotations + 1) / 2) + 1) * (1 << (s * layersCollapse)), M / 4);
        }

        for (int32_t i = 0; i < b; i++) {
            rot_out[s][i] = ReduceRotation((g * i) * (1 << (s * layersCollapse)), M / 4);
        }
    }

    if (flagRem) {
        int32_t s = levelBudget - flagRem;
        for (int32_t j = 0; j < gRem; j++) {
            rot_in[s][j] =
                ReduceRotation((j - int32_t((numRotationsRem + 1) / 2) + 1) * (1 << (s * layersCollapse)), M / 4);
        }

        for (int32_t i = 0; i < bRem; i++) {
            rot_out[s][i] = ReduceRotation((gRem * i) * (1 << (s * layersCollapse)), M / 4);
        }
    }

    //  No need for Encrypted Bit Reverse
    Ciphertext<DCRTPoly> result = ctxt->Clone();

    // hoisted automorphisms
    for (int32_t s = 0; s < levelBudget - flagRem; s++) {
        if (s != 0) {
            algo->ModReduceInternalInPlace(result, BASE_NUM_LEVELS_TO_DROP);
        }
        // computes the NTTs for each CRT limb (for the hoisted automorphisms used later on)
        auto digits = cc->EvalFastRotationPrecompute(result);

        std::vector<Ciphertext<DCRTPoly>> fastRotation(g);
#pragma omp parallel for
        for (int32_t j = 0; j < g; j++) {
            if (rot_in[s][j] != 0) {
                fastRotation[j] = cc->EvalFastRotationExt(result, rot_in[s][j], digits, true);
            }
            else {
                fastRotation[j] = cc->KeySwitchExt(result, true);
            }
        }

        Ciphertext<DCRTPoly> outer;
        DCRTPoly first;
        for (int32_t i = 0; i < b; i++) {
            Ciphertext<DCRTPoly> inner;
            // for the first iteration with j=0:
            int32_t G = g * i;
            inner     = EvalMultExt(fastRotation[0], A[s][G]);
            // continue the loop
            for (int32_t j = 1; j < g; j++) {
                if ((G + j) != int32_t(numRotations)) {
                    EvalAddExtInPlace(inner, EvalMultExt(fastRotation[j], A[s][G + j]));
                }
            }

            if (i == 0) {
                first         = cc->KeySwitchDownFirstElement(inner);
                auto elements = inner->GetElements();
                elements[0].SetValuesToZero();
                inner->SetElements(elements);
                outer = inner;
            }
            else {
                if (rot_out[s][i] != 0) {
                    inner = cc->KeySwitchDown(inner);
                    // Find the automorphism index that corresponds to rotation index index.
                    usint autoIndex = FindAutomorphismIndex2nComplex(rot_out[s][i], M);
                    std::vector<usint> map(N);
                    PrecomputeAutoMap(N, autoIndex, &map);
                    first += inner->GetElements()[0].AutomorphismTransform(autoIndex, map);
                    auto innerDigits = cc->EvalFastRotationPrecompute(inner);
                    EvalAddExtInPlace(outer, cc->EvalFastRotationExt(inner, rot_out[s][i], innerDigits, false));
                }
                else {
                    first += cc->KeySwitchDownFirstElement(inner);
                    auto elements = inner->GetElements();
                    elements[0].SetValuesToZero();
                    inner->SetElements(elements);
                    EvalAddExtInPlace(outer, inner);
                }
            }
        }

        result                          = cc->KeySwitchDown(outer);
        std::vector<DCRTPoly>& elements = result->GetElements();
        elements[0] += first;
    }

    if (flagRem) {
        algo->ModReduceInternalInPlace(result, BASE_NUM_LEVELS_TO_DROP);
        // computes the NTTs for each CRT limb (for the hoisted automorphisms used later on)
        auto digits = cc->EvalFastRotationPrecompute(result);
        std::vector<Ciphertext<DCRTPoly>> fastRotation(gRem);

        int32_t s = levelBudget - flagRem;
#pragma omp parallel for
        for (int32_t j = 0; j < gRem; j++) {
            if (rot_in[s][j] != 0) {
                fastRotation[j] = cc->EvalFastRotationExt(result, rot_in[s][j], digits, true);
            }
            else {
                fastRotation[j] = cc->KeySwitchExt(result, true);
            }
        }

        Ciphertext<DCRTPoly> outer;
        DCRTPoly first;
        for (int32_t i = 0; i < bRem; i++) {
            Ciphertext<DCRTPoly> inner;
            // for the first iteration with j=0:
            int32_t GRem = gRem * i;
            inner        = EvalMultExt(fastRotation[0], A[s][GRem]);
            // continue the loop
            for (int32_t j = 1; j < gRem; j++) {
                if ((GRem + j) != int32_t(numRotationsRem))
                    EvalAddExtInPlace(inner, EvalMultExt(fastRotation[j], A[s][GRem + j]));
            }

            if (i == 0) {
                first         = cc->KeySwitchDownFirstElement(inner);
                auto elements = inner->GetElements();
                elements[0].SetValuesToZero();
                inner->SetElements(elements);
                outer = inner;
            }
            else {
                if (rot_out[s][i] != 0) {
                    inner = cc->KeySwitchDown(inner);
                    // Find the automorphism index that corresponds to rotation index index.
                    usint autoIndex = FindAutomorphismIndex2nComplex(rot_out[s][i], M);
                    std::vector<usint> map(N);
                    PrecomputeAutoMap(N, autoIndex, &map);
                    first += inner->GetElements()[0].AutomorphismTransform(autoIndex, map);
                    auto innerDigits = cc->EvalFastRotationPrecompute(inner);
                    EvalAddExtInPlace(outer, cc->EvalFastRotationExt(inner, rot_out[s][i], innerDigits, false));
                }
                else {
                    first += cc->KeySwitchDownFirstElement(inner);
                    auto elements = inner->GetElements();
                    elements[0].SetValuesToZero();
                    inner->SetElements(elements);
                    EvalAddExtInPlace(outer, inner);
                }
            }
        }

        result                          = cc->KeySwitchDown(outer);
        std::vector<DCRTPoly>& elements = result->GetElements();
        elements[0] += first;
    }

    return result;
}

Ciphertext<DCRTPoly> EvalMultExt(ConstCiphertext<DCRTPoly> ciphertext, ConstPlaintext plaintext) {
    Ciphertext<DCRTPoly> result = ciphertext->Clone();
    std::vector<DCRTPoly>& cv   = result->GetElements();

    DCRTPoly pt = plaintext->GetElement<DCRTPoly>();
    pt.SetFormat(Format::EVALUATION);

    for (auto& c : cv) {
        c *= pt;
    }
    result->SetNoiseScaleDeg(result->GetNoiseScaleDeg() + plaintext->GetNoiseScaleDeg());
    result->SetScalingFactor(result->GetScalingFactor() * plaintext->GetScalingFactor());
    return result;
}

void EvalAddExtInPlace(Ciphertext<DCRTPoly>& ciphertext1, ConstCiphertext<DCRTPoly> ciphertext2) {
    std::vector<DCRTPoly>& cv1       = ciphertext1->GetElements();
    const std::vector<DCRTPoly>& cv2 = ciphertext2->GetElements();

    for (size_t i = 0; i < cv1.size(); ++i) {
        cv1[i] += cv2[i];
    }
}

Ciphertext<DCRTPoly> EvalAddExt(ConstCiphertext<DCRTPoly> ciphertext1, ConstCiphertext<DCRTPoly> ciphertext2) {
    Ciphertext<DCRTPoly> result = ciphertext1->Clone();
    EvalAddExtInPlace(result, ciphertext2);
    return result;
}

Ciphertext<DCRTPoly> Conjugate(ConstCiphertext<DCRTPoly> ciphertext) {
    const std::vector<DCRTPoly>& cv = ciphertext->GetElements();
    usint N                         = cv[0].GetRingDimension();

    std::vector<usint> vec(N);
    PrecomputeAutoMap(N, 2 * N - 1, &vec);

    auto algo = ciphertext->GetCryptoContext()->GetScheme();

    Ciphertext<DCRTPoly> result = ciphertext->Clone();

    algo->KeySwitchInPlace(result, m_conjKey);

    std::vector<DCRTPoly>& rcv = result->GetElements();

    rcv[0] = rcv[0].AutomorphismTransform(2 * N - 1, vec);
    rcv[1] = rcv[1].AutomorphismTransform(2 * N - 1, vec);

    return result;
}

void AdjustCiphertext(Ciphertext<DCRTPoly>& ciphertext, double correction) {
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(ciphertext->GetCryptoParameters());

    auto cc   = ciphertext->GetCryptoContext();
    auto algo = cc->GetScheme();

    if (cryptoParams->GetScalingTechnique() == FLEXIBLEAUTO || cryptoParams->GetScalingTechnique() == FLEXIBLEAUTOEXT) {
        uint32_t lvl       = cryptoParams->GetScalingTechnique() == FLEXIBLEAUTO ? 0 : 1;
        double targetSF    = cryptoParams->GetScalingFactorReal(lvl);
        double sourceSF    = ciphertext->GetScalingFactor();
        uint32_t numTowers = ciphertext->GetElements()[0].GetNumOfElements();
        double modToDrop = cryptoParams->GetElementParams()->GetParams()[numTowers - 1]->GetModulus().ConvertToDouble();

        // in the case of FLEXIBLEAUTO, we need to bring the ciphertext to the right scale using a
        // a scaling multiplication. Note the at currently FLEXIBLEAUTO is only supported for NATIVEINT = 64.
        // So the other branch is for future purposes (in case we decide to add add the FLEXIBLEAUTO support
        // for NATIVEINT = 128.
#if NATIVEINT != 128
        // Scaling down the message by a correction factor to emulate using a larger q0.
        // This step is needed so we could use a scaling factor of up to 2^59 with q9 ~= 2^60.
        double adjustmentFactor = (targetSF / sourceSF) * (modToDrop / sourceSF) * correction;
#else
        double adjustmentFactor = (targetSF / sourceSF) * (modToDrop / sourceSF);
#endif
        cc->EvalMultInPlace(ciphertext, adjustmentFactor);

        algo->ModReduceInternalInPlace(ciphertext, BASE_NUM_LEVELS_TO_DROP);
        ciphertext->SetScalingFactor(targetSF);
    }
    else {
#if NATIVEINT != 128
        // Scaling down the message by a correction factor to emulate using a larger q0.
        // This step is needed so we could use a scaling factor of up to 2^59 with q9 ~= 2^60.
        cc->EvalMultInPlace(ciphertext, correction);
        algo->ModReduceInternalInPlace(ciphertext, BASE_NUM_LEVELS_TO_DROP);
#endif
    }
}

Ciphertext<DCRTPoly> EvalFuncBT(ConstCiphertext<DCRTPoly> ciphertext, uint32_t digitBitSize, BigInteger initialScaling,
                                uint64_t postScaling, bool step) {
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(ciphertext->GetCryptoParameters());

    if (cryptoParams->GetKeySwitchTechnique() != HYBRID)
        OPENFHE_THROW("CKKS Bootstrapping is only supported for the Hybrid key switching method.");
#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
    if (cryptoParams->GetScalingTechnique() == FLEXIBLEAUTO || cryptoParams->GetScalingTechnique() == FLEXIBLEAUTOEXT)
        OPENFHE_THROW("128-bit CKKS Bootstrapping is supported for FIXEDMANUAL and FIXEDAUTO methods only.");
#endif

#ifdef BOOTSTRAPTIMING
    TimeVar t;
    double timeEncode(0.0);
    double timeModReduce(0.0);
    double timeDecode(0.0);
#endif

    auto cc     = ciphertext->GetCryptoContext();
    uint32_t M  = cc->GetCyclotomicOrder();
    uint32_t L0 = cryptoParams->GetElementParams()->GetParams().size();
    // auto initSizeQ = ciphertext->GetElements()[0].GetNumOfElements();

    uint32_t slots = ciphertext->GetSlots();
    size_t N       = cc->GetRingDimension();

    auto elementParamsRaised = *(cryptoParams->GetElementParams());

    // // For FLEXIBLEAUTOEXT the raised ciphertext does not include extra modulus
    // // as it is multiplied by auxiliary plaintext
    // if (cryptoParams->GetScalingTechnique() == FLEXIBLEAUTOEXT) {
    //     elementParamsRaised.PopLastParam();
    // }

    auto paramsQ = elementParamsRaised.GetParams();
    usint sizeQ  = paramsQ.size();

    std::vector<NativeInteger> moduli(sizeQ);
    std::vector<NativeInteger> roots(sizeQ);
    for (size_t i = 0; i < sizeQ; i++) {
        moduli[i] = paramsQ[i]->GetModulus();
        roots[i]  = paramsQ[i]->GetRootOfUnity();
    }
    auto elementParamsRaisedPtr = std::make_shared<ILDCRTParams<DCRTPoly::Integer>>(M, moduli, roots);

    // Andreea: we don't need the type of scaling and correction as in the standard CKKS bootstrapping
    /*
    NativeInteger q = elementParamsRaisedPtr->GetParams()[0]->GetModulus().ConvertToInt();
    double qDouble  = q.ConvertToDouble();
    const auto p = cryptoParams->GetPlaintextModulus();
    double powP  = pow(2, p);
    std::cerr << "q: " << qDouble << std::endl;
    std::cerr << std::setprecision(15) << "p: " << powP << std::endl;

    int32_t deg = std::round(std::log2(qDouble / powP));
    // #if NATIVEINT != 128
    //     if (deg > static_cast<int32_t>(m_correctionFactor)) {
    //         OPENFHE_THROW("Degree [" + std::to_string(deg) + "] must be less than or equal to the correction factor [" +
    //                       std::to_string(m_correctionFactor) + "].");
    //     }
    // #endif
    double post = std::pow(2, static_cast<double>(deg));
    double pre = 1. / post;
    uint64_t scalar = std::llround(post);
    */
    double correction = cryptoParams->GetScalingFactorReal(0) / initialScaling.ConvertToDouble();  // 1.0;
    std::cerr << "correction: " << correction << std::endl;

    //------------------------------------------------------------------------------
    // RAISING THE MODULUS
    //------------------------------------------------------------------------------

    // In FLEXIBLEAUTO, raising the ciphertext to a larger number
    // of towers is a bit more complex, because we need to adjust
    // it's scaling factor to the one that corresponds to the level
    // it's being raised to.
    // Increasing the modulus

    Ciphertext<DCRTPoly> raised = ciphertext->Clone();
    auto algo                   = cc->GetScheme();
    algo->ModReduceInternalInPlace(raised, raised->GetNoiseScaleDeg() - 1);

    // Andreea: if correction ~ 1, we should not do this adjustment and save a level
    if (std::llround(correction) != 1.0) {
        AdjustCiphertext(raised, correction);
    }
    auto ctxtDCRT = raised->GetElements();

    // We only use the level 0 ciphertext here. All other towers are automatically ignored to make
    // CKKS bootstrapping faster.
    for (size_t i = 0; i < ctxtDCRT.size(); i++) {
        DCRTPoly temp(elementParamsRaisedPtr, COEFFICIENT);
        ctxtDCRT[i].SetFormat(COEFFICIENT);
        temp = ctxtDCRT[i].GetElementAtIndex(0);
        temp.SetFormat(EVALUATION);
        ctxtDCRT[i] = temp;
    }

    raised->SetElements(ctxtDCRT);
    raised->SetLevel(L0 - ctxtDCRT[0].GetNumOfElements());

    // std::cerr << "after raising" << std::endl;

#ifdef BOOTSTRAPTIMING
    std::cerr << "\nNumber of levels at the beginning of bootstrapping: "
              << raised->GetElements()[0].GetNumOfElements() - 1 << std::endl;
#endif

    //------------------------------------------------------------------------------
    // SETTING PARAMETERS FOR APPROXIMATE MODULAR REDUCTION
    //------------------------------------------------------------------------------

    // Coefficients of the Chebyshev series interpolating 0.5(1-cos(2 Pi K x)) for mod 2
    std::vector<double> coefficients1;
    std::vector<double> coefficients2;
    double k = 0;

    if (cryptoParams->GetSecretKeyDist() == SPARSE_TERNARY) {
        switch (digitBitSize) {
            case 1:
                coefficients1 = coeff_cos_12_mod2;  // g_coefficientsSparse;
                break;
            case 2:
                coefficients1 = coeff_cos_12_mod4;
                coefficients2 = coeff_sin_12_mod4;
                break;
            case 3:
                coefficients1 = coeff_cos_12_mod4;
                coefficients2 = coeff_sin_12_mod4;
                break;
            default:
                OPENFHE_THROW("Digit sizes of more than 3 bits are not currently allowed.");
                break;
        }

        if (step) {
            std::transform(coefficients1.begin(), coefficients1.end(), coefficients1.begin(),
                           std::bind(std::multiplies<double>(), std::placeholders::_1, 1.0 / (1 << digitBitSize)));
            std::transform(coefficients2.begin(), coefficients2.end(), coefficients2.begin(),
                           std::bind(std::multiplies<double>(), std::placeholders::_1, 1.0 / (1 << digitBitSize)));
        }

        // k = K_SPARSE;
        k = 1.0;  // do not divide by k as we already did it during precomputation
    }
    else {
        coefficients1 = g_coefficientsUniform;
        k             = K_UNIFORM;
    }

    // double constantEvalMult = pre * (1.0 / (k * N));
    double constantEvalMult = 1.0 / (k * N);

    cc->EvalMultInPlace(raised, constantEvalMult);

    // no linear transformations are needed for Chebyshev series as the range has been normalized to [-1,1]
    double coeffLowerBound = -1;
    double coeffUpperBound = 1;

    Ciphertext<DCRTPoly> ctxtDec;

    bool isLTBootstrap = (m_levelEnc == 1) && (m_levelDec == 1);
    if (slots == M / 4) {
        //------------------------------------------------------------------------------
        // FULLY PACKED CASE
        //------------------------------------------------------------------------------

#ifdef BOOTSTRAPTIMING
        TIC(t);
#endif

        //------------------------------------------------------------------------------
        // Running CoeffToSlot
        //------------------------------------------------------------------------------

        // need to call internal modular reduction so it also works for FLEXIBLEAUTO
        algo->ModReduceInternalInPlace(raised, BASE_NUM_LEVELS_TO_DROP);

        // only one linear transform is needed as the other one can be derived
        auto ctxtEnc =
            (isLTBootstrap) ? EvalLinearTransform(m_U0hatTPre, raised) : EvalCoeffsToSlots(m_U0hatTPreFFT, raised);

        // auto evalKeyMap = cc->GetEvalAutomorphismKeyMap(ctxtEnc->GetKeyTag());
        auto conj = Conjugate(ctxtEnc);

        auto ctxtEncI = cc->EvalSub(ctxtEnc, conj);
        cc->EvalAddInPlace(ctxtEnc, conj);
        algo->MultByMonomialInPlace(ctxtEncI, 3 * M / 4);

        if (cryptoParams->GetScalingTechnique() == FIXEDMANUAL) {
            while (ctxtEnc->GetNoiseScaleDeg() > 1) {
                cc->ModReduceInPlace(ctxtEnc);
                cc->ModReduceInPlace(ctxtEncI);
            }
        }
        else {
            if (ctxtEnc->GetNoiseScaleDeg() == 2) {
                algo->ModReduceInternalInPlace(ctxtEnc, BASE_NUM_LEVELS_TO_DROP);
                algo->ModReduceInternalInPlace(ctxtEncI, BASE_NUM_LEVELS_TO_DROP);
            }
        }

        // std::cerr << "Done with coeff to slot full" << std::endl;

        //------------------------------------------------------------------------------
        // Running Approximate Mod Reduction
        //------------------------------------------------------------------------------

        if (digitBitSize == 1) {  // Andreea: repeating code atm to not make extra copies of ciphertexts
            // Evaluate Chebyshev series for the cosine wave
            ctxtEnc  = cc->EvalChebyshevSeries(ctxtEnc, coefficients1, coeffLowerBound, coeffUpperBound);
            ctxtEncI = cc->EvalChebyshevSeries(ctxtEncI, coefficients1, coeffLowerBound, coeffUpperBound);

            // Double-angle iterations
            if ((cryptoParams->GetSecretKeyDist() == UNIFORM_TERNARY) ||
                (cryptoParams->GetSecretKeyDist() == SPARSE_TERNARY)) {
                if (cryptoParams->GetScalingTechnique() != FIXEDMANUAL) {
                    algo->ModReduceInternalInPlace(ctxtEnc, BASE_NUM_LEVELS_TO_DROP);
                    algo->ModReduceInternalInPlace(ctxtEncI, BASE_NUM_LEVELS_TO_DROP);
                }
                uint32_t numIter;
                if (cryptoParams->GetSecretKeyDist() == UNIFORM_TERNARY)
                    numIter = R_UNIFORM;
                else
                    numIter = R_SPARSE;
                ApplyDoubleAngleIterations(ctxtEnc, numIter);
                ApplyDoubleAngleIterations(ctxtEncI, numIter);
            }
        }

        if (digitBitSize == 2) {
            // Evaluate Chebyshev series for the sine wave
            auto ctxtEnc2  = cc->EvalChebyshevSeries(ctxtEnc, coefficients2, coeffLowerBound, coeffUpperBound);
            auto ctxtEncI2 = cc->EvalChebyshevSeries(ctxtEncI, coefficients2, coeffLowerBound, coeffUpperBound);

            // Double-angle iterations
            if ((cryptoParams->GetSecretKeyDist() == UNIFORM_TERNARY) ||
                (cryptoParams->GetSecretKeyDist() == SPARSE_TERNARY)) {
                if (cryptoParams->GetScalingTechnique() != FIXEDMANUAL) {
                    algo->ModReduceInternalInPlace(ctxtEnc2, BASE_NUM_LEVELS_TO_DROP);
                    algo->ModReduceInternalInPlace(ctxtEncI2, BASE_NUM_LEVELS_TO_DROP);
                }
                uint32_t numIter;
                if (cryptoParams->GetSecretKeyDist() == UNIFORM_TERNARY)
                    numIter = R_UNIFORM;
                else
                    numIter = R_SPARSE;
                ApplyDoubleAngleIterations(ctxtEnc2, numIter);
                ApplyDoubleAngleIterations(ctxtEncI2, numIter);
            }

            // Evaluate Chebyshev series for the cosine wave
            ctxtEnc  = cc->EvalChebyshevSeries(ctxtEnc, coefficients1, coeffLowerBound, coeffUpperBound);
            ctxtEncI = cc->EvalChebyshevSeries(ctxtEncI, coefficients1, coeffLowerBound, coeffUpperBound);

            // Double-angle iterations
            if ((cryptoParams->GetSecretKeyDist() == UNIFORM_TERNARY) ||
                (cryptoParams->GetSecretKeyDist() == SPARSE_TERNARY)) {
                if (cryptoParams->GetScalingTechnique() != FIXEDMANUAL) {
                    algo->ModReduceInternalInPlace(ctxtEnc, BASE_NUM_LEVELS_TO_DROP);
                    algo->ModReduceInternalInPlace(ctxtEncI, BASE_NUM_LEVELS_TO_DROP);
                }
                uint32_t numIter;
                if (cryptoParams->GetSecretKeyDist() == UNIFORM_TERNARY)
                    numIter = R_UNIFORM;
                else
                    numIter = R_SPARSE;
                ApplyDoubleAngleIterations(ctxtEnc, numIter);
                ApplyDoubleAngleIterations(ctxtEncI, numIter);
            }

            // Post-process cos(2pi x) and sin(2pi x) to get the mod 4 approximation or the step 4 approximation
            if (!step) {
                auto result = cc->EvalAdd(cc->EvalSub(ctxtEnc, ctxtEnc2), 1.0);
                cc->EvalSquareInPlace(ctxtEnc);
                cc->ModReduceInPlace(ctxtEnc);
                result = cc->EvalMult(result, ctxtEnc);
                cc->ModReduceInPlace(result);
                ctxtEnc2 = cc->EvalSub(2.0, ctxtEnc2);
                ctxtEnc  = cc->EvalSub(ctxtEnc2, result);

                result = cc->EvalAdd(cc->EvalSub(ctxtEncI, ctxtEncI2), 1.0);
                cc->EvalSquareInPlace(ctxtEncI);
                cc->ModReduceInPlace(ctxtEncI);
                result = cc->EvalMult(result, ctxtEncI);
                cc->ModReduceInPlace(result);
                ctxtEncI2 = cc->EvalSub(2.0, ctxtEncI2);
                ctxtEncI  = cc->EvalSub(ctxtEncI2, result);
            }
            else {
                auto result = cc->EvalAdd(ctxtEnc, ctxtEnc2);
                cc->EvalSquareInPlace(ctxtEnc);
                algo->MultByIntegerInPlace(ctxtEnc, 4);
                cc->ModReduceInPlace(ctxtEnc);
                result = cc->EvalMult(result, ctxtEnc);
                cc->ModReduceInPlace(result);
                ctxtEnc2 = cc->EvalSub(0.5, ctxtEnc2);
                ctxtEnc  = cc->EvalSub(ctxtEnc2, result);

                result = cc->EvalAdd(ctxtEncI, ctxtEncI2);
                cc->EvalSquareInPlace(ctxtEncI);
                algo->MultByIntegerInPlace(ctxtEncI, 4);
                cc->ModReduceInPlace(ctxtEncI);
                result = cc->EvalMult(result, ctxtEncI);
                cc->ModReduceInPlace(result);
                ctxtEnc2 = cc->EvalSub(0.5, ctxtEncI2);
                ctxtEncI = cc->EvalSub(ctxtEncI2, result);
            }
        }

        algo->MultByMonomialInPlace(ctxtEncI, M / 4);
        cc->EvalAddInPlace(ctxtEnc, ctxtEncI);

        // scale the message back up after Chebyshev interpolation
        // algo->MultByIntegerInPlace(ctxtEnc, scalar); // Andreea: no need, scalar = 1.0

#ifdef BOOTSTRAPTIMING
        timeModReduce = TOC(t);

        std::cerr << "Approximate modular reduction time: " << timeModReduce / 1000.0 << " s" << std::endl;

        // Running SlotToCoeff

        TIC(t);
#endif

        // std::cerr << "Done with modular approx full" << std::endl;

        //------------------------------------------------------------------------------
        // Running SlotToCoeff
        //------------------------------------------------------------------------------

        // In the case of FLEXIBLEAUTO, we need one extra tower
        // TODO: See if we can remove the extra level in FLEXIBLEAUTO
        if (cryptoParams->GetScalingTechnique() != FIXEDMANUAL) {
            algo->ModReduceInternalInPlace(ctxtEnc, BASE_NUM_LEVELS_TO_DROP);
        }

        // Only one linear transform is needed
        ctxtDec = (isLTBootstrap) ? EvalLinearTransform(m_U0Pre, ctxtEnc) : EvalSlotsToCoeffs(m_U0PreFFT, ctxtEnc);

        // Because the linear transform might be scaled differently, we might need to scale up the result separately
        if (postScaling > 1) {
            algo->MultByIntegerInPlace(ctxtDec, postScaling);
        }

        std::cerr << "ctxtDec levels: " << ctxtDec->GetLevel() << " and depth: " << ctxtDec->GetNoiseScaleDeg()
                  << std::endl;
        // std::cerr << "Done with slot to coeff full" << std::endl;
    }
    else {
        //------------------------------------------------------------------------------
        // SPARSELY PACKED CASE
        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------
        // Running PartialSum
        //------------------------------------------------------------------------------

        for (uint32_t j = 1; j < N / (2 * slots); j <<= 1) {
            auto temp = cc->EvalRotate(raised, j * slots);
            cc->EvalAddInPlace(raised, temp);
        }
        std::cerr << "after partial sum sparse" << std::endl;

#ifdef BOOTSTRAPTIMING
        TIC(t);
#endif

        //------------------------------------------------------------------------------
        // Running CoeffsToSlots
        //------------------------------------------------------------------------------

        algo->ModReduceInternalInPlace(raised, BASE_NUM_LEVELS_TO_DROP);

        auto ctxtEnc =
            (isLTBootstrap) ? EvalLinearTransform(m_U0hatTPre, raised) : EvalCoeffsToSlots(m_U0hatTPreFFT, raised);

        auto evalKeyMap = cc->GetEvalAutomorphismKeyMap(ctxtEnc->GetKeyTag());
        auto conj       = Conjugate(ctxtEnc);
        cc->EvalAddInPlace(ctxtEnc, conj);

        if (cryptoParams->GetScalingTechnique() == FIXEDMANUAL) {
            while (ctxtEnc->GetNoiseScaleDeg() > 1) {
                cc->ModReduceInPlace(ctxtEnc);
            }
        }
        else {
            if (ctxtEnc->GetNoiseScaleDeg() == 2) {
                algo->ModReduceInternalInPlace(ctxtEnc, BASE_NUM_LEVELS_TO_DROP);
            }
        }

#ifdef BOOTSTRAPTIMING
        timeEncode = TOC(t);

        std::cerr << "\nEncoding time: " << timeEncode / 1000.0 << " s" << std::endl;

        // Running Approximate Mod Reduction

        TIC(t);
#endif

        std::cerr << "after coefftoslots sparse" << std::endl;

        //------------------------------------------------------------------------------
        // Running Approximate Mod Reduction
        //------------------------------------------------------------------------------

        if (digitBitSize == 1) {  // Andreea: repeating code atm to not make extra copies of ciphertexts
            // Evaluate Chebyshev series for the cosine wave
            ctxtEnc = cc->EvalChebyshevSeries(ctxtEnc, coefficients1, coeffLowerBound, coeffUpperBound);

            // Double-angle iterations
            if ((cryptoParams->GetSecretKeyDist() == UNIFORM_TERNARY) ||
                (cryptoParams->GetSecretKeyDist() == SPARSE_TERNARY)) {
                if (cryptoParams->GetScalingTechnique() != FIXEDMANUAL) {
                    algo->ModReduceInternalInPlace(ctxtEnc, BASE_NUM_LEVELS_TO_DROP);
                }
                uint32_t numIter;
                if (cryptoParams->GetSecretKeyDist() == UNIFORM_TERNARY)
                    numIter = R_UNIFORM;
                else
                    numIter = R_SPARSE;
                ApplyDoubleAngleIterations(ctxtEnc, numIter);
            }
        }

        if (digitBitSize == 2) {
            // Evaluate Chebyshev series for the sine wave
            auto ctxtEnc2 = cc->EvalChebyshevSeries(ctxtEnc, coefficients2, coeffLowerBound, coeffUpperBound);

            // Double-angle iterations
            if ((cryptoParams->GetSecretKeyDist() == UNIFORM_TERNARY) ||
                (cryptoParams->GetSecretKeyDist() == SPARSE_TERNARY)) {
                if (cryptoParams->GetScalingTechnique() != FIXEDMANUAL) {
                    algo->ModReduceInternalInPlace(ctxtEnc, BASE_NUM_LEVELS_TO_DROP);
                }
                uint32_t numIter;
                if (cryptoParams->GetSecretKeyDist() == UNIFORM_TERNARY)
                    numIter = R_UNIFORM;
                else
                    numIter = R_SPARSE;
                ApplyDoubleAngleIterations(ctxtEnc2, numIter);
            }

            // Evaluate Chebyshev series for the cosine wave
            ctxtEnc = cc->EvalChebyshevSeries(ctxtEnc, coefficients1, coeffLowerBound, coeffUpperBound);

            // Double-angle iterations
            if ((cryptoParams->GetSecretKeyDist() == UNIFORM_TERNARY) ||
                (cryptoParams->GetSecretKeyDist() == SPARSE_TERNARY)) {
                if (cryptoParams->GetScalingTechnique() != FIXEDMANUAL) {
                    algo->ModReduceInternalInPlace(ctxtEnc, BASE_NUM_LEVELS_TO_DROP);
                }
                uint32_t numIter;
                if (cryptoParams->GetSecretKeyDist() == UNIFORM_TERNARY)
                    numIter = R_UNIFORM;
                else
                    numIter = R_SPARSE;
                ApplyDoubleAngleIterations(ctxtEnc, numIter);
            }

            // Post-process cos(2pi x) and sin(2pi x) to get the mod 4 approximation or the step 4 approximation
            if (!step) {
                auto result = cc->EvalAdd(cc->EvalSub(ctxtEnc, ctxtEnc2), 1.0);
                cc->EvalSquareInPlace(ctxtEnc);
                cc->ModReduceInPlace(ctxtEnc);
                result = cc->EvalMult(result, ctxtEnc);
                cc->ModReduceInPlace(result);
                ctxtEnc2 = cc->EvalSub(2.0, ctxtEnc2);
                ctxtEnc  = cc->EvalSub(ctxtEnc2, result);
            }
            else {
                auto result = cc->EvalAdd(ctxtEnc, ctxtEnc2);
                cc->EvalSquareInPlace(ctxtEnc);
                algo->MultByIntegerInPlace(ctxtEnc, 4);
                cc->ModReduceInPlace(ctxtEnc);
                result = cc->EvalMult(result, ctxtEnc);
                cc->ModReduceInPlace(result);
                ctxtEnc2 = cc->EvalSub(0.5, ctxtEnc2);
                ctxtEnc  = cc->EvalSub(ctxtEnc2, result);
            }
        }

        // scale the message back up after Chebyshev interpolation
        // algo->MultByIntegerInPlace(ctxtEnc, scalar); // Andreea: no need, scalar = 1.0

#ifdef BOOTSTRAPTIMING
        timeModReduce = TOC(t);

        std::cerr << "Approximate modular reduction time: " << timeModReduce / 1000.0 << " s" << std::endl;

        // Running SlotToCoeff

        TIC(t);
#endif
        std::cerr << "after approx mod sparse" << std::endl;

        //------------------------------------------------------------------------------
        // Running SlotsToCoeffs
        //------------------------------------------------------------------------------

        // In the case of FLEXIBLEAUTO, we need one extra tower
        // TODO: See if we can remove the extra level in FLEXIBLEAUTO
        if (cryptoParams->GetScalingTechnique() != FIXEDMANUAL) {
            algo->ModReduceInternalInPlace(ctxtEnc, BASE_NUM_LEVELS_TO_DROP);
        }

        // linear transform for decoding
        ctxtDec = (isLTBootstrap) ? EvalLinearTransform(m_U0Pre, ctxtEnc) : EvalSlotsToCoeffs(m_U0PreFFT, ctxtEnc);

        cc->EvalAddInPlace(ctxtDec, cc->EvalRotate(ctxtDec, slots));

        // Because the linear transform might be scaled differently, we might need to scale up the result separately
        if (postScaling > 1) {
            algo->MultByIntegerInPlace(ctxtDec, postScaling);
        }

        std::cerr << "after slotstocoeff sparse" << std::endl;
    }

#if NATIVEINT != 128
    // Andreea: this is not necessary
    // // 64-bit only: scale back the message to its original scale.
    // // uint64_t corFactor = (uint64_t)1 << std::llround(correction);
    // algo->MultByIntegerInPlace(ctxtDec, 512);
#endif

#ifdef BOOTSTRAPTIMING
    timeDecode = TOC(t);

    std::cout << "Decoding time: " << timeDecode / 1000.0 << " s" << std::endl;
#endif

    return ctxtDec;
}

void ApplyDoubleAngleIterations(Ciphertext<DCRTPoly>& ciphertext, uint32_t numIter) {
    auto cc = ciphertext->GetCryptoContext();

    int32_t r = numIter;
    for (int32_t j = 1; j < r + 1; j++) {
        cc->EvalSquareInPlace(ciphertext);
        ciphertext    = cc->EvalAdd(ciphertext, ciphertext);
        double scalar = -1.0 / std::pow((2.0 * M_PI), std::pow(2.0, j - r));
        cc->EvalAddInPlace(ciphertext, scalar);
        cc->ModReduceInPlace(ciphertext);
    }
}

//------------------------------------------------------------------------------
// Main code
//------------------------------------------------------------------------------

void Floor() {
    std::cerr << "\n=======Floor evaluation=======\n";

    // BFV parameters
    BigInteger Q("1152921504606846976");  // 2^60
    BigInteger p("1048576");              // 2^20
    // BigInteger Bigq          = BigInteger("35184372088832"); // Mod 2^45
    // BigInteger Bigq = BigInteger("2199023255552");  // Mod 2^41
    // BigInteger pNew("2");                           // 2
    BigInteger Bigq = BigInteger("4398046511104");  // Mod 2^42
    BigInteger pNew("4");                           // 4

    uint32_t dcrtBits = Bigq.GetMSB() - 1;
    uint32_t firstMod = Bigq.GetMSB() - 1;
    uint32_t numSlots = 8;

    CCParams<CryptoContextCKKSRNS> parameters;
    SecretKeyDist secretKeyDist = SPARSE_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(FIXEDMANUAL);
    parameters.SetFirstModSize(firstMod);
    parameters.SetNumLargeDigits(3);
    parameters.SetBatchSize(numSlots);
    parameters.SetRingDim(32);

    std::vector<uint32_t> levelBudget = {1, 1};

    // Andreea: setting levelsAvailableAfterBootstrap different than 2 creates issues.
    // When 1, fresh CKKS ciphertexts (no modulus reduction is being done) evaluate to zeroes
    // When 3, the floor evaluation from BFV evaluates to zeroes, but fresh CKKS ciphertexts evaluate correctly
    uint32_t levelsAvailableAfterBootstrap = 2;
    usint depth                            = levelsAvailableAfterBootstrap + 8 + pNew.GetMSB() - 1;
    depth                                  = (pNew.GetMSB() > 2) ? depth + 1 : depth;
    parameters.SetMultiplicativeDepth(depth);

    CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);

    cryptoContext->Enable(PKE);
    cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);
    cryptoContext->Enable(ADVANCEDSHE);
    cryptoContext->Enable(FHE);

    usint ringDim = cryptoContext->GetRingDimension();
    std::cout << "CKKS scheme is using ring dimension " << ringDim << ", depth " << depth << ", firstModSize "
              << firstMod << " and scalingModSize " << dcrtBits << std::endl
              << std::endl;

    auto keyPair = cryptoContext->KeyGen();
    cryptoContext->EvalMultKeyGen(keyPair.secretKey);

    // More CKKS parameters
    // The scaling done in decoding needs to align the CKKS ciphertext after bootstrapping with the BFV ciphertext
    BigInteger QPrime = keyPair.publicKey->GetPublicElements()[0].GetParams()->GetParams()[0]->GetModulus();
    uint32_t cnt      = 1;
    auto levels       = levelsAvailableAfterBootstrap;
    while (levels > 1) {
        QPrime *= keyPair.publicKey->GetPublicElements()[0].GetParams()->GetParams()[cnt]->GetModulus();
        levels--;
        cnt++;
    }
    std::cerr << "QPrime: " << QPrime << std::endl;
    double scaleMod = QPrime.ConvertToDouble() / (Bigq.ConvertToDouble() * p.ConvertToDouble());
    // double scaleMod = 1.0; // For fresh CKKS ciphertexts, not originating from BFV. // Andreea: if we start on level 0, this gets problematic
    std::cerr << "scaleMod = " << scaleMod << std::endl;
    EvalFuncBTSetup(*cryptoContext, numSlots, pNew.GetMSB() - 1, {0, 0}, levelBudget, scaleMod);
    EvalFuncBTKeyGen(keyPair.secretKey, numSlots);

    // =======Step 1. Create a CKKS ciphertext with proper metadata=======
    // We create the ciphertext modulus to be used in BFV

    // std::vector<double> y = {0.0, 1.0, 0.0, 1.0};
    // std::vector<std::complex<double>> y = { std::exp(Pi*1i/8.0) + std::exp(Pi*3i/8.0), std::exp(Pi*5i/8.0) + std::exp(Pi*15i/8.0), std::exp(Pi*9i/8.0) + std::exp(Pi*27i/8.0), std::exp(Pi*13i/8.0) + std::exp(Pi*39i/8.0)};
    // std::vector<double> y = {0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0};
    // std::vector<std::complex<double>> y = {
    //     std::exp(Pi * 1i / 16.0) + std::exp(Pi * 3i / 16.0),   std::exp(Pi * 5i / 16.0) + std::exp(Pi * 15i / 16.0),
    //     std::exp(Pi * 25i / 16.0) + std::exp(Pi * 75i / 16.0), std::exp(Pi * 29i / 16.0) + std::exp(Pi * 87i / 16.0),
    //     std::exp(Pi * 17i / 16.0) + std::exp(Pi * 51i / 16.0), std::exp(Pi * 21i / 16.0) + std::exp(Pi * 63i / 16.0),
    //     std::exp(Pi * 9i / 16.0) + std::exp(Pi * 27i / 16.0),  std::exp(Pi * 13i / 16.0) + std::exp(Pi * 39i / 16.0)};
    std::vector<std::complex<double>> y = {
        std::exp(Pi * 1i / 16.0) + 2.0 * std::exp(Pi * 2i / 16.0) - std::exp(Pi * 3i / 16.0) +
            std::exp(Pi * 5i / 16.0) - std::exp(Pi * 7i / 16.0),
        std::exp(Pi * 5i / 16.0) + 2.0 * std::exp(Pi * 10i / 16.0) - std::exp(Pi * 15i / 16.0) +
            std::exp(Pi * 25i / 16.0) - std::exp(Pi * 35i / 16.0),
        std::exp(Pi * 25i / 16.0) + 2.0 * std::exp(Pi * 2 * 25i / 16.0) - std::exp(Pi * 3 * 25i / 16.0) +
            std::exp(Pi * 5 * 25i / 16.0) - std::exp(Pi * 7 * 25i / 16.0),
        std::exp(Pi * 29i / 16.0) + 2.0 * std::exp(Pi * 2 * 29i / 16.0) - std::exp(Pi * 3 * 29i / 16.0) +
            std::exp(Pi * 5 * 29i / 16.0) - std::exp(Pi * 7 * 29i / 16.0),
        std::exp(Pi * 17i / 16.0) + 2.0 * std::exp(Pi * 2 * 17i / 16.0) - std::exp(Pi * 3 * 17i / 16.0) +
            std::exp(Pi * 5 * 17i / 16.0) - std::exp(Pi * 7 * 17i / 16.0),
        std::exp(Pi * 21i / 16.0) + 2.0 * std::exp(Pi * 2 * 21i / 16.0) - std::exp(Pi * 3 * 21i / 16.0) +
            std::exp(Pi * 5 * 21i / 16.0) - std::exp(Pi * 7 * 21i / 16.0),
        std::exp(Pi * 9i / 16.0) + 2.0 * std::exp(Pi * 2 * 9i / 16.0) - std::exp(Pi * 3 * 9i / 16.0) +
            std::exp(Pi * 5 * 9i / 16.0) - std::exp(Pi * 7 * 9i / 16.0),
        std::exp(Pi * 13i / 16.0) + 2.0 * std::exp(Pi * 2 * 13i / 16.0) - std::exp(Pi * 3 * 13i / 16.0) +
            std::exp(Pi * 5 * 13i / 16.0) - std::exp(Pi * 7 * 13i / 16.0)};
    std::transform(
        y.begin(), y.end(), y.begin(),
        std::bind(std::multiplies<std::complex<double>>(), std::placeholders::_1, 1.0 / pNew.ConvertToDouble()));
    std::cerr << "y = " << y << std::endl;
    // depth - 1 means we have two RNS limbs here; we need to the second limb
    // for internal downscaling (scalar multiplication)
    // so that the sine wave approximation of modular reduction
    // could achieve reasonable precision. Is this true?
    Plaintext ptxt = cryptoContext->MakeCKKSPackedPlaintext(y, 1, depth - levelsAvailableAfterBootstrap + 1);
    ptxt->SetLength(numSlots);
    Ciphertext<DCRTPoly> ctxt = cryptoContext->Encrypt(keyPair.publicKey, ptxt);

    // Switch encryped digit from q to q'
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(ctxt->GetCryptoParameters());
    // auto elementParams = cryptoParams->GetElementParams(); // Andreea: this still takes the full chain of moduli
    auto elementParams = ptxt->GetElement<DCRTPoly>().GetParams();

    auto qPrimeCKKS = elementParams->GetModulus();
    std::cerr << "\nqPrimeCKKS: " << qPrimeCKKS << std::endl << std::endl;

    // =======Step 2. Encrypting and decrypting using BFV-like encryption=======
    std::cerr << "\nBigQBFV: " << Q << std::endl << std::endl;

    std::vector<int64_t> x = {524288, 1, 2, 3, 16, 33, 64, 1048575};

    std::cerr << "Plaintext before BFV encryption: " << x << std::endl;

    // // This creates BFV ciphertext with the full CKKS modulus - we don't need that
    // auto ctxtBFV = EncryptBFVCoeff(x, Q, p, keyPair.secretKey);
    // auto decrypted = DecryptBFVCoeff(encrypted, Q, p, keyPair.secretKey, numSlots);

    // This creates BFV ciphertext with the ciphertext modulus corresponding to the last two limbs
    auto ctxtBFV   = EncryptBFVCoeff(x, Q, p, keyPair.secretKey, elementParams);
    auto decrypted = DecryptBFVCoeff(ctxtBFV, Q, p, keyPair.secretKey, elementParams, numSlots);
    std::cerr << "Plaintext after BFV encryption + decryption: " << decrypted << std::endl;

    // =======Step 3. Changing (\log Q, \log p) from  (41,60) to (41,2), i.e., doing mod q=======

    auto encryptedDigit = ctxtBFV;

    // Apply mod q
    encryptedDigit[0].SwitchModulus(Bigq, 1, 0, 0);
    encryptedDigit[1].SwitchModulus(Bigq, 1, 0, 0);

    decrypted = DecryptBFVCoeff(encryptedDigit, Bigq, pNew, keyPair.secretKey, elementParams, numSlots);
    std::cerr << "Plaintext after BFV decryption of ciphertext mod q: " << decrypted << std::endl;

    // =======Step 4. Populate CKKS ciphertext from the BFV ciphertext=======
    std::vector<lbcrypto::DCRTPoly> elementsCKKS(2);
    if (qPrimeCKKS > Bigq) {
        elementsCKKS = ModSwitchUp(encryptedDigit, Bigq, qPrimeCKKS, elementParams);
    }
    else {
        elementsCKKS = ModSwitchDown(encryptedDigit, Bigq, qPrimeCKKS, elementParams);
    }
    auto ctxtNew = ctxt->Clone();
    ctxtNew->SetElements(elementsCKKS);
    // ctxt->SetElements(elementsCKKS);
    std::cerr << "Number of elements of populated ciphertext: " << ctxtNew->GetElements()[0].GetNumOfElements()
              << std::endl;
    std::cerr << "Populated CKKS ciphertext levels: " << ctxtNew->GetLevel()
              << " and depth: " << ctxtNew->GetNoiseScaleDeg() << std::endl;
    std::cerr << "Number of elements of initial ciphertext: " << ctxt->GetElements()[0].GetNumOfElements() << std::endl;
    std::cerr << "Initial CKKS ciphertext levels: " << ctxt->GetLevel() << " and depth: " << ctxt->GetNoiseScaleDeg()
              << std::endl;

    // =======Check decryption correctness=======

    Plaintext result;
    // cryptoContext->Decrypt(keyPair.secretKey, ctxt, &result);
    // result->SetLength(numSlots);
    // std::cout << "\n Decrypting CKKS initial ciphertext" << result << std::endl;
    // std::cout << "Decrypting BFV ciphertext switched to CKKS GetElement<Poly> = " << result->GetElement<Poly>() << std::endl;

    // auto initial = DecryptCKKSCoeff({ctxt->GetElements()[0], ctxt->GetElements()[1]}, keyPair.secretKey, numSlots);
    // std::cerr << "Decrypting CKKS initial ciphertext: " << initial << std::endl;

    auto output = DecryptCKKSCoeff(elementsCKKS, keyPair.secretKey, numSlots);
    std::cerr << "Decrypting BFV ciphertext switched to CKKS: " << output << std::endl;

    auto vec1 = DecryptWithoutDecode(*cryptoContext, ctxtNew, keyPair.secretKey, numSlots,
                                     cryptoContext->GetRingDimension(), false);
    std::cerr << "\nDecrypt without decoding switched to CKKS: " << std::setprecision(15) << vec1 << std::endl;

    std::cout << "\nInitial number of levels remaining: " << depth - ctxtNew->GetLevel() << std::endl;

    // =======Step 5. Bootstrap the digit=======

    auto ctxtAfterFuncBT = EvalFuncBT(ctxtNew, pNew.GetMSB() - 1, qPrimeCKKS, 1.0);
    // auto ctxtAfterFuncBT = EvalFuncBT(ctxt, pNew.GetMSB() - 1, Bigq, 1.0); // Test for fresh ciphertext

    cryptoContext->ModReduceInPlace(ctxtAfterFuncBT);

    auto vec2 = DecryptWithoutDecode(*cryptoContext, ctxtAfterFuncBT, keyPair.secretKey, numSlots,
                                     cryptoContext->GetRingDimension(), true);

    std::cerr << "\nDecrypted func bootstrapped digit without decoding: " << vec2 << std::endl;

    elementsCKKS[0] = ctxtAfterFuncBT->GetElements()[0];
    elementsCKKS[1] = ctxtAfterFuncBT->GetElements()[1];

    output = DecryptCKKSCoeff(elementsCKKS, keyPair.secretKey, numSlots);
    std::cerr << "Decrypting func bootstrapped digit without scaling: " << output << std::endl;
    std::transform(output.begin(), output.end(), output.begin(),
                   std::bind(std::multiplies<double>(), std::placeholders::_1,
                             1.0 / cryptoParams->GetScalingFactorReal(0)));  // 1.0 / qPrimeCKKS.ConvertToDouble()));
    std::cerr << "Decrypting func bootstrapped digit scaled: " << output << std::endl;

    cryptoContext->Decrypt(keyPair.secretKey, ctxtAfterFuncBT, &result);
    result->SetLength(numSlots);
    std::cout << "\nFull decryption func bootstrapped digit " << result << std::endl;

    // vec2 = DecryptWithoutDecode(*cryptoContext, ctxtNew, keyPair.secretKey, numSlots, cryptoContext->GetRingDimension(),
    //                             true);
    // std::cerr << "\nDecrypted initial without decoding: " << vec2 << std::endl << std::endl;

    std::cout << "Number of levels remaining: "
              << depth - ctxtAfterFuncBT->GetLevel() + 1 - ctxtAfterFuncBT->GetNoiseScaleDeg() << std::endl;
    std::cout << "CKKS modulus after func bootstrapping: " << ctxtAfterFuncBT->GetElements()[0].GetModulus() << ", "
              << ctxtAfterFuncBT->GetElements()[0].GetWorkingModulus() << std::endl;

    if (QPrime != ctxtAfterFuncBT->GetElements()[0].GetModulus()) {
        OPENFHE_THROW("The ciphertext modulus after bootstrapping is not as expected.");
    }

    // =======Step 6. Modulus switch to align with the initial BFV ciphertext=======

    // Go from Q' (larger) to Q. This changes the message scaling, so we need to set scaleMod above by Q'/(q'*p)
    auto b = ctxtAfterFuncBT->GetElements()[0];
    auto a = ctxtAfterFuncBT->GetElements()[1];
    a.SetFormat(Format::COEFFICIENT);
    b.SetFormat(Format::COEFFICIENT);

    auto aPoly = a.CRTInterpolate();
    auto bPoly = b.CRTInterpolate();

    // Do modulus switching from Q' to Q
    if (Q < QPrime) {
        bPoly = bPoly.MultiplyAndRound(Q, QPrime);
        bPoly.SwitchModulus(Q, 1, 0, 0);

        aPoly = aPoly.MultiplyAndRound(Q, QPrime);
        aPoly.SwitchModulus(Q, 1, 0, 0);
    }
    else {
        bPoly.SwitchModulus(Q, 1, 0, 0);
        bPoly = bPoly.MultiplyAndRound(Q, QPrime);

        aPoly.SwitchModulus(Q, 1, 0, 0);
        aPoly = aPoly.MultiplyAndRound(Q, QPrime);
    }

    ctxtBFV[0] = ctxtBFV[0] - bPoly;
    ctxtBFV[1] = ctxtBFV[1] - aPoly;

    decrypted = DecryptBFVCoeff({bPoly, aPoly}, Q, p, keyPair.secretKey, elementParams, numSlots);
    std::cerr << "\nPlaintext after BFV decryption of digit mod Q: " << decrypted << std::endl;

    decrypted = DecryptBFVCoeff(ctxtBFV, Q, p, keyPair.secretKey, elementParams, numSlots);
    std::cerr << "\nPlaintext after BFV decryption of ciphertext mod Q: " << decrypted << std::endl;
}

void Sign() {
    std::cerr << "\n=======Sign evaluation=======\n";

    // BFV parameters
    BigInteger Q("1152921504606846976");  // 2^60
    BigInteger p("1048576");              // 2^20
    // BigInteger Bigq          = BigInteger("35184372088832"); // Mod 2^45
    // BigInteger Bigq = BigInteger("2199023255552");  // Mod 2^41
    // BigInteger pNew("2");                           // 2
    BigInteger Bigq = BigInteger("4398046511104");  // Mod 2^42
    BigInteger pNew("4");                           // 4

    uint32_t dcrtBits = Bigq.GetMSB() - 1;
    uint32_t firstMod = Bigq.GetMSB() - 1;
    uint32_t numSlots = 8;

    CCParams<CryptoContextCKKSRNS> parameters;
    SecretKeyDist secretKeyDist = SPARSE_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(FIXEDMANUAL);
    parameters.SetFirstModSize(firstMod);
    parameters.SetNumLargeDigits(3);
    parameters.SetBatchSize(numSlots);
    parameters.SetRingDim(32);

    std::vector<uint32_t> levelBudget = {1, 1};

    // Andreea: setting levelsAvailableAfterBootstrap different than 2 creates issues.
    // When 1, fresh CKKS ciphertexts (no modulus reduction is being done) evaluate to zeroes
    // When 3, the floor evaluation from BFV evaluates to zeroes, but fresh CKKS ciphertexts evaluate correctly
    uint32_t levelsAvailableAfterBootstrap = 2;
    usint depth                            = levelsAvailableAfterBootstrap + 8 + pNew.GetMSB() - 1;
    depth                                  = (pNew.GetMSB() > 2) ? depth + 1 : depth;
    parameters.SetMultiplicativeDepth(depth);

    CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);

    cryptoContext->Enable(PKE);
    cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);
    cryptoContext->Enable(ADVANCEDSHE);
    cryptoContext->Enable(FHE);

    usint ringDim = cryptoContext->GetRingDimension();
    std::cout << "CKKS scheme is using ring dimension " << ringDim << ", depth " << depth << ", firstModSize "
              << firstMod << " and scalingModSize " << dcrtBits << std::endl
              << std::endl;

    auto keyPair = cryptoContext->KeyGen();
    cryptoContext->EvalMultKeyGen(keyPair.secretKey);

    // More CKKS parameters
    // The scaling done in decoding needs to align the CKKS ciphertext after bootstrapping with the BFV ciphertext
    BigInteger QPrime = keyPair.publicKey->GetPublicElements()[0].GetParams()->GetParams()[0]->GetModulus();
    uint32_t cnt      = 1;
    auto levels       = levelsAvailableAfterBootstrap;
    while (levels > 1) {
        QPrime *= keyPair.publicKey->GetPublicElements()[0].GetParams()->GetParams()[cnt]->GetModulus();
        levels--;
        cnt++;
    }
    std::cerr << "QPrime: " << QPrime << std::endl;
    double scaleMod = QPrime.ConvertToDouble() / (Bigq.ConvertToDouble() * p.ConvertToDouble());
    // double scaleMod = 1.0; // For fresh CKKS ciphertexts, not originating from BFV
    std::cerr << "scaleMod = " << scaleMod << std::endl;
    EvalFuncBTSetup(*cryptoContext, numSlots, pNew.GetMSB() - 1, {0, 0}, levelBudget, scaleMod);
    EvalFuncBTKeyGen(keyPair.secretKey, numSlots);

    // =======Step 1. Create a CKKS ciphertext with proper metadata=======
    // We create the ciphertext modulus to be used in BFV

    // std::vector<double> y = {0.0, 1.0, 0.0, 1.0};
    // std::vector<std::complex<double>> y = { std::exp(Pi*1i/8.0) + std::exp(Pi*3i/8.0), std::exp(Pi*5i/8.0) + std::exp(Pi*15i/8.0), std::exp(Pi*9i/8.0) + std::exp(Pi*27i/8.0), std::exp(Pi*13i/8.0) + std::exp(Pi*39i/8.0)};
    // std::vector<double> y = {0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0};
    std::vector<std::complex<double>> y = {
        std::exp(Pi * 1i / 16.0) + std::exp(Pi * 3i / 16.0),   std::exp(Pi * 5i / 16.0) + std::exp(Pi * 15i / 16.0),
        std::exp(Pi * 25i / 16.0) + std::exp(Pi * 75i / 16.0), std::exp(Pi * 29i / 16.0) + std::exp(Pi * 87i / 16.0),
        std::exp(Pi * 17i / 16.0) + std::exp(Pi * 51i / 16.0), std::exp(Pi * 21i / 16.0) + std::exp(Pi * 63i / 16.0),
        std::exp(Pi * 9i / 16.0) + std::exp(Pi * 27i / 16.0),  std::exp(Pi * 13i / 16.0) + std::exp(Pi * 39i / 16.0)};
    std::cerr << "y = " << y << std::endl;
    // depth - 1 means we have two RNS limbs here; we need to the second limb
    // for internal downscaling (scalar multiplication)
    // so that the sine wave approximation of modular reduction
    // could achieve reasonable precision
    Plaintext ptxt = cryptoContext->MakeCKKSPackedPlaintext(y, 1, depth - levelsAvailableAfterBootstrap + 1);
    ptxt->SetLength(numSlots);
    Ciphertext<DCRTPoly> ctxt = cryptoContext->Encrypt(keyPair.publicKey, ptxt);
    std::cerr << "Initial CKKS ciphertext levels: " << ctxt->GetLevel() << " and depth: " << ctxt->GetNoiseScaleDeg()
              << std::endl;

    // Switch encryped digit from q to q'
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(ctxt->GetCryptoParameters());
    // auto elementParams = cryptoParams->GetElementParams(); // Andreea: this still takes the full chain of moduli
    auto elementParams = ptxt->GetElement<DCRTPoly>().GetParams();

    auto qPrimeCKKS = elementParams->GetModulus();
    std::cerr << "\nqPrimeCKKS: " << qPrimeCKKS << std::endl << std::endl;

    // =======Step 2. Encrypting and decrypting using BFV-like encryption=======
    std::cerr << "\nBigQBFV: " << Q << std::endl << std::endl;

    std::vector<int64_t> x = {524288, 1, 2, 3, 16, 33, 64, 1048575};

    std::cerr << "Plaintext before BFV encryption: " << x << std::endl;

    // // This creates BFV ciphertext with the full CKKS modulus - we don't need that
    // auto ctxtBFV = EncryptBFVCoeff(x, Q, p, keyPair.secretKey);
    // auto decrypted = DecryptBFVCoeff(encrypted, Q, p, keyPair.secretKey, numSlots);

    // This creates BFV ciphertext with the ciphertext modulus corresponding to the last two limbs
    auto ctxtBFV   = EncryptBFVCoeff(x, Q, p, keyPair.secretKey, elementParams);
    auto decrypted = DecryptBFVCoeff(ctxtBFV, Q, p, keyPair.secretKey, elementParams, numSlots);
    std::cerr << "Plaintext after BFV encryption + decryption: " << decrypted << std::endl;

    // =======Start sign loop=======

    double QBFVDouble   = Q.ConvertToDouble();
    double pBFVDouble   = p.ConvertToDouble();
    double pDigitDouble = pNew.ConvertToDouble();
    double qDigitDouble = Bigq.ConvertToDouble();
    BigInteger pOrig    = p;

    uint32_t iter = 0;
    bool step     = false;
    bool go       = QBFVDouble > qDigitDouble;

    // For arbitrary digit size, pNew > 2, the last iteration needs to evaluate step pNew not mod pNew
    while (go) {
        std::cerr << "\n=======Iteration " << iter << "=======\n";
        std::cerr << "step = " << step << std::endl;

        // =======Step 3. Changing (\log Q, \log p) from  (41,60) to (41,2), i.e., doing mod q=======

        auto encryptedDigit = ctxtBFV;

        // Apply mod q
        encryptedDigit[0].SwitchModulus(Bigq, 1, 0, 0);
        encryptedDigit[1].SwitchModulus(Bigq, 1, 0, 0);

        decrypted = DecryptBFVCoeff(encryptedDigit, Bigq, pNew, keyPair.secretKey, elementParams, numSlots);
        std::cerr << "Plaintext after BFV decryption of ciphertext mod q: " << decrypted << std::endl;

        // =======Step 4. Populate CKKS ciphertext from the BFV ciphertext=======
        std::vector<lbcrypto::DCRTPoly> elementsCKKS(2);
        if (qPrimeCKKS > Bigq) {
            elementsCKKS = ModSwitchUp(encryptedDigit, Bigq, qPrimeCKKS, elementParams);
        }
        else {
            elementsCKKS = ModSwitchDown(encryptedDigit, Bigq, qPrimeCKKS, elementParams);
        }
        auto ctxtNew = ctxt->Clone();
        ctxtNew->SetElements(elementsCKKS);
        // ctxt->SetElements(elementsCKKS);
        std::cerr << "Number of elements of populated ciphertext: " << ctxtNew->GetElements()[0].GetNumOfElements()
                  << std::endl;
        std::cerr << "Populated CKKS ciphertext levels: " << ctxtNew->GetLevel()
                  << " and depth: " << ctxtNew->GetNoiseScaleDeg() << std::endl;

        // =======Check decryption correctness=======

        Plaintext result;
        // cryptoContext->Decrypt(keyPair.secretKey, ctxt, &result);
        // result->SetLength(numSlots);
        // std::cout << "\n Decrypting CKKS initial ciphertext" << result << std::endl;
        // std::cout << "Decrypting BFV ciphertext switched to CKKS GetElement<Poly> = " << result->GetElement<Poly>() << std::endl;

        auto initial = DecryptCKKSCoeff({ctxt->GetElements()[0], ctxt->GetElements()[1]}, keyPair.secretKey, numSlots);
        std::cerr << "Decrypting CKKS initial ciphertext: " << initial << std::endl;

        auto output = DecryptCKKSCoeff(elementsCKKS, keyPair.secretKey, numSlots);
        std::cerr << "Decrypting BFV ciphertext switched to CKKS: " << output << std::endl;

        auto vec1 = DecryptWithoutDecode(*cryptoContext, ctxtNew, keyPair.secretKey, numSlots,
                                         cryptoContext->GetRingDimension(), false);
        std::cerr << "\nDecrypt without decoding switched to CKKS: " << std::setprecision(15) << vec1 << std::endl;

        std::cout << "\nInitial number of levels remaining: " << depth - ctxtNew->GetLevel() << std::endl;

        // =======Step 5. Bootstrap the digit=======

        // scaleMod = QPrime.ConvertToDouble() / (Bigq.ConvertToDouble() * p.ConvertToDouble());
        // double scaleMod = 1.0 / pNew.ConvertToDouble(); // For fresh CKKS ciphertexts, not originating from BFV
        // std::cerr << "scaleMod = " << scaleMod << std::endl;

        // Andreea: to not have to modify the decoding matrix, we multiply the BFV ciphertext by pOrig/pNew
        auto ctxtAfterFuncBT =
            EvalFuncBT(ctxtNew, pNew.GetMSB() - 1, qPrimeCKKS, pOrig.ConvertToDouble() / pBFVDouble, step);

        cryptoContext->ModReduceInPlace(ctxtAfterFuncBT);
        std::cout << "Number of levels remaining: "
                  << depth - ctxtAfterFuncBT->GetLevel() + 1 - ctxtAfterFuncBT->GetNoiseScaleDeg() << std::endl;
        std::cout << "CKKS modulus after func bootstrapping: " << ctxtAfterFuncBT->GetElements()[0].GetModulus() << ", "
                  << ctxtAfterFuncBT->GetElements()[0].GetWorkingModulus() << std::endl;

        if (QPrime != ctxtAfterFuncBT->GetElements()[0].GetModulus()) {
            OPENFHE_THROW("The ciphertext modulus after bootstrapping is not as expected.");
        }

        auto vec2 = DecryptWithoutDecode(*cryptoContext, ctxtAfterFuncBT, keyPair.secretKey, numSlots,
                                         cryptoContext->GetRingDimension(), true);

        std::cerr << "\nDecrypted func bootstrapped digit without decoding: " << vec2 << std::endl;

        elementsCKKS[0] = ctxtAfterFuncBT->GetElements()[0];
        elementsCKKS[1] = ctxtAfterFuncBT->GetElements()[1];

        output = DecryptCKKSCoeff(elementsCKKS, keyPair.secretKey, numSlots);
        std::cerr << "Decrypting func bootstrapped digit without scaling: " << output << std::endl;

        cryptoContext->Decrypt(keyPair.secretKey, ctxtAfterFuncBT, &result);
        result->SetLength(numSlots);
        std::cout << "\nFull decryption func bootstrapped digit " << result << std::endl;

        // vec2 = DecryptWithoutDecode(*cryptoContext, ctxtNew, keyPair.secretKey, numSlots, cryptoContext->GetRingDimension(),
        //                             true);
        // std::cerr << "\nDecrypted initial without decoding: " << vec2 << std::endl << std::endl;

        // =======Step 6. Modulus switch to align with the initial BFV ciphertext=======

        // Go from Q' (larger) to Q. This changes the message scaling, so we need to set scaleMod above by Q'/(q'*p)
        auto b = ctxtAfterFuncBT->GetElements()[0];
        auto a = ctxtAfterFuncBT->GetElements()[1];
        a.SetFormat(Format::COEFFICIENT);
        b.SetFormat(Format::COEFFICIENT);

        auto aPoly = a.CRTInterpolate();
        auto bPoly = b.CRTInterpolate();

        // Do modulus switching from Q' to Q
        if (Q < QPrime) {
            bPoly = bPoly.MultiplyAndRound(Q, QPrime);
            bPoly.SwitchModulus(Q, 1, 0, 0);

            aPoly = aPoly.MultiplyAndRound(Q, QPrime);
            aPoly.SwitchModulus(Q, 1, 0, 0);
        }
        else {
            bPoly.SwitchModulus(Q, 1, 0, 0);
            bPoly = bPoly.MultiplyAndRound(Q, QPrime);

            aPoly.SwitchModulus(Q, 1, 0, 0);
            aPoly = aPoly.MultiplyAndRound(Q, QPrime);
        }

        BigInteger QNew(std::to_string(static_cast<uint64_t>(QBFVDouble / pDigitDouble)));
        BigInteger PNew(std::to_string(static_cast<uint64_t>(pBFVDouble / pDigitDouble)));
        std::cerr << "QNew = " << QNew << std::endl;

        if (!step) {
            ctxtBFV[0] = ctxtBFV[0] - bPoly;
            ctxtBFV[1] = ctxtBFV[1] - aPoly;

            decrypted = DecryptBFVCoeff({bPoly, aPoly}, Q, p, keyPair.secretKey, elementParams, numSlots);
            std::cerr << "\nPlaintext after BFV decryption of digit mod Q: " << decrypted << std::endl;

            decrypted = DecryptBFVCoeff(ctxtBFV, Q, p, keyPair.secretKey, elementParams, numSlots);
            std::cerr << "\nPlaintext after BFV decryption of ciphertext mod Q: " << decrypted << std::endl;

            // Do modulus switching from Q to QNew for the BFV ciphertext
            ctxtBFV[0] = ctxtBFV[0].MultiplyAndRound(QNew, Q);
            ctxtBFV[0].SwitchModulus(QNew, 1, 0, 0);

            ctxtBFV[1] = ctxtBFV[1].MultiplyAndRound(QNew, Q);
            ctxtBFV[1].SwitchModulus(QNew, 1, 0, 0);

            decrypted = DecryptBFVCoeff(ctxtBFV, QNew, PNew, keyPair.secretKey, elementParams, numSlots);
            std::cerr << "\nPlaintext after BFV decryption of digit mod Q/pDigit: " << decrypted << std::endl;

            QBFVDouble /= pDigitDouble;
            pBFVDouble /= pDigitDouble;
            Q = QNew;
            p = PNew;
            iter++;
        }
        else {
            ctxtBFV[0] = bPoly;
            ctxtBFV[1] = aPoly;

            decrypted = DecryptBFVCoeff(ctxtBFV, Q, p, keyPair.secretKey, elementParams, numSlots);
            std::cerr << "\nPlaintext after BFV decryption of ciphertext mod Q: " << decrypted << std::endl;
        }

        // if (iter == 1) {
        //     return;
        // }
        go = QBFVDouble > qDigitDouble;

        if (pNew > 2 && !go && !step) {
            step = true;
            go   = true;
        }
    }
}

int main(int argc, char* argv[]) {
    // SimpleBootstrapExample();
    // TestModApprox();
    // Floor();
    Sign();
}
