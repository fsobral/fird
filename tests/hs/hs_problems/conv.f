      SUBROUTINE CONV(MODE)
C
      COMMON/L8/NTP
C
      IF (NTP.GT.200) GOTO 200
      GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,
     *      20,21,22,23,24,25,26,27,28,29,
     *      30,31,32,33,34,35,36,37,38,39,
     *      40,41,42,43,44,45,46,47,48,49,
     *      50,51,52,53,54,55,56,57,58,59,
     *      60,61,62,63,64,65,66,67,68,69,
     *      70,71,72,73,74,75,76,77,78,79,
     *      80,81,82,83,84,85,86,87,88,89,
     *      90,91,92,93,94,95,96,97,98,99,
     *      100,101,102,103,104,105,106,107,108,109,
     *      110,111,112,113,114,115,116,117,118,119), NTP
C
    1 CALL TP1(MODE)
      RETURN
    2 CALL TP2(MODE)
      RETURN
    3 CALL TP3(MODE)
      RETURN
    4 CALL TP4(MODE)
      RETURN
    5 CALL TP5(MODE)
      RETURN
    6 CALL TP6(MODE)
      RETURN
    7 CALL TP7(MODE)
      RETURN
    8 CALL TP8(MODE)
      RETURN
    9 CALL TP9(MODE)
      RETURN
   10 CALL TP10(MODE)
      RETURN
   11 CALL TP11(MODE)
      RETURN
   12 CALL TP12(MODE)
      RETURN
   13 CALL TP13(MODE)
      RETURN
   14 CALL TP14(MODE)
      RETURN
   15 CALL TP15(MODE)
      RETURN
   16 CALL TP16(MODE)
      RETURN
   17 CALL TP17(MODE)
      RETURN
   18 CALL TP18(MODE)
      RETURN
   19 CALL TP19(MODE)
      RETURN
   20 CALL TP20(MODE)
      RETURN
   21 CALL TP21(MODE)
      RETURN
   22 CALL TP22(MODE)
      RETURN
   23 CALL TP23(MODE)
      RETURN
   24 CALL TP24(MODE)
      RETURN
   25 CALL TP25(MODE)
      RETURN
   26 CALL TP26(MODE)
      RETURN
   27 CALL TP27(MODE)
      RETURN
   28 CALL TP28(MODE)
      RETURN
   29 CALL TP29(MODE)
      RETURN
   30 CALL TP30(MODE)
      RETURN
   31 CALL TP31(MODE)
      RETURN
   32 CALL TP32(MODE)
      RETURN
   33 CALL TP33(MODE)
      RETURN
   34 CALL TP34(MODE)
      RETURN
   35 CALL TP35(MODE)
      RETURN
   36 CALL TP36(MODE)
      RETURN
   37 CALL TP37(MODE)
      RETURN
   38 CALL TP38(MODE)
      RETURN
   39 CALL TP39(MODE)
      RETURN
   40 CALL TP40(MODE)
      RETURN
   41 CALL TP41(MODE)
      RETURN
   42 CALL TP42(MODE)
      RETURN
   43 CALL TP43(MODE)
      RETURN
   44 CALL TP44(MODE)
      RETURN
   45 CALL TP45(MODE)
      RETURN
   46 CALL TP46(MODE)
      RETURN
   47 CALL TP47(MODE)
      RETURN
   48 CALL TP48(MODE)
      RETURN
   49 CALL TP49(MODE)
      RETURN
   50 CALL TP50(MODE)
      RETURN
   51 CALL TP51(MODE)
      RETURN
   52 CALL TP52(MODE)
      RETURN
   53 CALL TP53(MODE)
      RETURN
   54 CALL TP54(MODE)
      RETURN
   55 CALL TP55(MODE)
      RETURN
   56 CALL TP56(MODE)
      RETURN
   57 CALL TP57(MODE)
      RETURN
   58 CALL TP58(MODE)
      RETURN
   59 CALL TP59(MODE)
      RETURN
   60 CALL TP60(MODE)
      RETURN
   61 CALL TP61(MODE)
      RETURN
   62 CALL TP62(MODE)
      RETURN
   63 CALL TP63(MODE)
      RETURN
   64 CALL TP64(MODE)
      RETURN
   65 CALL TP65(MODE)
      RETURN
   66 CALL TP66(MODE)
      RETURN
   67 CALL TP67(MODE)
      RETURN
   68 CALL TP68(MODE)
      RETURN
   69 CALL TP69(MODE)
      RETURN
   70 CALL TP70(MODE)
      RETURN
   71 CALL TP71(MODE)
      RETURN
   72 CALL TP72(MODE)
      RETURN
   73 CALL TP73(MODE)
      RETURN
   74 CALL TP74(MODE)
      RETURN
   75 CALL TP75(MODE)
      RETURN
   76 CALL TP76(MODE)
      RETURN
   77 CALL TP77(MODE)
      RETURN
   78 CALL TP78(MODE)
      RETURN
   79 CALL TP79(MODE)
      RETURN
   80 CALL TP80(MODE)
      RETURN
   81 CALL TP81(MODE)
      RETURN
   82 RETURN
   83 CALL TP83(MODE)
      RETURN
   84 CALL TP84(MODE)
      RETURN
   85 CALL TP85(MODE)
      RETURN
   86 CALL TP86(MODE)
      RETURN
   87 CALL TP87(MODE)
      RETURN
   88 CALL TP88(MODE)
      RETURN
   89 CALL TP89(MODE)
      RETURN
   90 CALL TP90(MODE)
      RETURN
   91 CALL TP91(MODE)
      RETURN
   92 CALL TP92(MODE)
      RETURN
   93 CALL TP93(MODE)
      RETURN
   94 RETURN
   95 CALL TP95(MODE)
      RETURN
   96 CALL TP96(MODE)
      RETURN
   97 CALL TP97(MODE)
      RETURN
   98 CALL TP98(MODE)
      RETURN
   99 CALL TP99(MODE)
      RETURN
  100 CALL TP100(MODE)
      RETURN
  101 CALL TP101(MODE)
      RETURN
  102 CALL TP102(MODE)
      RETURN
  103 CALL TP103(MODE)
      RETURN
  104 CALL TP104(MODE)
      RETURN
  105 CALL TP105(MODE)
      RETURN
  106 CALL TP106(MODE)
      RETURN
  107 CALL TP107(MODE)
      RETURN
  108 CALL TP108(MODE)
      RETURN
  109 CALL TP109(MODE)
      RETURN
  110 CALL TP110(MODE)
      RETURN
  111 CALL TP111(MODE)
      RETURN
  112 CALL TP112(MODE)
      RETURN
  113 CALL TP113(MODE)
      RETURN
  114 CALL TP114(MODE)
      RETURN
  115 RETURN
  116 CALL TP116(MODE)
      RETURN
  117 CALL TP117(MODE)
      RETURN
  118 CALL TP118(MODE)
      RETURN
  119 CALL TP119(MODE)
      RETURN
C
  200 CONTINUE
      NO=NTP-200
C
      GOTO (201,202,203,204,205,206,207,208,209,
     *      210,211,212,213,214,215,216,217,218,219,
     *      220,221,222,223,224,225,226,227,228,229,
     *      230,231,232,233,234,235,236,237,238,239,
     *      240,241,242,243,244,245,246,247,248,249,
     *      250,251,252,253,254,255,256,257,258,259,
     *      260,261,262,263,264,265,266,267,268,269,
     *      270,271,272,273,274,275,276,277,278,279,
     *      280,281,282,283,284,285,286,287,288,289,
     *      290,291,292,293,294,295,296,297,298,299,
     *      300,301,302,303,304,305,306,307,308,309,
     *      310,311,312,313,314,315,316,317,318,319,
     *      320,321,322,323,324,325,326,327,328,329,
     *      330,331,332,333,334,335,336,337,338,339,
     *      340,341,342,343,344,345,346,347,348,349,
     *      350,351,352,353,354,355,356,357,358,359,
     *      360,361,362,363,364,365,366,367,368,369,
     *      370,371,372,373,374,375,376,377,378,379,
     *      380,381,382,383,384,385,386,387,388,389,
     *      390,391,392,393,394,395), NO
C
  201 CALL TP201(MODE)
      RETURN
  202 CALL TP202(MODE)
      RETURN
  203 CALL TP203(MODE)
      RETURN
  204 CALL TP204(MODE)
      RETURN
  205 CALL TP205(MODE)
      RETURN
  206 CALL TP206(MODE)
      RETURN
  207 CALL TP207(MODE)
      RETURN
  208 CALL TP208(MODE)
      RETURN
  209 CALL TP209(MODE)
      RETURN
  210 CALL TP210(MODE)
      RETURN
  211 CALL TP211(MODE)
      RETURN
  212 CALL TP212(MODE)
      RETURN
  213 CALL TP213(MODE)
      RETURN
  214 CALL TP214(MODE)
      RETURN
  215 CALL TP215(MODE)
      RETURN
  216 CALL TP216(MODE)
      RETURN
  217 CALL TP217(MODE)
      RETURN
  218 CALL TP218(MODE)
      RETURN
  219 CALL TP219(MODE)
      RETURN
  220 CALL TP220(MODE)
      RETURN
  221 CALL TP221(MODE)
      RETURN
  222 CALL TP222(MODE)
      RETURN
  223 CALL TP223(MODE)
      RETURN
  224 CALL TP224(MODE)
      RETURN
  225 CALL TP225(MODE)
      RETURN
  226 CALL TP226(MODE)
      RETURN
  227 CALL TP227(MODE)
      RETURN
  228 CALL TP228(MODE)
      RETURN
  229 CALL TP229(MODE)
      RETURN
  230 CALL TP230(MODE)
      RETURN
  231 CALL TP231(MODE)
      RETURN
  232 CALL TP232(MODE)
      RETURN
  233 CALL TP233(MODE)
      RETURN
  234 CALL TP234(MODE)
      RETURN
  235 CALL TP235(MODE)
      RETURN
  236 CALL TP236(MODE)
      RETURN
  237 CALL TP237(MODE)
      RETURN
  238 CALL TP238(MODE)
      RETURN
  239 CALL TP239(MODE)
      RETURN
  240 CALL TP240(MODE)
      RETURN
  241 CALL TP241(MODE)
      RETURN
  242 CALL TP242(MODE)
      RETURN
  243 CALL TP243(MODE)
      RETURN
  244 CALL TP244(MODE)
      RETURN
  245 CALL TP245(MODE)
      RETURN
  246 CALL TP246(MODE)
      RETURN
  247 CALL TP247(MODE)
      RETURN
  248 CALL TP248(MODE)
      RETURN
  249 CALL TP249(MODE)
      RETURN
  250 CALL TP250(MODE)
      RETURN
  251 CALL TP251(MODE)
      RETURN
  252 CALL TP252(MODE)
      RETURN
  253 CALL TP253(MODE)
      RETURN
  254 CALL TP254(MODE)
      RETURN
  255 CALL TP255(MODE)
      RETURN
  256 CALL TP256(MODE)
      RETURN
  257 CALL TP257(MODE)
      RETURN
  258 CALL TP258(MODE)
      RETURN
  259 CALL TP259(MODE)
      RETURN
  260 CALL TP260(MODE)
      RETURN
  261 CALL TP261(MODE)
      RETURN
  262 CALL TP262(MODE)
      RETURN
  263 CALL TP263(MODE)
      RETURN
  264 CALL TP264(MODE)
      RETURN
  265 CALL TP265(MODE)
      RETURN
  266 CALL TP266(MODE)
      RETURN
  267 CALL TP267(MODE)
      RETURN
  268 CALL TP268(MODE)
      RETURN
  269 CALL TP269(MODE)
      RETURN
  270 CALL TP270(MODE)
      RETURN
  271 CALL TP271(MODE)
      RETURN
  272 CALL TP272(MODE)
      RETURN
  273 CALL TP273(MODE)
      RETURN
  274 CALL TP274(MODE)
      RETURN
  275 CALL TP275(MODE)
      RETURN
  276 CALL TP276(MODE)
      RETURN
  277 CALL TP277(MODE)
      RETURN
  278 CALL TP278(MODE)
      RETURN
  279 CALL TP279(MODE)
      RETURN
  280 CALL TP280(MODE)
      RETURN
  281 CALL TP281(MODE)
      RETURN
  282 CALL TP282(MODE)
      RETURN
  283 CALL TP283(MODE)
      RETURN
  284 CALL TP284(MODE)
      RETURN
  285 CALL TP285(MODE)
      RETURN
  286 CALL TP286(MODE)
      RETURN
  287 CALL TP287(MODE)
      RETURN
  288 CALL TP288(MODE)
      RETURN
  289 CALL TP289(MODE)
      RETURN
  290 CALL TP290(MODE)
      RETURN
  291 CALL TP291(MODE)
      RETURN
  292 CALL TP292(MODE)
      RETURN
  293 CALL TP293(MODE)
      RETURN
  294 CALL TP294(MODE)
      RETURN
  295 CALL TP295(MODE)
      RETURN
  296 CALL TP296(MODE)
      RETURN
  297 CALL TP297(MODE)
      RETURN
  298 CALL TP298(MODE)
      RETURN
  299 CALL TP299(MODE)
      RETURN
  300 CALL TP300(MODE)
      RETURN
  301 CALL TP301(MODE)
      RETURN
  302 CALL TP302(MODE)
      RETURN
  303 CALL TP303(MODE)
      RETURN
  304 CALL TP304(MODE)
      RETURN
  305 CALL TP305(MODE)
      RETURN
  306 CALL TP306(MODE)
      RETURN
  307 CALL TP307(MODE)
      RETURN
  308 CALL TP308(MODE)
      RETURN
  309 CALL TP309(MODE)
      RETURN
  310 CALL TP310(MODE)
      RETURN
  311 CALL TP311(MODE)
      RETURN
  312 CALL TP312(MODE)
      RETURN
  313 CALL TP313(MODE)
      RETURN
  314 CALL TP314(MODE)
      RETURN
  315 CALL TP315(MODE)
      RETURN
  316 CALL TP316(MODE)
      RETURN
  317 CALL TP317(MODE)
      RETURN
  318 CALL TP318(MODE)
      RETURN
  319 CALL TP319(MODE)
      RETURN
  320 CALL TP320(MODE)
      RETURN
  321 CALL TP321(MODE)
      RETURN
  322 CALL TP322(MODE)
      RETURN
  323 CALL TP323(MODE)
      RETURN
  324 CALL TP324(MODE)
      RETURN
  325 CALL TP325(MODE)
      RETURN
  326 CALL TP326(MODE)
      RETURN
  327 CALL TP327(MODE)
      RETURN
  328 CALL TP328(MODE)
      RETURN
  329 CALL TP329(MODE)
      RETURN
  330 CALL TP330(MODE)
      RETURN
  331 CALL TP331(MODE)
      RETURN
  332 CALL TP332(MODE)
      RETURN
  333 CALL TP333(MODE)
      RETURN
  334 CALL TP334(MODE)
      RETURN
  335 CALL TP335(MODE)
      RETURN
  336 CALL TP336(MODE)
      RETURN
  337 CALL TP337(MODE)
      RETURN
  338 CALL TP338(MODE)
      RETURN
  339 CALL TP339(MODE)
      RETURN
  340 CALL TP340(MODE)
      RETURN
  341 CALL TP341(MODE)
      RETURN
  342 CALL TP342(MODE)
      RETURN
  343 CALL TP343(MODE)
      RETURN
  344 CALL TP344(MODE)
      RETURN
  345 CALL TP345(MODE)
      RETURN
  346 CALL TP346(MODE)
      RETURN
  347 CALL TP347(MODE)
      RETURN
  348 CALL TP348(MODE)
      RETURN
  349 CALL TP349(MODE)
      RETURN
  350 CALL TP350(MODE)
      RETURN
  351 CALL TP351(MODE)
      RETURN
  352 CALL TP352(MODE)
      RETURN
  353 CALL TP353(MODE)
      RETURN
  354 CALL TP354(MODE)
      RETURN
  355 CALL TP355(MODE)
      RETURN
  356 CALL TP356(MODE)
      RETURN
  357 CALL TP357(MODE)
      RETURN
  358 CALL TP358(MODE)
      RETURN
  359 CALL TP359(MODE)
      RETURN
  360 CALL TP360(MODE)
      RETURN
  361 CALL TP361(MODE)
      RETURN
  362 CALL TP362(MODE)
      RETURN
  363 CALL TP363(MODE)
      RETURN
  364 CALL TP364(MODE)
      RETURN
  365 CALL TP365(MODE)
      RETURN
  366 CALL TP366(MODE)
      RETURN
  367 CALL TP367(MODE)
      RETURN
  368 CALL TP368(MODE)
      RETURN
  369 CALL TP369(MODE)
      RETURN
  370 CALL TP370(MODE)
      RETURN
  371 CALL TP371(MODE)
      RETURN
  372 CALL TP372(MODE)
      RETURN
  373 CALL TP373(MODE)
      RETURN
  374 CALL TP374(MODE)
      RETURN
  375 CALL TP375(MODE)
      RETURN
  376 CALL TP376(MODE)
      RETURN
  377 CALL TP377(MODE)
      RETURN
  378 CALL TP378(MODE)
      RETURN
  379 CALL TP379(MODE)
      RETURN
  380 CALL TP380(MODE)
      RETURN
  381 CALL TP381(MODE)
      RETURN
  382 CALL TP382(MODE)
      RETURN
  383 CALL TP383(MODE)
      RETURN
  384 CALL TP384(MODE)
      RETURN
  385 CALL TP385(MODE)
      RETURN
  386 CALL TP386(MODE)
      RETURN
  387 CALL TP387(MODE)
      RETURN
  388 CALL TP388(MODE)
      RETURN
  389 CALL TP389(MODE)
      RETURN
  390 CALL TP390(MODE)
      RETURN
  391 CALL TP391(MODE)
      RETURN
  392 CALL TP392(MODE)
      RETURN
  393 CALL TP393(MODE)
      RETURN
  394 CALL TP394(MODE)
      RETURN
  395 CALL TP395(MODE)
      RETURN
      END