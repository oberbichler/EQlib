import json

data = json.loads("""
{"vertices": [[1.4369117364812567, 8.1753614070371672, 2.7755575615628889e-17], [-3.5527136788005009e-15, 8.1753614070371672, 2.7755575615628889e-17], [1.8353189009063493e-15, 1.0038720130629026e-30, 0.0], [1.4369117364812611, 0.0, 2.7369110631344083e-48], [0.3225033664864505, 7.9558136499935159, -0.0053725694203199765], [0.70840373143978042, 7.8978612538643702, -0.0067907267974224583], [1.0886204972024132, 7.9919461045168685, -0.0044883694500389964], [0.048750995113513829, 0.50446428539037957, 0.15120936754855027], [0.091077893406520422, 1.028687480235406, 0.28249332403549199], [0.1267340113119656, 1.5696143959346622, 0.39308673910661185], [0.1555115243615996, 2.1240922109769884, 0.48234501040391714], [0.17724271247775714, 2.6888893123827984, 0.54974792604639322], [0.19180092089590775, 3.2607139248003492, 0.59490264509209001], [0.19910129951286432, 3.8362332686673302, 0.61754599075025518], [0.19910129972124813, 4.4120930831128504, 0.61754599139659205], [0.19180092103786253, 4.9849370707267946, 0.59490264553239025], [0.17724271264507793, 5.5514265150124817, 0.54974792656537064], [0.15551152356473666, 6.1082597750590368, 0.48234500793232615], [0.12673400899020748, 6.6521914692745279, 0.3930867319053063], [0.09107789375294148, 7.1800513861428348, 0.28249332510998026], [0.048750991662628168, 7.6887630159440796, 0.15120935684506459], [1.0886204972024178, 0.18341530252029922, 0.0044883694500390181], [0.70840373143978574, 0.27750015317280052, 0.0067907267974224757], [0.32250336648645661, 0.21954775704364921, 0.0053725694203200147], [1.3881607448186266, 7.6887630159440796, 0.15120935684506459], [1.3458338427283165, 7.1800513861428348, 0.28249332510998026], [1.3101777274910467, 6.6521914692745279, 0.3930867319053063], [1.2814002129165172, 6.1082597750590368, 0.48234500793232615], [1.2596690238361792, 5.5514265150124817, 0.54974792656537064], [1.245110815443391, 4.9849370707267946, 0.59490264553239025], [1.2378104367600107, 4.4120930831128504, 0.61754599139659205], [1.2378104369683971, 3.8362332686673302, 0.61754599075025518], [1.245110815585349, 3.2607139248003492, 0.59490264509209001], [1.259669024003502, 2.6888893123827984, 0.54974792604639322], [1.281400212119661, 2.1240922109769884, 0.48234501040391714], [1.3101777251692965, 1.569614395934662, 0.39308673910661185], [1.3458338430747396, 1.028687480235406, 0.28249332403549199], [1.3881607413677477, 0.50446428539037957, 0.15120936754855027], [1.0642450013710985, 7.5298030870931534, 0.14731943665503078], [1.0430815503259441, 7.0455468309612819, 0.27920185417995164], [1.0252534927073083, 6.5421422877623483, 0.39039371023528285], [1.0108647354200433, 6.0226659672162306, 0.48025043552230795], [0.99999914087987518, 5.4902880808390488, 0.54825180341535762], [0.99272003668348008, 4.948254010222735, 0.59400497164238242], [0.98906984734179126, 4.3998653962781651, 0.61724676676658941], [0.98906984744598525, 3.8484609555020177, 0.61784521538025783], [0.99272003675446019, 3.2973969853044096, 0.59580031898209784], [0.99999914096353737, 2.7500277465562317, 0.55124404919640624], [1.0108647350216171, 2.2096860188197951, 0.48443958281393534], [1.0252534915464353, 1.679663577446842, 0.39577976077663524], [1.0430815504991564, 1.1631920354169591, 0.28578479496552062], [1.0642449996456607, 0.66342421424130571, 0.15509928773858409], [0.7084037314397813, 7.4482628831943227, 0.14532406028729844], [0.70840373143978286, 6.9765512738161188, 0.27751345879187045], [0.70840373143978097, 6.485691377370852, 0.38901229582685282], [0.70840373143978086, 5.9787597035783993, 0.47917600209352901], [0.70840373143978241, 5.4589264639548842, 0.54748435096622983], [0.70840373143978064, 4.9294370400922372, 0.59354450017290572], [0.70840373143978319, 4.3935930729013331, 0.61709327627676402], [0.70840373143978452, 3.8547332788788524, 0.61799870587008332], [0.7084037314397823, 3.3162139554349106, 0.59626079045157454], [0.70840373143978363, 2.7813893634403994, 0.55201150164553414], [0.7084037314397843, 2.2535922824576295, 0.48551401624271434], [0.70840373143978519, 1.7361144878383434, 0.39716117518506533], [0.70840373143978419, 1.2321875925621271, 0.2874731903536018], [0.70840373143978497, 0.74496441814014036, 0.15709466410631645], [0.34687886231776688, 7.4984882931729153, 0.1465531300141206], [0.36804231336292431, 7.0190496976441583, 0.27855344086841227], [0.38587037098155641, 6.5204628150483384, 0.38986319025311428], [0.40025912826882099, 6.0058041551053334, 0.47983780886951016], [0.41112472280899248, 5.4782439293312652, 0.54795707009193062], [0.41840382700538387, 4.9410275193180642, 0.59382813164832626], [0.42205401634707801, 4.3974565659766069, 0.61718782010190421], [0.42205401624288674, 3.8508697858035732, 0.61790416204494314], [0.41840382693440742, 3.3046234762090787, 0.595977158976154], [0.41112472272533285, 2.7620718980640144, 0.55153878251983335], [0.40025912866725455, 2.2265478309306914, 0.48485220946673319], [0.38587037214243802, 1.701343050160852, 0.39631028075880387], [0.36804231318971503, 1.1896891687340825, 0.28643320827705998], [0.34687886404321217, 0.69473900816154244, 0.15586559437949432]], "trusses": [{"node_a": [1.0886204972024132, 7.9919461045168685, -0.0044883694500389964], "node_b": [1.0642450013710985, 7.5298030870931534, 0.14731943665503078]}, {"node_a": [0.70840373143978042, 7.8978612538643702, -0.0067907267974224583], "node_b": [0.7084037314397813, 7.4482628831943227, 0.14532406028729844]}, {"node_a": [0.3225033664864505, 7.9558136499935159, -0.0053725694203199765], "node_b": [0.34687886231776688, 7.4984882931729153, 0.1465531300141206]}, {"node_a": [1.3881607448186266, 7.6887630159440796, 0.15120935684506459], "node_b": [1.0642450013710985, 7.5298030870931534, 0.14731943665503078]}, {"node_a": [1.0642450013710985, 7.5298030870931534, 0.14731943665503078], "node_b": [1.0430815503259441, 7.0455468309612819, 0.27920185417995164]}, {"node_a": [1.3458338427283165, 7.1800513861428348, 0.28249332510998026], "node_b": [1.0430815503259441, 7.0455468309612819, 0.27920185417995164]}, {"node_a": [1.0430815503259441, 7.0455468309612819, 0.27920185417995164], "node_b": [1.0252534927073083, 6.5421422877623483, 0.39039371023528285]}, {"node_a": [1.3101777274910467, 6.6521914692745279, 0.3930867319053063], "node_b": [1.0252534927073083, 6.5421422877623483, 0.39039371023528285]}, {"node_a": [1.0252534927073083, 6.5421422877623483, 0.39039371023528285], "node_b": [1.0108647354200433, 6.0226659672162306, 0.48025043552230795]}, {"node_a": [1.2814002129165172, 6.1082597750590368, 0.48234500793232615], "node_b": [1.0108647354200433, 6.0226659672162306, 0.48025043552230795]}, {"node_a": [1.0108647354200433, 6.0226659672162306, 0.48025043552230795], "node_b": [0.99999914087987518, 5.4902880808390488, 0.54825180341535762]}, {"node_a": [1.2596690238361792, 5.5514265150124817, 0.54974792656537064], "node_b": [0.99999914087987518, 5.4902880808390488, 0.54825180341535762]}, {"node_a": [0.99999914087987518, 5.4902880808390488, 0.54825180341535762], "node_b": [0.99272003668348008, 4.948254010222735, 0.59400497164238242]}, {"node_a": [1.245110815443391, 4.9849370707267946, 0.59490264553239025], "node_b": [0.99272003668348008, 4.948254010222735, 0.59400497164238242]}, {"node_a": [0.99272003668348008, 4.948254010222735, 0.59400497164238242], "node_b": [0.98906984734179126, 4.3998653962781651, 0.61724676676658941]}, {"node_a": [1.2378104367600107, 4.4120930831128504, 0.61754599139659205], "node_b": [0.98906984734179126, 4.3998653962781651, 0.61724676676658941]}, {"node_a": [0.98906984734179126, 4.3998653962781651, 0.61724676676658941], "node_b": [0.98906984744598525, 3.8484609555020177, 0.61784521538025783]}, {"node_a": [1.2378104369683971, 3.8362332686673302, 0.61754599075025518], "node_b": [0.98906984744598525, 3.8484609555020177, 0.61784521538025783]}, {"node_a": [0.98906984744598525, 3.8484609555020177, 0.61784521538025783], "node_b": [0.99272003675446019, 3.2973969853044096, 0.59580031898209784]}, {"node_a": [1.245110815585349, 3.2607139248003492, 0.59490264509209001], "node_b": [0.99272003675446019, 3.2973969853044096, 0.59580031898209784]}, {"node_a": [0.99272003675446019, 3.2973969853044096, 0.59580031898209784], "node_b": [0.99999914096353737, 2.7500277465562317, 0.55124404919640624]}, {"node_a": [1.259669024003502, 2.6888893123827984, 0.54974792604639322], "node_b": [0.99999914096353737, 2.7500277465562317, 0.55124404919640624]}, {"node_a": [0.99999914096353737, 2.7500277465562317, 0.55124404919640624], "node_b": [1.0108647350216171, 2.2096860188197951, 0.48443958281393534]}, {"node_a": [1.281400212119661, 2.1240922109769884, 0.48234501040391714], "node_b": [1.0108647350216171, 2.2096860188197951, 0.48443958281393534]}, {"node_a": [1.0108647350216171, 2.2096860188197951, 0.48443958281393534], "node_b": [1.0252534915464353, 1.679663577446842, 0.39577976077663524]}, {"node_a": [1.3101777251692965, 1.569614395934662, 0.39308673910661185], "node_b": [1.0252534915464353, 1.679663577446842, 0.39577976077663524]}, {"node_a": [1.0252534915464353, 1.679663577446842, 0.39577976077663524], "node_b": [1.0430815504991564, 1.1631920354169591, 0.28578479496552062]}, {"node_a": [1.3458338430747396, 1.028687480235406, 0.28249332403549199], "node_b": [1.0430815504991564, 1.1631920354169591, 0.28578479496552062]}, {"node_a": [1.0430815504991564, 1.1631920354169591, 0.28578479496552062], "node_b": [1.0642449996456607, 0.66342421424130571, 0.15509928773858409]}, {"node_a": [1.3881607413677477, 0.50446428539037957, 0.15120936754855027], "node_b": [1.0642449996456607, 0.66342421424130571, 0.15509928773858409]}, {"node_a": [1.0642449996456607, 0.66342421424130571, 0.15509928773858409], "node_b": [1.0886204972024178, 0.18341530252029922, 0.0044883694500390181]}, {"node_a": [1.0642450013710985, 7.5298030870931534, 0.14731943665503078], "node_b": [0.7084037314397813, 7.4482628831943227, 0.14532406028729844]}, {"node_a": [0.7084037314397813, 7.4482628831943227, 0.14532406028729844], "node_b": [0.70840373143978286, 6.9765512738161188, 0.27751345879187045]}, {"node_a": [1.0430815503259441, 7.0455468309612819, 0.27920185417995164], "node_b": [0.70840373143978286, 6.9765512738161188, 0.27751345879187045]}, {"node_a": [0.70840373143978286, 6.9765512738161188, 0.27751345879187045], "node_b": [0.70840373143978097, 6.485691377370852, 0.38901229582685282]}, {"node_a": [1.0252534927073083, 6.5421422877623483, 0.39039371023528285], "node_b": [0.70840373143978097, 6.485691377370852, 0.38901229582685282]}, {"node_a": [0.70840373143978097, 6.485691377370852, 0.38901229582685282], "node_b": [0.70840373143978086, 5.9787597035783993, 0.47917600209352901]}, {"node_a": [1.0108647354200433, 6.0226659672162306, 0.48025043552230795], "node_b": [0.70840373143978086, 5.9787597035783993, 0.47917600209352901]}, {"node_a": [0.70840373143978086, 5.9787597035783993, 0.47917600209352901], "node_b": [0.70840373143978241, 5.4589264639548842, 0.54748435096622983]}, {"node_a": [0.99999914087987518, 5.4902880808390488, 0.54825180341535762], "node_b": [0.70840373143978241, 5.4589264639548842, 0.54748435096622983]}, {"node_a": [0.70840373143978241, 5.4589264639548842, 0.54748435096622983], "node_b": [0.70840373143978064, 4.9294370400922372, 0.59354450017290572]}, {"node_a": [0.99272003668348008, 4.948254010222735, 0.59400497164238242], "node_b": [0.70840373143978064, 4.9294370400922372, 0.59354450017290572]}, {"node_a": [0.70840373143978064, 4.9294370400922372, 0.59354450017290572], "node_b": [0.70840373143978319, 4.3935930729013331, 0.61709327627676402]}, {"node_a": [0.98906984734179126, 4.3998653962781651, 0.61724676676658941], "node_b": [0.70840373143978319, 4.3935930729013331, 0.61709327627676402]}, {"node_a": [0.70840373143978319, 4.3935930729013331, 0.61709327627676402], "node_b": [0.70840373143978452, 3.8547332788788524, 0.61799870587008332]}, {"node_a": [0.98906984744598525, 3.8484609555020177, 0.61784521538025783], "node_b": [0.70840373143978452, 3.8547332788788524, 0.61799870587008332]}, {"node_a": [0.70840373143978452, 3.8547332788788524, 0.61799870587008332], "node_b": [0.7084037314397823, 3.3162139554349106, 0.59626079045157454]}, {"node_a": [0.99272003675446019, 3.2973969853044096, 0.59580031898209784], "node_b": [0.7084037314397823, 3.3162139554349106, 0.59626079045157454]}, {"node_a": [0.7084037314397823, 3.3162139554349106, 0.59626079045157454], "node_b": [0.70840373143978363, 2.7813893634403994, 0.55201150164553414]}, {"node_a": [0.99999914096353737, 2.7500277465562317, 0.55124404919640624], "node_b": [0.70840373143978363, 2.7813893634403994, 0.55201150164553414]}, {"node_a": [0.70840373143978363, 2.7813893634403994, 0.55201150164553414], "node_b": [0.7084037314397843, 2.2535922824576295, 0.48551401624271434]}, {"node_a": [1.0108647350216171, 2.2096860188197951, 0.48443958281393534], "node_b": [0.7084037314397843, 2.2535922824576295, 0.48551401624271434]}, {"node_a": [0.7084037314397843, 2.2535922824576295, 0.48551401624271434], "node_b": [0.70840373143978519, 1.7361144878383434, 0.39716117518506533]}, {"node_a": [1.0252534915464353, 1.679663577446842, 0.39577976077663524], "node_b": [0.70840373143978519, 1.7361144878383434, 0.39716117518506533]}, {"node_a": [0.70840373143978519, 1.7361144878383434, 0.39716117518506533], "node_b": [0.70840373143978419, 1.2321875925621271, 0.2874731903536018]}, {"node_a": [1.0430815504991564, 1.1631920354169591, 0.28578479496552062], "node_b": [0.70840373143978419, 1.2321875925621271, 0.2874731903536018]}, {"node_a": [0.70840373143978419, 1.2321875925621271, 0.2874731903536018], "node_b": [0.70840373143978497, 0.74496441814014036, 0.15709466410631645]}, {"node_a": [1.0642449996456607, 0.66342421424130571, 0.15509928773858409], "node_b": [0.70840373143978497, 0.74496441814014036, 0.15709466410631645]}, {"node_a": [0.70840373143978497, 0.74496441814014036, 0.15709466410631645], "node_b": [0.70840373143978574, 0.27750015317280052, 0.0067907267974224757]}, {"node_a": [0.7084037314397813, 7.4482628831943227, 0.14532406028729844], "node_b": [0.34687886231776688, 7.4984882931729153, 0.1465531300141206]}, {"node_a": [0.34687886231776688, 7.4984882931729153, 0.1465531300141206], "node_b": [0.36804231336292431, 7.0190496976441583, 0.27855344086841227]}, {"node_a": [0.70840373143978286, 6.9765512738161188, 0.27751345879187045], "node_b": [0.36804231336292431, 7.0190496976441583, 0.27855344086841227]}, {"node_a": [0.36804231336292431, 7.0190496976441583, 0.27855344086841227], "node_b": [0.38587037098155641, 6.5204628150483384, 0.38986319025311428]}, {"node_a": [0.70840373143978097, 6.485691377370852, 0.38901229582685282], "node_b": [0.38587037098155641, 6.5204628150483384, 0.38986319025311428]}, {"node_a": [0.38587037098155641, 6.5204628150483384, 0.38986319025311428], "node_b": [0.40025912826882099, 6.0058041551053334, 0.47983780886951016]}, {"node_a": [0.70840373143978086, 5.9787597035783993, 0.47917600209352901], "node_b": [0.40025912826882099, 6.0058041551053334, 0.47983780886951016]}, {"node_a": [0.40025912826882099, 6.0058041551053334, 0.47983780886951016], "node_b": [0.41112472280899248, 5.4782439293312652, 0.54795707009193062]}, {"node_a": [0.70840373143978241, 5.4589264639548842, 0.54748435096622983], "node_b": [0.41112472280899248, 5.4782439293312652, 0.54795707009193062]}, {"node_a": [0.41112472280899248, 5.4782439293312652, 0.54795707009193062], "node_b": [0.41840382700538387, 4.9410275193180642, 0.59382813164832626]}, {"node_a": [0.70840373143978064, 4.9294370400922372, 0.59354450017290572], "node_b": [0.41840382700538387, 4.9410275193180642, 0.59382813164832626]}, {"node_a": [0.41840382700538387, 4.9410275193180642, 0.59382813164832626], "node_b": [0.42205401634707801, 4.3974565659766069, 0.61718782010190421]}, {"node_a": [0.70840373143978319, 4.3935930729013331, 0.61709327627676402], "node_b": [0.42205401634707801, 4.3974565659766069, 0.61718782010190421]}, {"node_a": [0.42205401634707801, 4.3974565659766069, 0.61718782010190421], "node_b": [0.42205401624288674, 3.8508697858035732, 0.61790416204494314]}, {"node_a": [0.70840373143978452, 3.8547332788788524, 0.61799870587008332], "node_b": [0.42205401624288674, 3.8508697858035732, 0.61790416204494314]}, {"node_a": [0.42205401624288674, 3.8508697858035732, 0.61790416204494314], "node_b": [0.41840382693440742, 3.3046234762090787, 0.595977158976154]}, {"node_a": [0.7084037314397823, 3.3162139554349106, 0.59626079045157454], "node_b": [0.41840382693440742, 3.3046234762090787, 0.595977158976154]}, {"node_a": [0.41840382693440742, 3.3046234762090787, 0.595977158976154], "node_b": [0.41112472272533285, 2.7620718980640144, 0.55153878251983335]}, {"node_a": [0.70840373143978363, 2.7813893634403994, 0.55201150164553414], "node_b": [0.41112472272533285, 2.7620718980640144, 0.55153878251983335]}, {"node_a": [0.41112472272533285, 2.7620718980640144, 0.55153878251983335], "node_b": [0.40025912866725455, 2.2265478309306914, 0.48485220946673319]}, {"node_a": [0.7084037314397843, 2.2535922824576295, 0.48551401624271434], "node_b": [0.40025912866725455, 2.2265478309306914, 0.48485220946673319]}, {"node_a": [0.40025912866725455, 2.2265478309306914, 0.48485220946673319], "node_b": [0.38587037214243802, 1.701343050160852, 0.39631028075880387]}, {"node_a": [0.70840373143978519, 1.7361144878383434, 0.39716117518506533], "node_b": [0.38587037214243802, 1.701343050160852, 0.39631028075880387]}, {"node_a": [0.38587037214243802, 1.701343050160852, 0.39631028075880387], "node_b": [0.36804231318971503, 1.1896891687340825, 0.28643320827705998]}, {"node_a": [0.70840373143978419, 1.2321875925621271, 0.2874731903536018], "node_b": [0.36804231318971503, 1.1896891687340825, 0.28643320827705998]}, {"node_a": [0.36804231318971503, 1.1896891687340825, 0.28643320827705998], "node_b": [0.34687886404321217, 0.69473900816154244, 0.15586559437949432]}, {"node_a": [0.70840373143978497, 0.74496441814014036, 0.15709466410631645], "node_b": [0.34687886404321217, 0.69473900816154244, 0.15586559437949432]}, {"node_a": [0.34687886404321217, 0.69473900816154244, 0.15586559437949432], "node_b": [0.32250336648645661, 0.21954775704364921, 0.0053725694203200147]}, {"node_a": [0.34687886231776688, 7.4984882931729153, 0.1465531300141206], "node_b": [0.048750991662628168, 7.6887630159440796, 0.15120935684506459]}, {"node_a": [0.36804231336292431, 7.0190496976441583, 0.27855344086841227], "node_b": [0.09107789375294148, 7.1800513861428348, 0.28249332510998026]}, {"node_a": [0.38587037098155641, 6.5204628150483384, 0.38986319025311428], "node_b": [0.12673400899020748, 6.6521914692745279, 0.3930867319053063]}, {"node_a": [0.40025912826882099, 6.0058041551053334, 0.47983780886951016], "node_b": [0.15551152356473666, 6.1082597750590368, 0.48234500793232615]}, {"node_a": [0.41112472280899248, 5.4782439293312652, 0.54795707009193062], "node_b": [0.17724271264507793, 5.5514265150124817, 0.54974792656537064]}, {"node_a": [0.41840382700538387, 4.9410275193180642, 0.59382813164832626], "node_b": [0.19180092103786253, 4.9849370707267946, 0.59490264553239025]}, {"node_a": [0.42205401634707801, 4.3974565659766069, 0.61718782010190421], "node_b": [0.19910129972124813, 4.4120930831128504, 0.61754599139659205]}, {"node_a": [0.42205401624288674, 3.8508697858035732, 0.61790416204494314], "node_b": [0.19910129951286432, 3.8362332686673302, 0.61754599075025518]}, {"node_a": [0.41840382693440742, 3.3046234762090787, 0.595977158976154], "node_b": [0.19180092089590775, 3.2607139248003492, 0.59490264509209001]}, {"node_a": [0.41112472272533285, 2.7620718980640144, 0.55153878251983335], "node_b": [0.17724271247775714, 2.6888893123827984, 0.54974792604639322]}, {"node_a": [0.40025912866725455, 2.2265478309306914, 0.48485220946673319], "node_b": [0.1555115243615996, 2.1240922109769884, 0.48234501040391714]}, {"node_a": [0.38587037214243802, 1.701343050160852, 0.39631028075880387], "node_b": [0.1267340113119656, 1.5696143959346622, 0.39308673910661185]}, {"node_a": [0.36804231318971503, 1.1896891687340825, 0.28643320827705998], "node_b": [0.091077893406520422, 1.028687480235406, 0.28249332403549199]}, {"node_a": [0.34687886404321217, 0.69473900816154244, 0.15586559437949432], "node_b": [0.048750995113513829, 0.50446428539037957, 0.15120936754855027]}], "supports": [{"directions": ["x", "y", "z"], "node": [1.4369117364812567, 8.1753614070371672, 2.7755575615628889e-17]}, {"directions": ["x", "y", "z"], "node": [-3.5527136788005009e-15, 8.1753614070371672, 2.7755575615628889e-17]}, {"directions": ["x", "y", "z"], "node": [1.8353189009063493e-15, 1.0038720130629026e-30, 0.0]}, {"directions": ["x", "y", "z"], "node": [1.4369117364812611, 0.0, 2.7369110631344083e-48]}, {"directions": ["x", "y", "z"], "node": [0.3225033664864505, 7.9558136499935159, -0.0053725694203199765]}, {"directions": ["x", "y", "z"], "node": [0.70840373143978042, 7.8978612538643702, -0.0067907267974224583]}, {"directions": ["x", "y", "z"], "node": [1.0886204972024132, 7.9919461045168685, -0.0044883694500389964]}, {"directions": ["x", "y", "z"], "node": [0.048750995113513829, 0.50446428539037957, 0.15120936754855027]}, {"directions": ["x", "y", "z"], "node": [0.091077893406520422, 1.028687480235406, 0.28249332403549199]}, {"directions": ["x", "y", "z"], "node": [0.1267340113119656, 1.5696143959346622, 0.39308673910661185]}, {"directions": ["x", "y", "z"], "node": [0.1555115243615996, 2.1240922109769884, 0.48234501040391714]}, {"directions": ["x", "y", "z"], "node": [0.17724271247775714, 2.6888893123827984, 0.54974792604639322]}, {"directions": ["x", "y", "z"], "node": [0.19180092089590775, 3.2607139248003492, 0.59490264509209001]}, {"directions": ["x", "y", "z"], "node": [0.19910129951286432, 3.8362332686673302, 0.61754599075025518]}, {"directions": ["x", "y", "z"], "node": [0.19910129972124813, 4.4120930831128504, 0.61754599139659205]}, {"directions": ["x", "y", "z"], "node": [0.19180092103786253, 4.9849370707267946, 0.59490264553239025]}, {"directions": ["x", "y", "z"], "node": [0.17724271264507793, 5.5514265150124817, 0.54974792656537064]}, {"directions": ["x", "y", "z"], "node": [0.15551152356473666, 6.1082597750590368, 0.48234500793232615]}, {"directions": ["x", "y", "z"], "node": [0.12673400899020748, 6.6521914692745279, 0.3930867319053063]}, {"directions": ["x", "y", "z"], "node": [0.09107789375294148, 7.1800513861428348, 0.28249332510998026]}, {"directions": ["x", "y", "z"], "node": [0.048750991662628168, 7.6887630159440796, 0.15120935684506459]}, {"directions": ["x", "y", "z"], "node": [1.0886204972024178, 0.18341530252029922, 0.0044883694500390181]}, {"directions": ["x", "y", "z"], "node": [0.70840373143978574, 0.27750015317280052, 0.0067907267974224757]}, {"directions": ["x", "y", "z"], "node": [0.32250336648645661, 0.21954775704364921, 0.0053725694203200147]}, {"directions": ["x", "y", "z"], "node": [1.3881607448186266, 7.6887630159440796, 0.15120935684506459]}, {"directions": ["x", "y", "z"], "node": [1.3458338427283165, 7.1800513861428348, 0.28249332510998026]}, {"directions": ["x", "y", "z"], "node": [1.3101777274910467, 6.6521914692745279, 0.3930867319053063]}, {"directions": ["x", "y", "z"], "node": [1.2814002129165172, 6.1082597750590368, 0.48234500793232615]}, {"directions": ["x", "y", "z"], "node": [1.2596690238361792, 5.5514265150124817, 0.54974792656537064]}, {"directions": ["x", "y", "z"], "node": [1.245110815443391, 4.9849370707267946, 0.59490264553239025]}, {"directions": ["x", "y", "z"], "node": [1.2378104367600107, 4.4120930831128504, 0.61754599139659205]}, {"directions": ["x", "y", "z"], "node": [1.2378104369683971, 3.8362332686673302, 0.61754599075025518]}, {"directions": ["x", "y", "z"], "node": [1.245110815585349, 3.2607139248003492, 0.59490264509209001]}, {"directions": ["x", "y", "z"], "node": [1.259669024003502, 2.6888893123827984, 0.54974792604639322]}, {"directions": ["x", "y", "z"], "node": [1.281400212119661, 2.1240922109769884, 0.48234501040391714]}, {"directions": ["x", "y", "z"], "node": [1.3101777251692965, 1.569614395934662, 0.39308673910661185]}, {"directions": ["x", "y", "z"], "node": [1.3458338430747396, 1.028687480235406, 0.28249332403549199]}, {"directions": ["x", "y", "z"], "node": [1.3881607413677477, 0.50446428539037957, 0.15120936754855027]}], "loads": [{"node": [1.4369117364812567, 8.1753614070371672, 2.7755575615628889e-17], "load": [0, 0, 1]}, {"node": [-3.5527136788005009e-15, 8.1753614070371672, 2.7755575615628889e-17], "load": [0, 0, 1]}, {"node": [1.8353189009063493e-15, 1.0038720130629026e-30, 0.0], "load": [0, 0, 1]}, {"node": [1.4369117364812611, 0.0, 2.7369110631344083e-48], "load": [0, 0, 1]}, {"node": [0.3225033664864505, 7.9558136499935159, -0.0053725694203199765], "load": [0, 0, 1]}, {"node": [0.70840373143978042, 7.8978612538643702, -0.0067907267974224583], "load": [0, 0, 1]}, {"node": [1.0886204972024132, 7.9919461045168685, -0.0044883694500389964], "load": [0, 0, 1]}, {"node": [0.048750995113513829, 0.50446428539037957, 0.15120936754855027], "load": [0, 0, 1]}, {"node": [0.091077893406520422, 1.028687480235406, 0.28249332403549199], "load": [0, 0, 1]}, {"node": [0.1267340113119656, 1.5696143959346622, 0.39308673910661185], "load": [0, 0, 1]}, {"node": [0.1555115243615996, 2.1240922109769884, 0.48234501040391714], "load": [0, 0, 1]}, {"node": [0.17724271247775714, 2.6888893123827984, 0.54974792604639322], "load": [0, 0, 1]}, {"node": [0.19180092089590775, 3.2607139248003492, 0.59490264509209001], "load": [0, 0, 1]}, {"node": [0.19910129951286432, 3.8362332686673302, 0.61754599075025518], "load": [0, 0, 1]}, {"node": [0.19910129972124813, 4.4120930831128504, 0.61754599139659205], "load": [0, 0, 1]}, {"node": [0.19180092103786253, 4.9849370707267946, 0.59490264553239025], "load": [0, 0, 1]}, {"node": [0.17724271264507793, 5.5514265150124817, 0.54974792656537064], "load": [0, 0, 1]}, {"node": [0.15551152356473666, 6.1082597750590368, 0.48234500793232615], "load": [0, 0, 1]}, {"node": [0.12673400899020748, 6.6521914692745279, 0.3930867319053063], "load": [0, 0, 1]}, {"node": [0.09107789375294148, 7.1800513861428348, 0.28249332510998026], "load": [0, 0, 1]}, {"node": [0.048750991662628168, 7.6887630159440796, 0.15120935684506459], "load": [0, 0, 1]}, {"node": [1.0886204972024178, 0.18341530252029922, 0.0044883694500390181], "load": [0, 0, 1]}, {"node": [0.70840373143978574, 0.27750015317280052, 0.0067907267974224757], "load": [0, 0, 1]}, {"node": [0.32250336648645661, 0.21954775704364921, 0.0053725694203200147], "load": [0, 0, 1]}, {"node": [1.3881607448186266, 7.6887630159440796, 0.15120935684506459], "load": [0, 0, 1]}, {"node": [1.3458338427283165, 7.1800513861428348, 0.28249332510998026], "load": [0, 0, 1]}, {"node": [1.3101777274910467, 6.6521914692745279, 0.3930867319053063], "load": [0, 0, 1]}, {"node": [1.2814002129165172, 6.1082597750590368, 0.48234500793232615], "load": [0, 0, 1]}, {"node": [1.2596690238361792, 5.5514265150124817, 0.54974792656537064], "load": [0, 0, 1]}, {"node": [1.245110815443391, 4.9849370707267946, 0.59490264553239025], "load": [0, 0, 1]}, {"node": [1.2378104367600107, 4.4120930831128504, 0.61754599139659205], "load": [0, 0, 1]}, {"node": [1.2378104369683971, 3.8362332686673302, 0.61754599075025518], "load": [0, 0, 1]}, {"node": [1.245110815585349, 3.2607139248003492, 0.59490264509209001], "load": [0, 0, 1]}, {"node": [1.259669024003502, 2.6888893123827984, 0.54974792604639322], "load": [0, 0, 1]}, {"node": [1.281400212119661, 2.1240922109769884, 0.48234501040391714], "load": [0, 0, 1]}, {"node": [1.3101777251692965, 1.569614395934662, 0.39308673910661185], "load": [0, 0, 1]}, {"node": [1.3458338430747396, 1.028687480235406, 0.28249332403549199], "load": [0, 0, 1]}, {"node": [1.3881607413677477, 0.50446428539037957, 0.15120936754855027], "load": [0, 0, 1]}, {"node": [1.0642450013710985, 7.5298030870931534, 0.14731943665503078], "load": [0, 0, 1]}, {"node": [1.0430815503259441, 7.0455468309612819, 0.27920185417995164], "load": [0, 0, 1]}, {"node": [1.0252534927073083, 6.5421422877623483, 0.39039371023528285], "load": [0, 0, 1]}, {"node": [1.0108647354200433, 6.0226659672162306, 0.48025043552230795], "load": [0, 0, 1]}, {"node": [0.99999914087987518, 5.4902880808390488, 0.54825180341535762], "load": [0, 0, 1]}, {"node": [0.99272003668348008, 4.948254010222735, 0.59400497164238242], "load": [0, 0, 1]}, {"node": [0.98906984734179126, 4.3998653962781651, 0.61724676676658941], "load": [0, 0, 1]}, {"node": [0.98906984744598525, 3.8484609555020177, 0.61784521538025783], "load": [0, 0, 1]}, {"node": [0.99272003675446019, 3.2973969853044096, 0.59580031898209784], "load": [0, 0, 1]}, {"node": [0.99999914096353737, 2.7500277465562317, 0.55124404919640624], "load": [0, 0, 1]}, {"node": [1.0108647350216171, 2.2096860188197951, 0.48443958281393534], "load": [0, 0, 1]}, {"node": [1.0252534915464353, 1.679663577446842, 0.39577976077663524], "load": [0, 0, 1]}, {"node": [1.0430815504991564, 1.1631920354169591, 0.28578479496552062], "load": [0, 0, 1]}, {"node": [1.0642449996456607, 0.66342421424130571, 0.15509928773858409], "load": [0, 0, 1]}, {"node": [0.7084037314397813, 7.4482628831943227, 0.14532406028729844], "load": [0, 0, 1]}, {"node": [0.70840373143978286, 6.9765512738161188, 0.27751345879187045], "load": [0, 0, 1]}, {"node": [0.70840373143978097, 6.485691377370852, 0.38901229582685282], "load": [0, 0, 1]}, {"node": [0.70840373143978086, 5.9787597035783993, 0.47917600209352901], "load": [0, 0, 1]}, {"node": [0.70840373143978241, 5.4589264639548842, 0.54748435096622983], "load": [0, 0, 1]}, {"node": [0.70840373143978064, 4.9294370400922372, 0.59354450017290572], "load": [0, 0, 1]}, {"node": [0.70840373143978319, 4.3935930729013331, 0.61709327627676402], "load": [0, 0, 1]}, {"node": [0.70840373143978452, 3.8547332788788524, 0.61799870587008332], "load": [0, 0, 1]}, {"node": [0.7084037314397823, 3.3162139554349106, 0.59626079045157454], "load": [0, 0, 1]}, {"node": [0.70840373143978363, 2.7813893634403994, 0.55201150164553414], "load": [0, 0, 1]}, {"node": [0.7084037314397843, 2.2535922824576295, 0.48551401624271434], "load": [0, 0, 1]}, {"node": [0.70840373143978519, 1.7361144878383434, 0.39716117518506533], "load": [0, 0, 1]}, {"node": [0.70840373143978419, 1.2321875925621271, 0.2874731903536018], "load": [0, 0, 1]}, {"node": [0.70840373143978497, 0.74496441814014036, 0.15709466410631645], "load": [0, 0, 1]}, {"node": [0.34687886231776688, 7.4984882931729153, 0.1465531300141206], "load": [0, 0, 1]}, {"node": [0.36804231336292431, 7.0190496976441583, 0.27855344086841227], "load": [0, 0, 1]}, {"node": [0.38587037098155641, 6.5204628150483384, 0.38986319025311428], "load": [0, 0, 1]}, {"node": [0.40025912826882099, 6.0058041551053334, 0.47983780886951016], "load": [0, 0, 1]}, {"node": [0.41112472280899248, 5.4782439293312652, 0.54795707009193062], "load": [0, 0, 1]}, {"node": [0.41840382700538387, 4.9410275193180642, 0.59382813164832626], "load": [0, 0, 1]}, {"node": [0.42205401634707801, 4.3974565659766069, 0.61718782010190421], "load": [0, 0, 1]}, {"node": [0.42205401624288674, 3.8508697858035732, 0.61790416204494314], "load": [0, 0, 1]}, {"node": [0.41840382693440742, 3.3046234762090787, 0.595977158976154], "load": [0, 0, 1]}, {"node": [0.41112472272533285, 2.7620718980640144, 0.55153878251983335], "load": [0, 0, 1]}, {"node": [0.40025912866725455, 2.2265478309306914, 0.48485220946673319], "load": [0, 0, 1]}, {"node": [0.38587037214243802, 1.701343050160852, 0.39631028075880387], "load": [0, 0, 1]}, {"node": [0.36804231318971503, 1.1896891687340825, 0.28643320827705998], "load": [0, 0, 1]}, {"node": [0.34687886404321217, 0.69473900816154244, 0.15586559437949432], "load": [0, 0, 1]}]}
""")

import eqlib as eq
import numpy as np
import hyperjet as hj

class Truss(eq.Objective):
  def __init__(self, node_a, node_b, youngs_modulus, area, prestress=0, ref_length=None):
    super().__init__()

    self.node_a = node_a
    self.node_b = node_b

    self.ref_length = ref_length or np.linalg.norm(node_b.location - node_a.location)
    self.youngs_modulus = youngs_modulus
    self.area = area
    self.prestress = prestress

    self.parameters = [
      node_a.x, node_a.y, node_a.z,
      node_b.x, node_b.y, node_b.z,
    ]

  def compute(self, df, hm, request):
    ax, ay, az, bx, by, bz = hj.variables(self.parameter_values)

    a = np.array([ax, ay, az])
    b = np.array([bx, by, bz])

    act_length = np.linalg.norm(b - a)
    ref_length = self.ref_length

    eps = 0.5 * (act_length**2 / ref_length**2 - 1)

    E = self.youngs_modulus
    A = self.area

    pi = (0.5 * E * eps + self.prestress) * eps * A * ref_length

    df[:] = hj.d(pi)
    hm[:] = hj.dd(pi)
    return hj.f(pi)

class NodeLoad(eq.Objective):
  def __init__(self, node, load):
    super().__init__()

    self.node = node
    self.load = np.asarray(load)

    self.parameters = [node.x, node.y, node.z]

  def compute(self, df, hm, request):
    df[:] = self.load
    hm[:] = 0
    return 0

from scipy.spatial import cKDTree


# --- nodes

vertices = np.array(data['vertices'])

index = cKDTree(vertices)

nodes = [eq.Node(*vertex) for vertex in vertices]

def get_node(location):
  _, idx = index.query(location)
  return nodes[idx]


# --- supports

for support_data in data['supports']:
  node = get_node(support_data['node'])
  directions = support_data['directions']

  if 'x' in directions: node.x.is_active = False
  if 'y' in directions: node.y.is_active = False
  if 'z' in directions: node.z.is_active = False


# --- elements

elements = []


load_factor = 1.0

for load in data['loads']:
  node = get_node(load['node'])
  load = np.multiply(load['load'], load_factor)
  
  element = NodeLoad(node, load)

  elements.append(element)


youngs_modulus = 1
area = 1
prestress = 0

for truss_data in data['trusses']:
  node_a = get_node(truss_data['node_a'])
  node_b = get_node(truss_data['node_b'])

  l = np.linalg.norm(node_b.location - node_a.location) * 0.9

  element = Truss(node_a, node_b, youngs_modulus, area, prestress, l)

  elements.append(element)


# --- solve

problem = eq.Problem(elements)

for i in range(100):
  problem.compute()

  rnorm = np.linalg.norm(problem.df)

  print(rnorm)

  if rnorm < 1e-8:
    break

  problem.newton_step()

# # --- output

# import plotly.graph_objs as go

# lines = []

# for element in elements:
#   if not isinstance(element, Truss):
#     continue

#   x, y, z = np.transpose([element.node_a.location, element.node_b.location])

#   line = go.Scatter3d(x=x, y=y, z=z, mode='lines', line=dict(color='black'))

#   lines.append(line)

# fig = go.Figure(data=lines)
# fig.show()

# # for node in nodes:
# #   x, y, z = node.location
# #   print(x, y, z)