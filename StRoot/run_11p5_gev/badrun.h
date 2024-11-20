const Int_t n_badrun= 374+98+94;
const Int_t badrun[n_badrun]={
  // Au+Au 19.6 GeV
20057007, 20057025, 20057026, 20057050, 20058001, 20058002, 20058003, 20058004, 20058005, 20060012, 20060022, 20060025, 20060060, 20060061, 20060062, 20062010,
20062011, 20062012, 20062036, 20063011, 20063034, 20063035, 20063036, 20063039, 20064008, 20064009, 20064011, 20064012, 20064040, 20065018, 20067014, 20067023,
20067024, 20067029, 20067030, 20067045, 20067046, 20069030, 20069032, 20069054, 20070042, 20070043, 20070044, 20070047, 20071001, 20071004, 20071005, 20071006,
20071027, 20071037, 20072034, 20072035, 20072036, 20072039, 20072041, 20072045, 20072047, 20073071, 20073072, 20073076, 20074001, 20074003, 20074004, 20074005,
20074007, 20074008, 20074009, 20074012, 20074014, 20074017, 20074018, 20074020, 20074021, 20074026, 20074027, 20074029, 20074032, 20074033, 20074034, 20074044,
20074045, 20075001, 20075002, 20075006, 20075007, 20075009, 20075011, 20075013, 20081002, 20081014, 20082060, 20082065, 20083024, 20086012, 20087007, 20089008,
20090024, 20091011, 20092054, 20062007, 20062009, 20065017, 20065056, 20065060, 20066001, 20066008, 20066015, 20066019, 20066023, 20066026, 20066066, 20066067,
20066068, 20066073, 20066078, 20067001, 20067004, 20067009, 20067012, 20067015, 20067016, 20067019, 20067028, 20067038, 20067041, 20067047, 20068001, 20068004,
20068008, 20068012, 20068019, 20068026, 20068034, 20068051, 20068055, 20068058, 20068060, 20068064, 20069001, 20069004, 20069007, 20069010, 20069020, 20069023,
20069026, 20069031, 20069033, 20069042, 20069050, 20069053, 20069057, 20069060, 20070002, 20070005, 20070010, 20070013, 20070016, 20070019, 20070041, 20070045,
20071003, 20071007, 20071010, 20071013, 20071016, 20071019, 20071024, 20071029, 20071036, 20071041, 20071044, 20071047, 20071050, 20071053, 20071056, 20071059,
20071063, 20072002, 20072005, 20072009, 20072012, 20072016, 20072037, 20072038, 20072046, 20072050, 20072055, 20073002, 20073006, 20073013, 20073017, 20073022,
20073025, 20073074, 20074002, 20074006, 20074010, 20074011, 20074016, 20074019, 20074023, 20074030, 20074043, 20074046, 20075004, 20075008, 20075014, 20075010,
20075015, 20075020, 20075025, 20075031, 20075035, 20075039, 20075043, 20075048, 20075054, 20075057, 20075060, 20075066, 20076001, 20076004, 20076007, 20076010,
20076013, 20076017, 20076021, 20076025, 20076028, 20076031, 20076034, 20076037, 20076040, 20076045, 20076048, 20076051, 20076054, 20076059, 20077002, 20077005,
20077008, 20077011, 20077014, 20077017, 20077018, 20078001, 20078007, 20078013, 20078016, 20078019, 20078022, 20078028, 20078032, 20078035, 20078040, 20078043,
20078046, 20078051, 20078054, 20078057, 20078060, 20078063, 20078067, 20079006, 20079009, 20079013, 20079017, 20079020, 20079023, 20079044, 20080006, 20080009,
20080012, 20080016, 20080020, 20081001, 20081004, 20081007, 20081012, 20081015, 20081018, 20081025, 20082002, 20082005, 20082010, 20082013, 20082016, 20082019,
20082024, 20082029, 20082034, 20082038, 20082047, 20082050, 20082053, 20082056, 20082059, 20082063, 20082066, 20083001, 20083004, 20083019, 20083022, 20083025,
20083029, 20083032, 20083074, 20083077, 20083079, 20084001, 20084002, 20084005, 20084009, 20084013, 20084016, 20084022, 20085006, 20085009, 20085017, 20086002,
20086005, 20086056, 20086011, 20086015, 20087008, 20087012, 20087021, 20088005, 20088009, 20088012, 20088030, 20088033, 20088037, 20089003, 20089006, 20089009,
20089012, 20089015, 20089018, 20089028, 20090002, 20090005, 20090008, 20090011, 20090014, 20090017, 20090021, 20090031, 20090048, 20091003, 20091006, 20091009,
20091012, 20091016, 20091019, 20091020, 20092005, 20092012, 20092015, 20092018, 20092021, 20092024, 20092027, 20092030, 20092033, 20092038, 20092053, 20092057,
20093001, 20093005, 20093010, 20093016, 20093025, 20093035,

//9.2 GeV 2020
21036018, 21036022, 21036025, 21036027, 21036028, 21036033, 21037063, 21038005, 21038006, 21038007, 21038020, 21038021, 21038037, 21038043, 21038048, 21038050,
21039005, 21039016, 21039025, 21039029, 21039032, 21039035, 21040002, 21040003, 21040005, 21040006, 21056026, 21056032, 21057028, 21062015, 21062020, 21062021,
21064048, 21069006, 21069038, 21069042, 21070028, 21072016, 21073008, 21073042, 21077024, 21080019, 21169035, 21169039, 21171031, 21171032, 21171033, 21174050,
21175009, 21177019, 21177021, 21177022, 21177032, 21178013, 21179001, 21184025, 21186026, 21193009, 21193027, 21194002, 21198001, 21198008, 21203017, 21203018,
21203019, 21204001, 21204002, 21204003, 21204004, 21204006, 21204007, 21204008, 21204009, 21204010, 21204013, 21205002, 21205020, 21209006, 21209007, 21213020,
21218001, 21218013, 21221012, 21225035, 21225040, 21225041, 21225042, 21225043, 21231020, 21233008, 21235015, 21235035, 21236006, 21243038,

////Au+Au 7.7 Gev
22033014,
22035002,
22038017,
22038018,
22038019,
22039010,
22039013,
22039028,
22040030,
22040038,
22044004,
22044005,
22044006,
22046006,
22047008,
22050003,
22050040,
22050044,
22052032,
22052033,
22052035,
22052036,
22052048,
22053022,
22054007,
22054022,
22054028,
22054030,
22054031,
22054042,
22055023,
22056036,
22059005,
22061012,
22062034,
22062035,
22062036,
22065014,
22065015,
22067039,
22069033,
22069034,
22069040,
22074042,
22076033,
22078016,
22079027,
22080021,
22081025,
22081038,
22081041,
22083027,
22085009,
22085037,
22091018,
22091020,
22091021,
22091022,
22092027,
22093029,
22096010,
22096037,
22097016,
22097030,
22097037,
22098054,
22099024,
22099042,
22100045,
22101043,
22102034,
22103027,
22103033,
22104027,
22105030,
22106032,
22108050,
22109032,
22110007,
22111047,
22112021,
22113001,
22113004,
22113029,
22113031,
22114030,
22115004,
22115008,
22115019,
22115032,
22115033,
22116008,
22116025,
22116026,
22117021,
22117023,
22118058,
22120031
};
  // end