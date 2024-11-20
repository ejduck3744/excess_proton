TString dataset ="11p5_gev_2020";
// TString ep_setting = "flow";
// TString ep_setting = "flatten";
// TString ep_setting = "recenter";
double y_cm = 0.0;
double beam_energy= 11.5;
bool include_etof =0;
bool do_tof_efficiency = 0;
int do_eta_weighting = 0;
TString re_file_name = "";
TString shift_file_name = "";
TString eta_weight_name = "";
TString efficiency_tof_file_name = "";
int centFull[9]={10,19,31,49,74,107,151,212,253};
const int nrun=1373;
int numbers[1373]={
20365001,
20365002,
20365003,
20365004,
20365005,
20365006,
20365007,
20365008,
20365009,
20365010,
20365011,
20365012,
20365013,
20365014,
20365015,
20365016,
20365017,
20365018,
20365019,
20365020,
20365021,
20365022,
20365023,
20365024,
20365025,
20365028,
20365029,
20365030,
20365031,
20365032,
20365033,
20365034,
20365035,
20365036,
20365039,
20365040,
20365041,
21001001,
21001002,
21001022,
21001023,
21001024,
21001025,
21001026,
21001027,
21001028,
21001029,
21001031,
21001032,
21001033,
21001034,
21001035,
21001036,
21001037,
21001038,
21001039,
21001040,
21001041,
21002001,
21002002,
21002003,
21002006,
21002007,
21002009,
21002010,
21002011,
21002012,
21002013,
21002014,
21002018,
21002019,
21002020,
21002021,
21002022,
21002023,
21002024,
21002025,
21002026,
21002027,
21002028,
21002029,
21002030,
21002031,
21002032,
21002033,
21002034,
21002035,
21002036,
21002037,
21002038,
21002040,
21003001,
21003003,
21003004,
21003005,
21003006,
21003007,
21003008,
21003011,
21003014,
21003015,
21003016,
21003017,
21003044,
21003045,
21003046,
21003047,
21004001,
21004002,
21004003,
21004004,
21004005,
21004007,
21004008,
21004011,
21004013,
21004015,
21004016,
21004017,
21004018,
21004019,
21004020,
21004021,
21004022,
21004023,
21004024,
21004026,
21004027,
21004028,
21004029,
21004030,
21004033,
21004034,
21004035,
21004036,
21004037,
21004038,
21004039,
21004040,
21004041,
21004042,
21004043,
21004044,
21004045,
21004046,
21004047,
21004048,
21005001,
21005004,
21005006,
21005007,
21005008,
21005009,
21005010,
21005011,
21005012,
21005013,
21005014,
21005015,
21005016,
21005017,
21005018,
21005019,
21005020,
21005021,
21005022,
21005023,
21005024,
21005025,
21005026,
21005027,
21005028,
21005029,
21005032,
21005033,
21005034,
21005035,
21005036,
21005037,
21005038,
21005039,
21005040,
21005042,
21005043,
21005044,
21005045,
21006001,
21006002,
21006003,
21006004,
21006005,
21006006,
21006007,
21006008,
21006009,
21006010,
21006011,
21006012,
21006013,
21006014,
21006015,
21006026,
21006027,
21006028,
21006029,
21006031,
21007001,
21007002,
21007003,
21007004,
21007005,
21007006,
21007008,
21007010,
21007011,
21007012,
21007013,
21007014,
21007015,
21007016,
21007020,
21007021,
21007022,
21007023,
21007024,
21007025,
21007026,
21007027,
21007028,
21007029,
21007030,
21007031,
21007032,
21007033,
21007034,
21007035,
21007036,
21007037,
21007038,
21007039,
21007040,
21007041,
21007042,
21007043,
21007044,
21008001,
21008002,
21008003,
21008004,
21008005,
21008006,
21008007,
21008008,
21008010,
21008011,
21008012,
21008014,
21008015,
21008016,
21008034,
21008042,
21008043,
21008044,
21008045,
21008046,
21008047,
21009001,
21009002,
21009003,
21009004,
21009005,
21009006,
21009007,
21009008,
21009009,
21009010,
21009011,
21009012,
21009013,
21009045,
21009046,
21009047,
21009048,
21009049,
21009050,
21009051,
21009052,
21009053,
21009054,
21010001,
21010002,
21010006,
21010007,
21010008,
21010009,
21010010,
21010011,
21010013,
21010014,
21010015,
21010016,
21010017,
21010018,
21010033,
21010034,
21010035,
21010036,
21011001,
21011002,
21011003,
21011004,
21011006,
21011007,
21011008,
21011009,
21011013,
21011015,
21011016,
21011017,
21011018,
21011021,
21011022,
21011023,
21011024,
21011036,
21011037,
21011039,
21011040,
21011041,
21011042,
21011044,
21011045,
21011047,
21011048,
21011049,
21012001,
21012002,
21012004,
21012011,
21012012,
21012013,
21012014,
21012015,
21012016,
21012017,
21012018,
21012019,
21012020,
21012021,
21012022,
21012023,
21012024,
21012025,
21012026,
21012027,
21012028,
21012029,
21012030,
21012031,
21012032,
21012033,
21012034,
21012035,
21012036,
21012037,
21012038,
21012039,
21013001,
21013002,
21013003,
21013004,
21013005,
21013006,
21013007,
21013008,
21013009,
21013011,
21013012,
21013013,
21013014,
21013015,
21013016,
21013020,
21013021,
21013022,
21013023,
21013024,
21014001,
21014002,
21014003,
21014005,
21014006,
21014007,
21014008,
21014009,
21014010,
21014011,
21014012,
21014014,
21014015,
21014016,
21014017,
21014018,
21014019,
21014022,
21014027,
21014029,
21014030,
21014031,
21014032,
21014033,
21014034,
21014035,
21014036,
21014039,
21014040,
21014041,
21014042,
21014043,
21014045,
21014046,
21014047,
21014048,
21014049,
21015001,
21015002,
21015003,
21015004,
21015005,
21015006,
21015007,
21015008,
21015009,
21015010,
21015011,
21015012,
21015013,
21015014,
21015015,
21015016,
21015017,
21015029,
21015030,
21015031,
21015032,
21015033,
21016001,
21016002,
21016003,
21016004,
21016005,
21016006,
21016007,
21016008,
21016009,
21016010,
21016011,
21016012,
21016013,
21016015,
21016024,
21016031,
21016032,
21016033,
21016034,
21016035,
21017001,
21017002,
21017003,
21017004,
21017005,
21017006,
21017007,
21017008,
21017009,
21017010,
21017011,
21017042,
21017043,
21017044,
21017045,
21017046,
21017047,
21017048,
21017050,
21017051,
21017052,
21017053,
21017054,
21017055,
21018001,
21018002,
21018003,
21018005,
21018031,
21018032,
21018033,
21018034,
21018035,
21018036,
21018037,
21018038,
21018039,
21018040,
21018041,
21019016,
21019017,
21019018,
21019019,
21019020,
21019021,
21019022,
21019023,
21019024,
21019031,
21019032,
21019033,
21019034,
21019035,
21019037,
21019038,
21019039,
21019040,
21019041,
21019042,
21019043,
21019044,
21019045,
21019046,
21019047,
21019048,
21019049,
21019050,
21019051,
21019052,
21019053,
21019054,
21019055,
21019056,
21019057,
21019058,
21019059,
21019060,
21019061,
21019062,
21019069,
21019072,
21019073,
21019074,
21019075,
21019076,
21019077,
21019078,
21020001,
21020002,
21020003,
21020004,
21020005,
21020006,
21020007,
21020013,
21020014,
21020015,
21020016,
21020017,
21020018,
21020019,
21020020,
21020023,
21020024,
21020026,
21020027,
21020037,
21020038,
21020039,
21020040,
21020041,
21020042,
21020043,
21020044,
21020045,
21020046,
21020047,
21020048,
21020049,
21020050,
21020051,
21020052,
21020053,
21020054,
21020055,
21020056,
21020057,
21020058,
21020059,
21020060,
21020061,
21020062,
21021001,
21021002,
21021008,
21021009,
21021010,
21021011,
21021012,
21021013,
21021014,
21021015,
21021016,
21021017,
21021019,
21021029,
21021030,
21021031,
21021032,
21021034,
21021042,
21021043,
21021044,
21021045,
21021046,
21021047,
21021048,
21021050,
21021051,
21021052,
21021053,
21021054,
21021055,
21021056,
21021057,
21021058,
21021059,
21021060,
21021062,
21021063,
21021064,
21021065,
21022005,
21022006,
21022007,
21022008,
21022009,
21022010,
21022011,
21022012,
21022013,
21022014,
21022015,
21022016,
21022017,
21022018,
21022019,
21022020,
21022053,
21022054,
21022055,
21023001,
21023002,
21023003,
21023004,
21023005,
21023006,
21023007,
21023008,
21023009,
21023010,
21023011,
21023012,
21023013,
21023014,
21023015,
21023016,
21023017,
21023018,
21023019,
21023026,
21023028,
21023029,
21023030,
21023031,
21024027,
21024028,
21024029,
21024030,
21024031,
21024032,
21024033,
21024034,
21024036,
21024037,
21024039,
21024040,
21024041,
21024043,
21024044,
21025002,
21025003,
21025004,
21025005,
21025006,
21025007,
21025008,
21025010,
21025011,
21025012,
21025013,
21025014,
21025015,
21025016,
21025017,
21025018,
21025019,
21025021,
21025022,
21025024,
21025026,
21025027,
21025028,
21025029,
21025030,
21025031,
21025032,
21025042,
21025043,
21025044,
21025046,
21025047,
21025048,
21025049,
21025050,
21025052,
21025054,
21025055,
21025056,
21025057,
21026001,
21026002,
21026003,
21026009,
21026010,
21026011,
21026012,
21026013,
21026014,
21026015,
21026016,
21026017,
21026018,
21026019,
21026020,
21026023,
21026024,
21026025,
21026026,
21026027,
21026028,
21026029,
21026031,
21026032,
21026033,
21026034,
21026035,
21026036,
21026037,
21026038,
21026039,
21026040,
21026041,
21026043,
21026044,
21026045,
21026047,
21026048,
21026050,
21026051,
21026052,
21026054,
21026055,
21026056,
21026059,
21026060,
21026062,
21027001,
21027002,
21027003,
21027004,
21027007,
21027008,
21027009,
21027011,
21027012,
21027013,
21027014,
21027015,
21027016,
21027018,
21041025,
21041027,
21041028,
21041029,
21041030,
21041032,
21041033,
21041034,
21041035,
21041036,
21041037,
21041038,
21041039,
21041040,
21042001,
21042002,
21042004,
21042005,
21042006,
21042007,
21042008,
21042009,
21042010,
21042011,
21042012,
21042014,
21042015,
21042016,
21042017,
21042020,
21042022,
21042023,
21042025,
21042026,
21042029,
21042030,
21042031,
21042032,
21042033,
21042034,
21042035,
21042036,
21042037,
21042038,
21042039,
21042040,
21042041,
21042042,
21043001,
21043002,
21043003,
21043004,
21043005,
21043006,
21043007,
21043008,
21043009,
21043010,
21043011,
21043012,
21043013,
21043014,
21043015,
21043016,
21043017,
21043018,
21043019,
21043020,
21043021,
21043022,
21043023,
21043024,
21043025,
21043026,
21043027,
21043028,
21043032,
21043033,
21043034,
21043035,
21043036,
21043037,
21043038,
21043039,
21043040,
21043041,
21043042,
21043043,
21043044,
21043045,
21043048,
21044001,
21044002,
21044003,
21044004,
21044005,
21044006,
21044007,
21044008,
21044009,
21044010,
21044011,
21044012,
21044013,
21044014,
21044015,
21044016,
21044017,
21045014,
21045015,
21045016,
21045017,
21045018,
21045019,
21045020,
21045021,
21045023,
21045024,
21045025,
21045028,
21045029,
21045030,
21045031,
21045032,
21045033,
21045034,
21045038,
21045039,
21045040,
21045041,
21045042,
21045043,
21045044,
21046001,
21046002,
21046003,
21046005,
21046006,
21046007,
21046008,
21046010,
21046011,
21046012,
21046013,
21046014,
21046015,
21046016,
21046017,
21046018,
21046019,
21046020,
21046021,
21046022,
21046023,
21046024,
21046025,
21046026,
21046027,
21046028,
21046029,
21046030,
21046036,
21046037,
21046038,
21046039,
21046040,
21046041,
21046042,
21046043,
21046044,
21046045,
21046046,
21046047,
21046048,
21046049,
21047001,
21047002,
21047003,
21047004,
21047005,
21047006,
21047007,
21047008,
21047009,
21047010,
21047011,
21047012,
21047013,
21047014,
21047015,
21047016,
21047017,
21047018,
21047019,
21047020,
21047021,
21047022,
21047023,
21047024,
21047025,
21047026,
21047027,
21047028,
21047029,
21047030,
21047031,
21047032,
21047036,
21047037,
21047038,
21047039,
21047040,
21047041,
21047042,
21047043,
21047044,
21047045,
21047046,
21047047,
21047048,
21047049,
21047050,
21047051,
21048001,
21048002,
21048003,
21048004,
21048005,
21048007,
21048008,
21048009,
21048010,
21048011,
21048012,
21048013,
21048014,
21048015,
21048017,
21048026,
21048027,
21048028,
21048029,
21048030,
21048031,
21048032,
21048033,
21048034,
21048035,
21048036,
21048039,
21048041,
21048042,
21048043,
21048044,
21048045,
21048046,
21048047,
21048048,
21048049,
21048050,
21048051,
21048052,
21048053,
21048054,
21048055,
21048056,
21048057,
21048058,
21048059,
21048060,
21048061,
21049001,
21049002,
21049003,
21049004,
21049005,
21049006,
21049007,
21049008,
21049009,
21049010,
21049011,
21049012,
21049014,
21049015,
21049016,
21049017,
21049019,
21049020,
21049021,
21049022,
21049023,
21049024,
21049025,
21049026,
21049027,
21049028,
21049029,
21049030,
21049031,
21049032,
21049033,
21049034,
21049035,
21049036,
21049037,
21049038,
21049039,
21049040,
21049041,
21049042,
21049043,
21049044,
21049045,
21049046,
21049051,
21049052,
21049053,
21049054,
21049055,
21049056,
21050001,
21050002,
21050003,
21050004,
21050005,
21050006,
21050007,
21050008,
21050009,
21050010,
21050011,
21050012,
21050013,
21050014,
21050015,
21050016,
21050017,
21050018,
21050047,
21050048,
21050049,
21050050,
21050055,
21050056,
21050057,
21050058,
21050059,
21050060,
21050061,
21050062,
21050063,
21050064,
21050065,
21050066,
21050067,
21050068,
21051001,
21051002,
21051003,
21051004,
21051005,
21051006,
21051007,
21051008,
21051009,
21051010,
21051011,
21051013,
21051014,
21051016,
21051017,
21051018,
21051019,
21051020,
21051021,
21051022,
21051023,
21051024,
21051025,
21051026,
21051027,
21051028,
21051029,
21051030,
21051031,
21051032,
21051033,
21051034,
21051035,
21051036,
21051037,
21051038,
21051039,
21051040,
21051041,
21051042,
21051043,
21051044,
21051045,
21051046,
21051047,
21051048,
21051049,
21051050,
21051051,
21051052,
21052001,
21052002,
21052003,
21052004,
21052005,
21052006,
21052007,
21052008,
21052009,
21052010,
21052011,
21052012,
21052013,
21052014,
21052015,
21052016,
21052017,
21052018,
21052020,
21052021,
21052028,
21052029,
21052030,
21052031,
21052032,
21052033,
21052034,
21052035,
21052036,
21052037,
21052038,
21052039,
21052040,
21052041,
21052042,
21052043,
21052044,
21052047,
21052048,
21052049,
21052050,
21052051,
21052052,
21052053,
21052054,
21052055,
21053001,
21053002,
21053003,
21053005,
21053006,
21053007,
21053008,
21053009,
21053010,
21053011,
21053012,
21053013,
21053014,
21053015,
21053016,
21053017,
21053018,
21053019,
21053020,
21053021,
21053022,
21053023,
21053024,
21053025,
21053026,
21053027,
21053028,
21053029,
21053030,
21053032,
21053033,
21053034,
21053035,
21053036,
21053037,
21053038,
21053039,
21053040,
21053041,
21053042,
21053043,
21053044,
21053045,
21053046,
21053047,
21053049,
21053050,
21053051,
21053054,
21053055,
21053056,
21053057,
21053058,
21053059,
21053060,
21053061,
21053062,
21053063,
21053064,
21053065,
21054001,
21054002,
21054003,
21054004,
21054005,
21054006,
21054007,
21054008,
21054010,
21054011,
21054012,
21054013,
21054014,
21054015,
21054016,
21054017,
21054018,
21054019,
21054020,
21054021,
21054022,
21054023,
21054024,
21054025,
21054026,
21054027,
21054028,
21054029,
21054030,
21054031,
21054032,
21054036,
21054037,
21054039,
21054041,
21054043,
21054044,
21054045,
21054046,
21054047,
21054048,
21054049,
21054050,
21054051,
21054052,
21054053,
21054054,
21054055,
21055001,
21055002,
21055003,
21055004,
21055005,
21055006,
21055007,
21055008,
21055009,
21055010,
21055011,
21055012,
21055015,
21055016,
21055017
};
