#TO REPLICATE THE RESULTS OF NALIN AND YAJIMA (2022) (INCLUDING THE APPENDIX)

#THIS VERSION: 02/03/2022

library(sfcr)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(kableExtra)
library(grid)
library(wbstats)


#help("sfcr_set") #for help on this step (NOTE, SCROLL DOWN THE HELP PAGE TO FIND THE INDEX OF THE MANUAL OF THE PACKAGE)
eqs <- sfcr_set(
  #model creation,
  #(1) arterial flows,
  #income,
  #yus~ cus+ gus +xus- imus,
  #yla ~ cla+ gla+ xla -imla,
  #disposable income,
  ydus~(yus + fbus +rus[-1]*busus_d[-1]+ rdus[-1]*dus_d[-1] +cba_d[-1]+cgus+cgbus+cgcus)*(1-thetaus),
  ydla~(yla + fbla + rla[-1]*blala_d[-1] + rdla[-1]*dla_d[-1]+der_d[-1]+cgla+cgbla+cgcla)*(1-thetala), 
  #disposable income (hs),
  #ydus_hs~ydus+ d(xrla)*busla_s[-1],
  #ydla_hs~ydla+ d(xrus)*blaus_s[-1],
  #wealth,
  vus ~ vus[-1]+  ydus-cus,
  vla ~ vla[-1]+ ydla-cla,
  #capital gains,
  cgus ~ (pcba-pcba[-1])*cba_s[-1],
  cgla ~ (pder-pder[-1])*der_s[-1],
  #tax,
  tla ~ thetala*(yla+ fbla+rla[-1]*blala_d[-1]+ rdla[-1]*dla_d[-1]+der_d[-1]),
  tus ~ thetaus*(yus+fbus +rus[-1]*busus_d[-1]+rdus[-1]*dus_d[-1] + cba_d[-1]) ,
  #cb profits,
  fcbla ~ rla[-1]*bcbla_d[-1] + rus[-1]*bcblaus_s[-1]*xrus,
  fcbus ~ rus[-1]*bcbus_d[-1], 
  #government budget constraint,
  bla_s ~ bla_s[-1]+gla-tla+ rla[-1]*bla_s[-1] - fcbla+cgbus, 
  bus_s ~ bus_s[-1]+gus-tus+ rus[-1]*bus_s[-1] - fcbus+cgbla,
  #current and capital account,
  cabla~xla - imla - rla[-1]*bbusla_s[-1] + rus[-1]*bblaus_s[-1]*xrus - rcusla[-1]*bcusla_s[-1] + rclaus[-1]*bclaus_s[-1]*xrus + rus[-1]*bcblaus_s[-1]*xrus ,
  kabla ~ (bbusla_s-bbusla_s[-1]) +(bcusla_s-bcusla_s[-1]) - (bcblaus_s-bcblaus_s[-1])*xrus- (bblaus_s-bblaus_s[-1])*xrus - (bclaus_s-bclaus_s[-1])*xrus,
  cabus~xus - imus + rla[-1]*bbusla_s[-1]*xrla - rus[-1]*bblaus_s[-1] + rcusla[-1]*bcusla_s[-1]*xrla - rclaus[-1]*bclaus_s[-1] - rus[-1]*bcblaus_s[-1],
  #kabus ~ -d(bbla_s)*xrla + d(blaus_s) + d(bcblaus_s),
  #(2) trade,
  #export and import,
  pmla~ exp(mi_0la+  mi_1la *log(r_pusy)  + (1- mi_1la)*log(r_play) - mi_1la *log(xrla)),
  pxla ~exp(chi_0la+  chi_1la *log(r_pusy)  + (1- chi_1la)*log(r_play)- chi_1la *log(xrla)),
  pxus ~ pmla*xrla,  
  pmus ~ pxla*xrla,
  r_xla~ exp(e_la - eta_la *log(pmus/r_pusy) + epsilon_la*log(r_yus)),
  r_imla ~ exp(p_la - psi_la *log(pmla[-1]/r_play[-1]) + pi_la*log(r_yla)),
  r_xus  ~r_imla,
  r_imus~r_xla,
  xla ~ r_xla*pxla,
  xus ~ r_xus*pxus,
  imla ~ r_imla*pmla,  
  imus ~ r_imus*pmus,
  #(3) income and expenditure,
  #real disposable income is of the haig-simon type,
  r_vla ~ vla/r_pdsla, 
  r_vus~ vus/r_pdsus,
  r_ydla ~ ydla/r_pdsla - r_vla[-1]* (r_pdsla-r_pdsla[-1])/r_pdsla ,
  r_ydus ~ ydus/r_pdsus - r_vus[-1]* (r_pdsus-r_pdsus[-1])/r_pdsus ,
  r_cla ~ alpha_1la*r_ydlae + alpha_2la*r_vla[-1],
  r_cus ~alpha_1us*r_yduse + alpha_2us*r_vus[-1],
  r_ydlae~ (r_ydla + r_ydla[-1])/2,
  r_yduse~ (r_ydus + r_ydus[-1])/2,
  r_sla ~r_cla+r_gla+r_xla ,
  r_sus ~r_cus+r_gus+r_xus ,
  sla ~r_sla*r_plas,
  sus ~r_sus*r_puss,
  r_plas ~((1+phila)*(wla*nla+imla+rcusla[-1]*bcusla_s[-1]*xrla)) /r_sla,
  r_puss ~((1+phius)*(wus*nus+imus+rclaus[-1]*bclaus_s[-1]*xrus)) / r_sus,
  r_pdsla ~(sla-xla) / (r_sla-r_xla),
  r_pdsus ~(sus-xus) / (r_sus-r_xus),
  dsla ~sla-xla ,
  dsus ~sus-xus ,
  r_dsla ~r_cla+r_gla ,
  r_dsus ~r_cus+r_gus  ,
  yla~sla-imla,
  yus ~sus-imus,
  r_yla ~r_sla-r_imla ,
  r_yus ~r_sus-r_imus ,
  r_play ~ yla /r_yla ,
  r_pusy ~yus/ r_yus,
  cla ~r_cla*r_pdsla ,
  cus ~r_cus*r_pdsus  ,
  gla ~r_gla*r_pdsla ,
  gus ~r_gus*r_pdsus ,
  nla_t ~ r_yla/ r_prla , 
  nus_t ~ r_yus/ r_prus ,
  #(4) financial intermediaries,
  #us financial sector,
  bbus_d~dus_s-bbusla_d-bcusla_d-hbus_d+pcba*cba_s+aus_d-lusus_s,
  #d(vbus)~fbus,
  bbusla_d~rho_2us*cba_s,
  hbus_d~rho_0us*dus_s,
  rdus~ rdus[-1] + rho_1us*(rus-rus[-1]),
  cgbus ~ (xrla-xrla[-1])*bbusla_s[-1],
  cgcus ~ (xrla-xrla[-1])*bcusla_s[-1],
  pcba~(1/rcba)+log(pxla*xrla),
  #d(rcba)~rho_0int*,
  #d(rlus)~rho_2us*d(rus),
  rcusla~(rho_3us+1)*rus,
  fbus~rus[-1]*bbus_d[-1]-rdus*dus_s[-1]+rla[-1]*bbusla_d[-1]*xrla-cba_d[-1]+rcusla[-1]*bcusla_s[-1]*xrla+rlus[-1]*lusus_s[-1],
  #la financial sector,
  bbla_d~dla_s-bblaus_d-bclaus_d-hbla_d+pder*der_s+ala_d-llala_s,
  #d(vbla)~fbla,
  bblaus_d~rho_2la*der_s,
  hbla_d~rho_0la*dla_s,
  rdla ~ rdla[-1]+ rho_1la*(rla-rla[-1]),
  cgbla ~ (xrus-xrus[-1])*bblaus_s[-1],
  cgcla ~ (xrus-xrus[-1])*bclaus_s[-1],
  pder~1/rder,
  #d(rder)~-rho_2int*(rla-rus),
  #d(rlla)~rho_2la*d(rla),
  rclaus~(rho_3la+1)*rla,
  fbla~rla[-1]*bbla_d[-1]-rdla*dla_s[-1]+rus[-1]*bblaus_d[-1]*xrus -der_d[-1]+rclaus[-1]*bclaus_s[-1]*xrus+rlla[-1]*llala_s[-1],
  #intermediaries,
  #d(vint)~rus[-1]*bintus_d[-1]+rla[-1]*bintla_d[-1]*xrla +cgint-rcba[-1]*cba_d[-1],
  #cba_s~cba_d,
  #cgint~d(xrla)*bintla_s[-1],
  #pcba~1/rcba,
  #d(rcba)~rho_0int*d(pxla)-rho_1int*(rla-rus),
  #fint~rus[-1]*bintus_d[-1]+rla[-1]*bintla_d[-1]*xrla +cgint-rcba[-1]*cba_d[-1],
  #bintla_d~cba_s-bintus_d,
  #(5) assets demand,
  #accumulation,
  #US
  # Wage bill
  fus~ yus-wus*nus-rclaus[-1]*bclaus_s[-1]*xrus -daus-cgcla-rlus[-1]*lusus_d[-1],
  influs ~ r_pusy/r_pusy[-1]-1,
  wus ~ omega0us + wus[-1]*(1+inflwus),
  inflwus ~  thetapus*(influs[-1]),
  phius ~ phi0us + phius[-1]*(1+inflphius),
  inflphius ~  thetaphius*(influs[-1]),
  nus ~ nus[-1] + etanus*(nus_t - nus[-1]),
  
  # thetapus=nus/nus_t
  
  # Fdus=rho_1int*Fus[-1]
  
  # Fuus=Fus-Fdus
  
  bclaus_s~lus_d-lusus_d,
  
  lusus_d~0.1*lus_d,
  
  
  # Demand for bank loans
  lus_d~ bclaus_s[-1]+ lusus_d[-1] + ius_d - daus-fus+cgcla,
  
  # Accumuustion of capital
  kus ~ kus[-1] + ius_d - daus,
  
  # Depreciation allowances
  daus ~ deltaus*kus[-1],
  
  # Capital stock target 
  kus_t ~ kappaus*yus[-1],
  kappaus ~0.99999+ 0.0001*(lusus_d[-1]+bclaus_s[-1]*xrus/kus[-1]),
  
  # Demand for investment goods
  ius_d~gamma0us +  gammaus*(kus_t - kus[-1]) + daus,
  
  #LA
  # Wage bill
  fla ~ yla -wla*nla - rcusla[-1]*bcusla_s[-1]*xrla-dala-cgcus-rlla[-1]*llala_s[-1],
  inflla ~ r_play/r_play[-1]-1,
  wla~omega0la + wla[-1]*(1+inflwla),
  inflwla ~  thetapla*(inflla[-1]),
  phila ~ phi0la + phila[-1]*(1+inflphila),
  inflphila ~  thetaphila*(inflla[-1]),
  nla ~ nla[-1] + etanla*(nla_t - nla[-1]),
  
  # thetapla=nla/nla_t
  
  # Fdla=rho_0int*Fla[-1]
  
  # Fula=Fla-Fdla
  
  bcusla_s~lla_d-llala_d,
  
  llala_d~0.5+rho_0int*lla_d,
  
  rho_0int~ if (rla>rus) {0} else {0.1},
  
  rho_1int~if (rla>rus) {0} else {1},
  
  # Demand for bank loans
  lla_d~ (bcusla_s[-1]+ rho_1int*llala_d[-1] + ila_d- dala -fla+cgcus),
  
  #(1-rho_1int)*bcusla_s~ (1-rho_1int)*(bcusla_s[-1]+ ila_d- dala -fla+cgcus),
  
  # Accumulation of capital
  kla ~ kla[-1] + ila_d - dala,
  
  # Depreciation allowances
  dala ~ deltala*kla[-1],
  
  # Capital stock target 
  kla_t ~ kappala*yla[-1],
  kappala ~ 0.99999+ 0.0001*(llala_d[-1]+bcusla_s[-1]*xrla/kla[-1]),
  
  # Demand for investment goods
  ila_d ~ gamma0la + gammala*(kla_t - kla[-1]) + dala,
  
  #asset demand for la resident,
  blala_d~vla*(lambda_10la+lambda_11la*rla-lambda_13la*rdla-lambda_14la*(rder+dxrlae)),
  dla_d~vla*(lambda_40la-lambda_41la*rla+lambda_43la*rdla-lambda_44la*(rder+dxrlae)),
  der_d~(vla/pder)*(lambda_50la-lambda_51la*rla-lambda_53la*rdla+lambda_54la*(rder+dxrlae)),
  #blaus_d~vla*(lambda_20la-lambda_21la*rla-lambda_23la*rdla+lambda_22la*(rus+dxruse)),
  #hla_d/ vla~lambda_30la - lambda_31la*(rus+dxruse) -lambda_32la*rla,
  #asset demand for us resident,
  busus_d~vus*(lambda_10us +lambda_11us*rus-lambda_13us*rdus-lambda_14us*(rcba+dxruse)),
  dus_d~vus*(lambda_40us-lambda_41us*rus+lambda_43us*rdus-lambda_44us*(rcba+dxruse)),
  cba_d~(vus/pcba)*(lambda_50us-lambda_51us*rus-lambda_53us*rdus+lambda_54us*(rcba+dxruse)),
  #busla_d~vus*(lambda_20us -lambda_21us*rus+lambda_22us*(rla+dxrlae)),
  #hus_d/ vus~lambda_30us -lambda_31us*rus -lambda_32us* (rla+dxrlae),
  #asset demand for intermediaries,
  #bintus_d~vint*(lambda_10int + lambda_11int*rus - lambda_12int*(rla+dxrlae)),
  #bintla_d~vint*(lambda_20int - lambda_21int*rus + lambda_22int*(rla+dxrlae)),
  #expected change in the exchange rate (expectations to update),
  dxrlae~d(pder)/ pder,
  dxruse~d(pcba)/ pcba,
  #(6) assets supply,
  #demand for cash,
  hus_d~ vus - busus_d- dus_d- pcba*cba_d ,
  hla_d~ vla - blala_d-dla_d -pder*der_d,
  #cb demand for b, h; d and cba,
  hus_s~ hus_d,
  hla_s~ hla_d,
  dus_s~dus_d,
  dla_s~dla_d,
  busus_s~ busus_d,
  blala_s ~blala_d,
  cba_s~cba_d,
  der_s~der_d,
  bbus_s~ bbus_d,
  bbla_s ~bbla_d,
  aus_s ~aus_d,
  ala_s ~ala_d,
  bcusla_d ~ bclaus_s*xrus,
  bclaus_d ~ bcusla_s*xrla,
  #bintus_s~bintus_d,
  lusus_s~lusus_d,
  llala_s~llala_d,
  hbus_s~hbus_d,
  hbla_s~hbla_d,
  #supply of domestic t bills to cb,
  bcbus_s~ bcbus_d,
  bcbla_s~ bcbla_d,
  bcbus_d ~ bcbus_d[-1]+ (hus_s-hus_s[-1])+(hbus_s-hbus_s[-1])-(aus_s -aus_s[-1]),
  bcbla_d ~ bcbla_d[-1]+ (hla_s-hla_s[-1])+(hbla_s-hbla_s[-1])-(bcblaus_s-bcblaus_s[-1])*xrus-(ala_s -ala_s[-1]),
  #supply of assets abroad ,
  #busla_s~ bla_s- blala_s- bcbla_s,
  #exchange rate,
  #xrus~ blaus_d /blaus_s,
  xrla ~ 1/xrus,
  bblaus_s ~ bblaus_d*xrla,
  bcblaus_d ~ bcblaus_s*xrus,
  xrus ~ (bbusla_s)/(bbusla_d),
  bbusla_s ~ bla_s - blala_s- bcbla_s - bbla_s,
  #final equation (not possible to appear twice),
  #blaus_s~ blaus_d /xrus,
  rerla~(r_play/r_pusy)*(xrla),
  rerus~1/rerla,
  tbla~xla - imla,
  tbus~xus - imus,
  psbrla~d(bla_s),
  psbrus~d(bus_s),
  prbla~cabla+d(bla_s),
  prbus~cabus+d(bus_s),
  nwla~((vla)-(bla_s) +(bbla_d-dla_s+bblaus_d+hbla_d-pder*der_s)+(bcblaus_d+bcbla_d-hla_s)),
  nwus~((vus)-(bus_s) +(bbus_d-dus_s+bbusla_d+hbus_d-pcba*cba_s)+(bcbus_d-hus_s)),
  rsh_cla~r_cla/yla,
  sh_cabla~cabla/yla,
  sh_tbla~tbla/yla,
  sh_govdef~-psbrla/yla,
  sh_prbla~prbla/yla,
  nip~rla[-1]*bla_s[-1]/yla,
  govdeb~bla_s/yla,
  sh_bbusla_s~bbusla_s/bla_s,
  sh_cba~pcba*cba_d/r_vus,
  totla~pxla/pmla,
  nwtla~kla-bcusla_s,
  ucsla~(wla*nla+imla)/r_sla
)

external <- sfcr_set(
  alpha_1la ~ 0.9,
  alpha_1us ~ 0.6,
  alpha_2la ~  0.053333276,
  alpha_2us ~  0.213326723889398,
  #vertical and horizontal constraints respected
  lambda_10la ~ 0.4,
  lambda_11la ~ 0.5,
  lambda_12la ~ 0.5,
  lambda_13la~ 0.25,
  lambda_14la~ 0.25,
  lambda_20la ~ 0.3,
  lambda_21la ~ 0.5,
  lambda_22la ~ 0.5,
  lambda_23la~ 0.5,
  lambda_40la ~ 0.4,
  lambda_41la ~ 0.25,
  lambda_42la ~ 5,
  lambda_43la~ 0.5,
  lambda_44la~ 0.25,
  lambda_50la ~ 0.1,
  lambda_51la ~ 0.25,
  lambda_52la ~ 5,
  lambda_53la~ 0.25,
  lambda_54la~ 0.5,
  #vertical and horizontal constraints respected
  lambda_24la~ 5,
  lambda_30la ~ 0.25,
  lambda_31la ~ 5,
  lambda_32la ~ 5,
  lambda_33la~ 5,
  lambda_34la~ 5,
  #vertical and horizontal constraints respected
  lambda_10us ~ 0.4,
  lambda_11us ~ 0.5,
  lambda_12us ~ 0.5,
  lambda_13us~ 0.25,
  lambda_14us~ 0.25,
  lambda_20us ~ 0.3,
  lambda_21us ~ 0.5,
  lambda_22us ~ 0.5,
  lambda_23us~ 0.5,
  lambda_40us ~ 0.4,
  lambda_41us ~ 0.25,
  lambda_42us ~ 5,
  lambda_43us~ 0.5,
  lambda_44us~ 0.25,
  lambda_50us ~ 0.1,
  lambda_51us ~ 0.25,
  lambda_52us ~ 5,
  lambda_53us~ 0.25,
  lambda_54us~ 0.5,
  #vertical and horizontal constraints respected
  lambda_24us~ 5,
  lambda_30us ~ 0.25,
  lambda_31us ~ 5,
  lambda_32us ~ 5,
  lambda_33us~ 5,
  lambda_34us~ 5,
  #vertical and horizontal constraints respected
  lambda_60us~0.5,
  lambda_61us~0.5,
  lambda_62us~0.5,
  lambda_70us~0.5,
  lambda_71us~0.5,
  lambda_72us~0.5,
  lambda_60la~0.5,
  lambda_61la~0.5,
  lambda_62la~0.5,
  lambda_70la~0.5,
  lambda_71la~0.5,
  lambda_72la~0.5,
  #vertical and horizontal constraints respected
  e_la ~  - 1.18,
  eta_la ~ 0.7,
  epsilon_la ~ 0.8,
  p_la ~ - 3.015,
  psi_la ~ 0.7,
  pi_la ~  1.2,
  mi_0la ~ - 0.00001,
  chi_0la ~ - 0.00001,
  mi_1la ~ 0.5,
  chi_1la ~ 0.49,
  thetala ~ 0.2,
  thetaus ~ 0.2,
  thetabla ~ 0,
  thetabus ~ 0,
  rho_0us~0.03,
  rho_1us~0.9,
  rho_2us~2.7,
  rho_3us~0.03,
  rho_0la~0.03,
  rho_1la~0.9,
  rho_2la~2.7,
  rho_3la~0.03,
  rho_2int~0.9,
  deltala ~ 0.2158,
  gammala ~ 1,
  gamma0la ~ 0,
  #kappala~1,
  deltaus ~ 0.2158,
  gammaus ~ 1,
  gamma0us ~ 0,
  #kappaus ~ 1,
  aus_d ~ 0,
  ala_d ~ 0,
  #lla_d~ 2.5,
  
  # exogenous variables
  bcblaus_s ~ 0.02031,
  r_gla ~ 16,
  r_gus ~ 16,
  r_prla ~ 0.8,
  r_prus ~ 1.3333,
  rla ~ 0.03,
  rus ~ 0.03,
  rlla ~ 0.03,
  rlus ~ 0.03,
  rcba~0.16,
  rder~0.16,
  thetapla ~ 0,
  thetapus ~ 0,
  omega0la ~ 0,
  omega0us ~ 0,
  thetaphius ~ 0,
  thetaphila ~ 0,
  phi0la ~ 0,
  phi0us ~ 0,
  etanla ~ 1,
  etanus ~ 1
  
)



initial <- sfcr_set(
  # starting values for stocks
  bcblaus_d ~ 0.02031,
  bla_s ~ 145.954295,
  blala_d ~ 102.18,
  blala_s ~ 102.18,
  bbus_d ~ 13.250255,
  bbus_s ~ 13.250255,
  bus_s ~ 146.005715,
  busus_d ~ 102.19,
  busus_s ~ 102.19,
  bcbla_d ~ 17.27839,
  bcbla_s ~ 17.27839,
  bcbus_d ~ 17.2995,
  bcbus_s ~ 17.2995,
  hla_d ~ 7.2987,
  hla_s ~ 7.2987,
  hus_d ~ 7.2995,
  hus_s ~ 7.2995,
  bblaus_d ~ 13.24565,
  bblaus_s ~ 13.24565,
  bbusla_d ~ 13.250255,
  bbusla_s ~ 13.250255,
  cba_d ~ 3.9178,
  cba_s ~ 3.9178,
  der_d ~ 3.9178,
  der_s ~ 3.9178,
  bbla_d ~ 13.24565,
  bbla_s ~ 13.24565,
  hbus_d~10,
  hbus_s~10,
  hbla_d~10,
  hbla_s~10,
  dus_s~12.01426,
  dla_s~12.01426,
  dus_d~12.01426,
  dla_d~12.01426,
  aus_s ~ 0,
  ala_s ~ 0,
  r_vla ~ 152.62,
  r_vus ~ 152.63,
  vla ~ 145.97921,
  vus ~ 145.99001,
  ila_d ~ 20,
  ius_d ~ 20,
  kla_t ~ 91,
  kus_t ~ 91,
  bcusla_d ~ 1.25,
  bcusla_s ~ 1.25,
  kla ~ 91,
  dala ~ 20,
  bclaus_d ~ 1.25,
  bclaus_s ~ 1.25,
  lla_d~ 2.5,
  llala_d~ 1.25,
  llala_s~ 1.25,
  lus_d~ 2.5,
  lusus_d~ 1.25,
  lusus_s~ 1.25,
  kus ~ 91,
  daus ~ 20,
  rho_0int~0.1,
  rho_1int~1,
  # other endogenous
  r_cla ~ 81.393,
  r_cus ~ 81.401,
  cabla ~ 0,
  cabus ~ 0,
  cla ~ 77.851,
  cus ~ 77.86,
  r_dsla ~ 97.393,
  r_dsus ~ 97.401,
  dsla ~ 93.154,
  dsus ~ 93.164,
  #dxrlae ~ 0,
  fcbla ~ 0.00869,
  fcbus ~ 0.00895,
  fbus~0.6,
  fbla~0.8,
  gla ~ 15.304,
  gus ~ 15.304,
  r_imla ~ 11.928,
  r_imus ~ 11.926,
  imla ~ 11.407,
  imus ~ 11.409,
  kabla ~ 0.00002,
  #kabus ~ - 0.00002,
  nla ~ 73.046,
  nus ~ 73.054,
  nla_t ~ 73.046, 
  nus_t ~ 73.054,
  r_pdsus ~ 0.95648,
  r_pdsla ~ 0.95649,
  pmla ~ 0.95628,
  pmus ~ 0.95661,
  r_plas ~ 0.95646,
  r_puss ~ 0.9565,
  pxla ~ 0.95634,
  pxus ~ 0.95656,
  r_play ~ 0.95648,
  r_pusy ~ 0.95649,
  r_sla ~ 109.32,
  r_sus ~ 109.33,
  sla ~ 104.56,
  sus ~ 104.57,
  tla ~ 19.463,
  tus ~ 19.465,
  r_xla ~ 11.926,
  r_xus ~ 11.928,
  xla ~ 11.406,
  xus ~ 11.41,
  xrla ~ 1.0003,
  xrus ~ 0.99971,
  #xrlae ~ 1.0003,
  #xruse ~ 0.99971,
  r_yla ~ 97.392,
  r_yus ~ 97.403,
  yla ~ 93.154,
  yus ~ 93.164,
  ydla ~ 77.851,
  ydus ~ 77.86,
  r_ydla ~ 81.394,
  r_ydus ~ 81.402,
  r_ydlae ~ 81.394,
  r_yduse ~ 81.402,
  phila ~ 0.2381,
  phius ~ 0.2381,
  rerus~0,
  rerla~0,
  rdus~0.01,
  rdla~0.01,
  rcusla~0.0315,
  rclaus~0.0315,
  pcba~ 6.25,
  pder~ 6.25,
  dxrlae~0,
  dxruse~0,
  tbla~0,
  tbus~0,
  psbrla~0,
  psbrus~0,
  prbla~0,
  prbus~0,
  nwla~0,
  nwus~0,
  rsh_cla~0,
  sh_cabla~0,
  sh_tbla~0,
  sh_govdef~0,
  sh_prbla~0,
  nip~0,
  govdeb~0,
  sh_bbusla_s~0,
  sh_cba~0,
  totla~0,
  inflwla ~ 0,
  inflla ~ 0,
  inflwus ~ 0,
  influs ~ 0,
  inflphila ~ 0,
  inflphius~ 0,
  wla ~ 0.6,
  wus ~ 1,
  fla ~ 0,
  fus ~ 0,
  cgus~ 0,
  cgbus~ 0,
  cgcus~ 0,
  cgla~ 0,
  cgbla~ 0,
  cgcla~ 0,
  kappala~1,
  kappaus ~ 1
  
)

rer <- sfcr_baseline(
  equations = eqs, 
  external = external,
  initial=initial,
  periods = 200,
  max_iter = 5000,
  
  
)

#help("sfcr_shock")
#EXPORT SHOCK
shock1 <- sfcr_shock(
  variables = sfcr_set(
    e_la ~ - 1.08
  ),
  start = 100,
  end = 112
  
)

shock2 <- sfcr_shock(
  variables = sfcr_set(
    rus ~ 0.02
  ),
  start = 106,
  end = 113
  
)

shock3 <- sfcr_shock(
  variables = sfcr_set(
    e_la ~ - 1.18
  ),
  start = 112,
  end = 200
  
)

shock4 <- sfcr_shock(
  variables = sfcr_set(
    rus ~ 0.03
  ),
  start = 113,
  end = 200
  
)

scenario1 <- sfcr_scenario(
  baseline = rer,
  scenario = list(shock1, shock2, shock3, shock4),
  periods = 200
)

# To replicate the Figure 1 (Resolution 1000*1000)

p1 <- scenario1 %>%
  ggplot( aes(x=period, y=((rerla-mean(rerla))/sd(rerla)))) +
  geom_line( aes(linetype="1.1: RER Mexico"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_linetype_manual(values=c("solid"))+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  xlim(90, 200)+
  ylim(-2, 5)

p2 <- scenario1 %>%
  ggplot( aes(x=period, y=((totla-mean(totla))/sd(totla)))) +
  geom_line( aes(linetype="1.2: TOT Mexico"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_linetype_manual(values=c("solid"))+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  xlim(90, 200)+
  ylim(-2, 5)

p3 <- scenario1 %>%
  ggplot( aes(x=period, y=((pcba-mean(pcba))/sd(pcba)))) +
  geom_line( aes(linetype="1.3: Price of CLN"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_linetype_manual(values=c("solid"))+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  xlim(90, 200)+
  ylim(-2, 5)


p4 <- scenario1 %>%
  ggplot( aes(x=period, y=((cba_s-mean(cba_s))/sd(cba_s)))) +
  geom_line( aes(linetype="1.4: Stock of CLN"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_linetype_manual(values=c("solid"))+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  xlim(90, 200)+
  ylim(-6, 6)

grid.arrange(p1, p2, p3, p4, nrow = 2)

# To replicate the Figure 2 (Resolution 1000*1000)

p15 <- scenario1 %>%
  ggplot( aes(x=period, y=((ila_d-mean(ila_d))/sd(ila_d)))) +
  geom_line( aes(linetype="2.1: Investment Mexico"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_linetype_manual(values=c("solid"))+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  xlim(90, 200)
p16 <- scenario1 %>%
  ggplot( aes(x=period, y=((ila_d/kla-mean(ila_d/kla))/sd(ila_d/kla)))) +
  geom_line( aes(linetype="2.2: Capital Accumulation Mexico"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_linetype_manual(values=c("solid"))+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  xlim(90, 200)

p17 <- scenario1 %>%
  ggplot( aes(x=period, y=((yla/kla-mean(yla/kla))/sd(yla/kla)))) +
  geom_line( aes(linetype="2.3: Capacity Utilization Mexico"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_linetype_manual(values=c("solid"))+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  xlim(90, 200)

p18 <- scenario1 %>%
  ggplot( aes(x=period, y=((kla-mean(kla))/sd(kla)))) +
  geom_line( aes(linetype="2.4: Capital Stock Mexico"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_linetype_manual(values=c("solid"))+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  xlim(90, 200)

grid.arrange(p15, p18, p16, p17, nrow = 2)



# To replicate the Figure 3 (Resolution 1000*1000)




p25 <- scenario1 %>%
  ggplot( aes(x=period, y=((bcusla_s*xrla-mean(bcusla_s*xrla))/sd(bcusla_s*xrla)))) +
  geom_line( aes(linetype="3.1: Balance Sheet Effect Mexico"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_linetype_manual(values=c("solid"))+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  xlim(90, 200)

p22 <- scenario1 %>%
  ggplot( aes(x=period, y=(((kla-bcusla_s-llala_d)-mean(kla-bcusla_s-llala_d))/sd(kla-bcusla_s-llala_d)))) +
  geom_line( aes(linetype="3.2: Net Worth Mexico"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_linetype_manual(values=c("solid"))+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  xlim(90, 200)+
  ylim(-4,2)

p24 <- scenario1 %>%
  ggplot( aes(x=period, y=((fla/kla-mean(fla/kla))/sd(fla/kla)))) +
  geom_line( aes(linetype="3.3: Profit over Capital Mexico"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_linetype_manual(values=c("solid"))+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  xlim(90, 200)

p23 <- scenario1 %>%
  ggplot( aes(x=period, y=((ucsla-mean(ucsla))/sd(ucsla)))) +
  geom_line( aes(linetype="3.4: Unit Costs Mexico"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_linetype_manual(values=c("solid"))+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  xlim(90, 200)



grid.arrange(p25, p22, p24, p23, nrow = 2)

# To replicate the Table 1

bs <- sfcr_matrix(
  columns = c("Households MX", "Firms MX", "Financial Sector MX", "Government MX", "Central bank MX", "FX", "Households RoW", "Firms RoW", "Financial Sector RoW", "Government RoW", "Central bank RoW", "sum"),
  codes = c("hMX", "fMX", "bMX", "gMX", "cbMX", "xr", "hRoW", "fRoW", "bRoW", "gRoW", "cbRoW", "s"),
  r1 = c("Deposits", hMX = "+D", bMX = "-D", xr = "xr", hRoW = "+D", bRoW = "-D", s = "0"),
  r2 = c("Bills MX", hMX = "+B_hMXMX", bMX = "+B_bMXMX", gMX = "-B_MX", cbMX = "+B_cbMX", xr = "xr", bRoW = "+B_hRoWMX", s = "0"),
  r3 = c("Bills RoW", bMX = "+B_bMXRoW", cbMX = "+B_cbMXRoW", xr = "xr", hRoW = "+B_hRoWRoW", bRoW = "+B_bRoWRoW", gRoW = "-B_RoW", cbRoW = "+B_cbRoW", s = "0"),
  r4 = c("Corporate Bills MX", bMX = "+B_cMXRoW", xr = "xr", fRoW = "-B_cMXRoW", s = "0"),
  r5 = c("Corporate Bills RoW", fMX = "-B_cRoWMX", xr = "xr", bRoW = "+B_cMXRoW", s = "0"),
  r6 = c("Loans MX", bMX = "+LMX", fMX = "-LMX", s = "0"),
  r7 = c("Loans RoW", fRoW = "-LRoW", bRoW = "+LRoW", s = "0"),
  r8 = c("Derivatives", hMX = "+pDER", bMX = "-pDER", xr = "xr", s = "0"),
  r9 = c("CLNs", xr = "xr", hRoW = "+pCLN", bRoW = "-pCLN", s = "0"),
  r10 = c("Capital", fMX = "+KMX", xr = "xr", fRoW = "+KRoW", s = "K"),
  r11 = c("HPM", hMX = "+H_h", bMX = "+H_b", cbMX = "-H", xr = "xr", hRoW = "+H_h", bRoW = "+H_b", cbRoW = "-H", s = "0"),
  r12 = c("Balance", hMX = "-V", bMX = "0", gMX = "-NW_g",  cbMX = "-NW_cb", hRoW = "-V", bRoW = "0", gRoW = "-NW_g",  cbRoW = "-NW_cb", s = "K")
)

#sfcr_matrix_display(bs, "bs")


# To replicate the Table 2

tfm <- sfcr_matrix(
  columns = c("Households MX", "Firms MX", "Financial Sector MX", "Government MX", "Central bank MX", "Households RoW", "Firms RoW", "Financial Sector RoW", "Government RoW", "Central bank RoW"),
  codes = c("hMX", "fMX", "bMX", "gMX", "cbMX", "hRoW", "fRoW", "bRoW", "gRoW", "cbRoW"),
  c("Exports", fMX = "+xla", fRoW = "-imus"),
  c("Imports", fMX = "-imla", fRoW = "+xus"),
  c("Consumption MX", hMX = "-cla", fMX = "+cla"),
  c("Govt. Expenditures MX", hMX = "+gla", gMX = "-gla"),
  c("Investment MX", fMX = "+ila_d", fMX = "-ila_d"),
  c("Firms Profit MX", fMX = "-fla", fMX = "+fla"),
  c("Taxes MX", hMX = "-tla", gMX = "+tla"),
  c("Income MX", hMX = "+yla", fMX = "-yla"),
  c("Banks Profits MX", hMX = "+fbla", bMX = "-fbla"),
  c("CB Profits MX", gMX = "+fcbla", cbMX = "-fcbla"),
  c("Consumption RoW",  hRoW = "-cus", fRoW = "+cus"),
  c("Govt. Expenditures RoW",  hRoW = "+gus", gRoW = "-gus"),
  c("Investment RoW", fRoW = "+ius_d", fRoW = "-ius_d"),
  c("Firms Profit RoW", fRoW = "-fus", fRoW = "+fus"),
  c("Taxes RoW",  fRoW = "-tus", gRoW = "+tus"),
  c("Income RoW",  hRoW = "+yus", fRoW = "-yus"),
  c("Banks Profits RoW",  hRoW = "+fbus", bRoW = "-fbus"),
  c("CB Profits RoW",  gRoW = "+fcbus", cbRoW = "-fcbus"),
  c("Int. Deposits MX", hMX = "+rdla[-1]*dla_d[-1]", bMX = "-rdla[-1]*dla_d[-1]"),
  c("Int. Bills MX", hMX = "+rla[-1]*blala_d[-1]", bMX = "+rla[-1]*bbla_d[-1]", gMX = "-rla[-1]*bla_s[-1]", cbMX = "+rla[-1]*bcbla_d[-1]", bRoW = "+rla[-1]*bbusla_d[-1]"),
  c("Int. Corporate Bills MX", fMX = "-rcusla[-1]*bcusla_d[-1]", bRoW = "+rcusla[-1]*bcusla_d[-1]"),
  c("Int. Loans MX", fMX = "-rlla[-1]*lla_d[-1]", bMX = "+rlla[-1]*lla_d[-1]"),
  c("Int. Derivatives", hMX = "+rder[-1]*der_d[-1]", bMX = "-rder[-1]*der_d[-1]"),
  c("Int. Deposits RoW", hRoW = "+rdus[-1]*dus_d[-1]", bRoW = "-rdus[-1]*dus_d[-1]"),
  c("Int. Bills RoW", bMX = "+rus[-1]*bblaus_d[-1]", cbMX = "+rus[-1]*bcblaus_d[-1]", hRoW = "+rus[-1]*busus_d[-1]", bRoW = "+rus[-1]*bbus_d[-1]", gRoW = "-rus[-1]*bus_s[-1]", cbRoW = "+rus[-1]*bcbus_d[-1]"),
  c("Int. Corporate Bills RoW", fRoW = "-rclaus[-1]*bclaus_d[-1]", bMX = "+rclaus[-1]*bclaus_d[-1]"),
  c("Int. Loans RoW", fRoW = "-rlus[-1]*lus_d[-1]", bRoW = "+rlus[-1]*lus_d[-1]"),
  c("Int. CLNs", hRoW = "+rcba[-1]*cba_d[-1]", bRoW = "-rcba[-1]*cba_d[-1]"),
  c("Ch. Deposits MX", hMX = "-(dla_s-dla_s[-1])", bMX = "+(dla_s-dla_s[-1])"),
  c("Ch. Bills MX", hMX = "-(blala_s-blala_s[-1])", bMX = "-(bbla_s-bbla_s[-1])", gMX = "+(bla_s-bla_s[-1])", cbMX = "-(bcbla_s-bcbla_s[-1])", bRoW = "-(bbusla_s-bbusla_s[-1])"),
  c("Ch. Corporate Bills MX", bRoW = "-(bcusla_s-bcusla_s[-1])", fMX = "+(bcusla_s-bcusla_s[-1])"),
  c("Ch. Loans MX", bMX = "-(lla_s-lla_s[-1])", fMX = "+(lla_s-lla_s[-1])"),
  c("Ch. Derivatives", hMX = "-(pder*der_s-pder[-1]*der_s[-1])", bMX = "+(pder*der_s-pder[-1]*der_s[-1])" ),
  c("Ch. HPM MX", hMX = "-(hla_s-hla_s[-1])", bMX = "-(hbla_s-hbla_s[-1])", cbMX = "+(hla_s-hla_s[-1])+(hbla_s-hbla_s[-1])"),
  c("Ch. Deposits RoW",  hRoW = "-(dus_s-dus_s[-1])", bRoW = "+(dus_s-dus_s[-1])"),
  c("Ch. Bills RoW", bMX = "-(bblaus_s-bblaus_s[-1])", cbMX = "-(bcblaus_s-bcblaus_s[-1])", hRoW = "-(busus_s-busus_s[-1])", bRoW = "-(bbus_s-bbus_s[-1])", gRoW = "+(bus_s-bus_s[-1])", cbRoW = "-(bcbus_s-bcbus_s[-1])"),
  c("Ch. Corporate Bills RoW", bRoW = "-(bclaus_s-bclaus_s[-1])", fMX = "+(bclaus_s-bclaus_s[-1])"),
  c("Ch. Loans RoW", bRoW = "-(lus_s-lus_s[-1])", fRoW = "+(lus_s-lus_s[-1])"),
  c("Ch. CLNs", hRoW = "-(pcba*cba_s-pcba[-1]*cba_s[-1])", bRoW = "+(pcba*cba_s-pcba[-1]*cba_s[-1])"),
  c("Ch. HPM RoW", hRoW = "-(hus_s-hus_s[-1])", bRoW = "-(hbus_s-hbus_s[-1])", cbRoW = "+(hus_s-hus_s[-1])+(hbus_s-hbus_s[-1])"),
)

#sfcr_matrix_display(tfm, "tfm")

#sfcr_dag_blocks(eqs)
#sfcr_dag_blocks_plot(eqs, title = NULL, size = 10)

#sfcr_dag_cycles(eqs)
#sfcr_dag_cycles_plot(eqs, title = NULL, size = 10)


# Robustness check


eqs1 <- sfcr_set(
  #model creation,
  #(1) arterial flows,
  #income,
  #yus~ cus+ gus +xus- imus,
  #yla ~ cla+ gla+ xla -imla,
  #disposable income,
  ydus~(yus + fbus +rus[-1]*busus_d[-1]+ rdus[-1]*dus_d[-1] +cba_d[-1]+cgus+cgbus+cgcus)*(1-thetaus),
  ydla~(yla + fbla + rla[-1]*blala_d[-1] + rdla[-1]*dla_d[-1]+der_d[-1]+cgla+cgbla+cgcla)*(1-thetala), 
  #disposable income (hs),
  #ydus_hs~ydus+ d(xrla)*busla_s[-1],
  #ydla_hs~ydla+ d(xrus)*blaus_s[-1],
  #wealth,
  vus ~ vus[-1]+  ydus-cus,
  vla ~ vla[-1]+ ydla-cla,
  #capital gains,
  cgus ~ (pcba-pcba[-1])*cba_s[-1],
  cgla ~ (pder-pder[-1])*der_s[-1],
  #tax,
  tla ~ thetala*(yla+ fbla+rla[-1]*blala_d[-1]+ rdla[-1]*dla_d[-1]+der_d[-1]),
  tus ~ thetaus*(yus+fbus +rus[-1]*busus_d[-1]+rdus[-1]*dus_d[-1] + cba_d[-1]) ,
  #cb profits,
  fcbla ~ rla[-1]*bcbla_d[-1] + rus[-1]*bcblaus_s[-1]*xrus,
  fcbus ~ rus[-1]*bcbus_d[-1], 
  #government budget constraint,
  bla_s ~ bla_s[-1]+gla-tla+ rla[-1]*bla_s[-1] - fcbla+cgbus, 
  bus_s ~ bus_s[-1]+gus-tus+ rus[-1]*bus_s[-1] - fcbus+cgbla,
  #current and capital account,
  cabla~xla - imla - rla[-1]*bbusla_s[-1] + rus[-1]*bblaus_s[-1]*xrus - rcusla[-1]*bcusla_s[-1] + rclaus[-1]*bclaus_s[-1]*xrus + rus[-1]*bcblaus_s[-1]*xrus ,
  kabla ~ (bbusla_s-bbusla_s[-1]) +(bcusla_s-bcusla_s[-1]) - (bcblaus_s-bcblaus_s[-1])*xrus- (bblaus_s-bblaus_s[-1])*xrus - (bclaus_s-bclaus_s[-1])*xrus,
  cabus~xus - imus + rla[-1]*bbusla_s[-1]*xrla - rus[-1]*bblaus_s[-1] + rcusla[-1]*bcusla_s[-1]*xrla - rclaus[-1]*bclaus_s[-1] - rus[-1]*bcblaus_s[-1],
  #kabus ~ -d(bbla_s)*xrla + d(blaus_s) + d(bcblaus_s),
  #(2) trade,
  #export and import,
  pmla~ exp(mi_0la+  mi_1la *log(r_pusy)  + (1- mi_1la)*log(r_play) - mi_1la *log(xrla)),
  pxla ~exp(chi_0la+  chi_1la *log(r_pusy)  + (1- chi_1la)*log(r_play)- chi_1la *log(xrla)),
  pxus ~ pmla*xrla,  
  pmus ~ pxla*xrla,
  r_xla~ exp(e_la - eta_la *log(pmus/r_pusy) + epsilon_la*log(r_yus)),
  r_imla ~ exp(p_la - psi_la *log(pmla[-1]/r_play[-1]) + pi_la*log(r_yla)),
  r_xus  ~r_imla,
  r_imus~r_xla,
  xla ~ r_xla*pxla,
  xus ~ r_xus*pxus,
  imla ~ r_imla*pmla,  
  imus ~ r_imus*pmus,
  #(3) income and expenditure,
  #real disposable income is of the haig-simon type,
  r_vla ~ vla/r_pdsla, 
  r_vus~ vus/r_pdsus,
  r_ydla ~ ydla/r_pdsla - r_vla[-1]* (r_pdsla-r_pdsla[-1])/r_pdsla ,
  r_ydus ~ ydus/r_pdsus - r_vus[-1]* (r_pdsus-r_pdsus[-1])/r_pdsus ,
  r_cla ~ alpha_1la*r_ydlae + alpha_2la*r_vla[-1],
  r_cus ~alpha_1us*r_yduse + alpha_2us*r_vus[-1],
  r_ydlae~ (r_ydla + r_ydla[-1])/2,
  r_yduse~ (r_ydus + r_ydus[-1])/2,
  r_sla ~r_cla+r_gla+r_xla ,
  r_sus ~r_cus+r_gus+r_xus ,
  sla ~r_sla*r_plas,
  sus ~r_sus*r_puss,
  r_plas ~((1+phila)*(wla*nla+imla+rcusla[-1]*bcusla_s[-1]*xrla)) /r_sla,
  r_puss ~((1+phius)*(wus*nus+imus+rclaus[-1]*bclaus_s[-1]*xrus)) / r_sus,
  r_pdsla ~(sla-xla) / (r_sla-r_xla),
  r_pdsus ~(sus-xus) / (r_sus-r_xus),
  dsla ~sla-xla ,
  dsus ~sus-xus ,
  r_dsla ~r_cla+r_gla ,
  r_dsus ~r_cus+r_gus  ,
  yla~sla-imla,
  yus ~sus-imus,
  r_yla ~r_sla-r_imla ,
  r_yus ~r_sus-r_imus ,
  r_play ~ yla /r_yla ,
  r_pusy ~yus/ r_yus,
  cla ~r_cla*r_pdsla ,
  cus ~r_cus*r_pdsus  ,
  gla ~r_gla*r_pdsla ,
  gus ~r_gus*r_pdsus ,
  nla_t ~ r_yla/ r_prla , 
  nus_t ~ r_yus/ r_prus ,
  #(4) financial intermediaries,
  #us financial sector,
  bbus_d~dus_s-bbusla_d-bcusla_d-hbus_d+pcba*cba_s+aus_d-lusus_s,
  #d(vbus)~fbus,
  bbusla_d~rho_2us*cba_s,
  hbus_d~rho_0us*dus_s,
  rdus~ rdus[-1] + rho_1us*(rus-rus[-1]),
  cgbus ~ (xrla-xrla[-1])*bbusla_s[-1],
  cgcus ~ (xrla-xrla[-1])*bcusla_s[-1],
  pcba~(1/rcba)+log(pxla*xrla),
  #d(rcba)~rho_0int*,
  #d(rlus)~rho_2us*d(rus),
  rcusla~(rho_3us+1)*rus,
  fbus~rus[-1]*bbus_d[-1]-rdus*dus_s[-1]+rla[-1]*bbusla_d[-1]*xrla-cba_d[-1]+rcusla[-1]*bcusla_s[-1]*xrla+rlus[-1]*lusus_s[-1],
  #la financial sector,
  bbla_d~dla_s-bblaus_d-bclaus_d-hbla_d+pder*der_s+ala_d-llala_s,
  #d(vbla)~fbla,
  bblaus_d~rho_2la*der_s,
  hbla_d~rho_0la*dla_s,
  rdla ~ rdla[-1]+ rho_1la*(rla-rla[-1]),
  cgbla ~ (xrus-xrus[-1])*bblaus_s[-1],
  cgcla ~ (xrus-xrus[-1])*bclaus_s[-1],
  pder~1/rder,
  #d(rder)~-rho_2int*(rla-rus),
  #d(rlla)~rho_2la*d(rla),
  rclaus~(rho_3la+1)*rla,
  fbla~rla[-1]*bbla_d[-1]-rdla*dla_s[-1]+rus[-1]*bblaus_d[-1]*xrus -der_d[-1]+rclaus[-1]*bclaus_s[-1]*xrus+rlla[-1]*llala_s[-1],
  #intermediaries,
  #d(vint)~rus[-1]*bintus_d[-1]+rla[-1]*bintla_d[-1]*xrla +cgint-rcba[-1]*cba_d[-1],
  #cba_s~cba_d,
  #cgint~d(xrla)*bintla_s[-1],
  #pcba~1/rcba,
  #d(rcba)~rho_0int*d(pxla)-rho_1int*(rla-rus),
  #fint~rus[-1]*bintus_d[-1]+rla[-1]*bintla_d[-1]*xrla +cgint-rcba[-1]*cba_d[-1],
  #bintla_d~cba_s-bintus_d,
  #(5) assets demand,
  #accumulation,
  #US
  # Wage bill
  fus~ yus-wus*nus-rclaus[-1]*bclaus_s[-1]*xrus -daus-cgcla-rlus[-1]*lusus_d[-1],
  influs ~ r_pusy/r_pusy[-1]-1,
  wus ~ omega0us + wus[-1]*(1+inflwus),
  inflwus ~  thetapus*(influs[-1]),
  phius ~ phi0us + phius[-1]*(1+inflphius),
  inflphius ~  thetaphius*(influs[-1]),
  nus ~ nus[-1] + etanus*(nus_t - nus[-1]),
  
  # thetapus=nus/nus_t
  
  # Fdus=rho_1int*Fus[-1]
  
  # Fuus=Fus-Fdus
  
  bclaus_s~lus_d-lusus_d,
  
  lusus_d~0.1*lus_d,
  
  
  # Demand for bank loans
  lus_d~ bclaus_s[-1]+ lusus_d[-1] + ius_d - daus-fus+cgcla,
  
  # Accumuustion of capital
  kus ~ kus[-1] + ius_d - daus,
  
  # Depreciation allowances
  daus ~ deltaus*kus[-1],
  
  # Capital stock target 
  kus_t ~ kappaus*yus[-1],
  kappaus ~0.99999+ 0.0001*(lusus_d[-1]+bclaus_s[-1]*xrus/kus[-1]),
  
  # Demand for investment goods
  ius_d~gamma0us +  gammaus*(kus_t - kus[-1]) + daus,
  
  #LA
  # Wage bill
  fla ~ yla -wla*nla - rcusla[-1]*bcusla_s[-1]*xrla-dala-cgcus-rlla[-1]*llala_s[-1],
  inflla ~ r_play/r_play[-1]-1,
  wla~omega0la + wla[-1]*(1+inflwla),
  inflwla ~  thetapla*(inflla[-1]),
  phila ~ phi0la + phila[-1]*(1+inflphila),
  inflphila ~  thetaphila*(inflla[-1]),
  nla ~ nla[-1] + etanla*(nla_t - nla[-1]),
  
  # thetapla=nla/nla_t
  
  # Fdla=rho_0int*Fla[-1]
  
  # Fula=Fla-Fdla
  
  bcusla_s~lla_d-llala_d,
  
  llala_d~rho_0int*lla_d,
  
  rho_0int~ 0.1,
  
  rho_1int~1,
  
  # Demand for bank loans
  lla_d~ (bcusla_s[-1]+ rho_1int*llala_d[-1] + ila_d- dala -fla+cgcus),
  
  #(1-rho_1int)*bcusla_s~ (1-rho_1int)*(bcusla_s[-1]+ ila_d- dala -fla+cgcus),
  
  # Accumulation of capital
  kla ~ kla[-1] + ila_d - dala,
  
  # Depreciation allowances
  dala ~ deltala*kla[-1],
  
  # Capital stock target 
  kla_t ~ kappala*yla[-1],
  kappala ~ 0.99999+ 0.0001*(llala_d[-1]+bcusla_s[-1]*xrla/kla[-1]),
  
  # Demand for investment goods
  ila_d ~ gamma0la + gammala*(kla_t - kla[-1]) + dala,
  
  #asset demand for la resident,
  blala_d~vla*(lambda_10la+lambda_11la*rla-lambda_13la*rdla-lambda_14la*(rder+dxrlae)),
  dla_d~vla*(lambda_40la-lambda_41la*rla+lambda_43la*rdla-lambda_44la*(rder+dxrlae)),
  der_d~(vla/pder)*(lambda_50la-lambda_51la*rla-lambda_53la*rdla+lambda_54la*(rder+dxrlae)),
  #blaus_d~vla*(lambda_20la-lambda_21la*rla-lambda_23la*rdla+lambda_22la*(rus+dxruse)),
  #hla_d/ vla~lambda_30la - lambda_31la*(rus+dxruse) -lambda_32la*rla,
  #asset demand for us resident,
  busus_d~vus*(lambda_10us +lambda_11us*rus-lambda_13us*rdus-lambda_14us*(rcba+dxruse)),
  dus_d~vus*(lambda_40us-lambda_41us*rus+lambda_43us*rdus-lambda_44us*(rcba+dxruse)),
  cba_d~(vus/pcba)*(lambda_50us-lambda_51us*rus-lambda_53us*rdus+lambda_54us*(rcba+dxruse)),
  #busla_d~vus*(lambda_20us -lambda_21us*rus+lambda_22us*(rla+dxrlae)),
  #hus_d/ vus~lambda_30us -lambda_31us*rus -lambda_32us* (rla+dxrlae),
  #asset demand for intermediaries,
  #bintus_d~vint*(lambda_10int + lambda_11int*rus - lambda_12int*(rla+dxrlae)),
  #bintla_d~vint*(lambda_20int - lambda_21int*rus + lambda_22int*(rla+dxrlae)),
  #expected change in the exchange rate (expectations to update),
  dxrlae~d(pder)/ pder,
  dxruse~d(pcba)/ pcba,
  #(6) assets supply,
  #demand for cash,
  hus_d~ vus - busus_d- dus_d- pcba*cba_d ,
  hla_d~ vla - blala_d-dla_d -pder*der_d,
  #cb demand for b, h; d and cba,
  hus_s~ hus_d,
  hla_s~ hla_d,
  dus_s~dus_d,
  dla_s~dla_d,
  busus_s~ busus_d,
  blala_s ~blala_d,
  cba_s~cba_d,
  der_s~der_d,
  bbus_s~ bbus_d,
  bbla_s ~bbla_d,
  aus_s ~aus_d,
  ala_s ~ala_d,
  bcusla_d ~ bclaus_s*xrus,
  bclaus_d ~ bcusla_s*xrla,
  #bintus_s~bintus_d,
  lusus_s~lusus_d,
  llala_s~llala_d,
  hbus_s~hbus_d,
  hbla_s~hbla_d,
  #supply of domestic t bills to cb,
  bcbus_s~ bcbus_d,
  bcbla_s~ bcbla_d,
  bcbus_d ~ bcbus_d[-1]+ (hus_s-hus_s[-1])+(hbus_s-hbus_s[-1])-(aus_s -aus_s[-1]),
  bcbla_d ~ bcbla_d[-1]+ (hla_s-hla_s[-1])+(hbla_s-hbla_s[-1])-(bcblaus_s-bcblaus_s[-1])*xrus-(ala_s -ala_s[-1]),
  #supply of assets abroad ,
  #busla_s~ bla_s- blala_s- bcbla_s,
  #exchange rate,
  #xrus~ blaus_d /blaus_s,
  xrla ~ 1/xrus,
  bblaus_s ~ bblaus_d*xrla,
  bcblaus_d ~ bcblaus_s*xrus,
  xrus ~ (bbusla_s)/(bbusla_d),
  bbusla_s ~ bla_s - blala_s- bcbla_s - bbla_s,
  #final equation (not possible to appear twice),
  #blaus_s~ blaus_d /xrus,
  rerla~(r_play/r_pusy)*(xrla),
  rerus~1/rerla,
  tbla~xla - imla,
  tbus~xus - imus,
  psbrla~d(bla_s),
  psbrus~d(bus_s),
  prbla~cabla+d(bla_s),
  prbus~cabus+d(bus_s),
  nwla~((vla)-(bla_s) +(bbla_d-dla_s+bblaus_d+hbla_d-pder*der_s)+(bcblaus_d+bcbla_d-hla_s)),
  nwus~((vus)-(bus_s) +(bbus_d-dus_s+bbusla_d+hbus_d-pcba*cba_s)+(bcbus_d-hus_s)),
  rsh_cla~r_cla/yla,
  sh_cabla~cabla/yla,
  sh_tbla~tbla/yla,
  sh_govdef~-psbrla/yla,
  sh_prbla~prbla/yla,
  nip~rla[-1]*bla_s[-1]/yla,
  govdeb~bla_s/yla,
  sh_bbusla_s~bbusla_s/bla_s,
  sh_cba~pcba*cba_d/r_vus,
  totla~pxla/pmla,
  nwtla~kla-bcusla_s,
  ucsla~(wla*nla+imla)/r_sla
)




rer1 <- sfcr_baseline(
  equations = eqs1, 
  external = external,
  initial=initial,
  periods = 200,
  max_iter = 5000,
  
  
)


scenario2 <- sfcr_scenario(
  baseline = rer1,
  scenario = list(shock1, shock2, shock3, shock4),
  periods = 200
)

# To replicate the Figure 1 (Resolution 1000*1000)

p1 <- scenario2 %>%
  ggplot( aes(x=period, y=((rerla-mean(rerla))/sd(rerla)))) +
  geom_line( aes(linetype="4.1: RER Mexico"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_linetype_manual(values=c("solid"))+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  xlim(90, 200)+
  ylim(-2, 5)

p2 <- scenario2 %>%
  ggplot( aes(x=period, y=((totla-mean(totla))/sd(totla)))) +
  geom_line( aes(linetype="4.2: TOT Mexico"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_linetype_manual(values=c("solid"))+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  xlim(90, 200)+
  ylim(-2, 5)

p3 <- scenario2 %>%
  ggplot( aes(x=period, y=((pcba-mean(pcba))/sd(pcba)))) +
  geom_line( aes(linetype="4.3: Price of CLN"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_linetype_manual(values=c("solid"))+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  xlim(90, 200)+
  ylim(-2, 5)


p4 <- scenario2 %>%
  ggplot( aes(x=period, y=((cba_s-mean(cba_s))/sd(cba_s)))) +
  geom_line( aes(linetype="4.4: Stock of CLN"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_linetype_manual(values=c("solid"))+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  xlim(90, 200)+
  ylim(-6, 6)

grid.arrange(p1, p2, p3, p4, nrow = 2)

# To replicate the Figure 2 (Resolution 1000*1000)

p15 <- scenario2 %>%
  ggplot( aes(x=period, y=((ila_d-mean(ila_d))/sd(ila_d)))) +
  geom_line( aes(linetype="5.1: Investment Mexico"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_linetype_manual(values=c("solid"))+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  xlim(90, 200)
p16 <- scenario2 %>%
  ggplot( aes(x=period, y=((ila_d/kla-mean(ila_d/kla))/sd(ila_d/kla)))) +
  geom_line( aes(linetype="5.2: Capital Accumulation Mexico"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_linetype_manual(values=c("solid"))+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  xlim(90, 200)

p17 <- scenario2 %>%
  ggplot( aes(x=period, y=((yla/kla-mean(yla/kla))/sd(yla/kla)))) +
  geom_line( aes(linetype="5.3: Capacity Utilization Mexico"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_linetype_manual(values=c("solid"))+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  xlim(90, 200)

p18 <- scenario2 %>%
  ggplot( aes(x=period, y=((kla-mean(kla))/sd(kla)))) +
  geom_line( aes(linetype="5.4: Capital Stock Mexico"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_linetype_manual(values=c("solid"))+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  xlim(90, 200)

grid.arrange(p15, p18, p16, p17, nrow = 2)



# To replicate the Figure 3 (Resolution 1000*1000)




p25 <- scenario2 %>%
  ggplot( aes(x=period, y=((bcusla_s*xrla-mean(bcusla_s*xrla))/sd(bcusla_s*xrla)))) +
  geom_line( aes(linetype="6.1: Balance Sheet Effect Mexico"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_linetype_manual(values=c("solid"))+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  xlim(90, 200)

p22 <- scenario2 %>%
  ggplot( aes(x=period, y=(((kla-bcusla_s-llala_d)-mean(kla-bcusla_s-llala_d))/sd(kla-bcusla_s-llala_d)))) +
  geom_line( aes(linetype="6.2: Net Worth Mexico"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_linetype_manual(values=c("solid"))+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  xlim(90, 200)+
  ylim(-4,2)

p24 <- scenario2 %>%
  ggplot( aes(x=period, y=((fla/kla-mean(fla/kla))/sd(fla/kla)))) +
  geom_line( aes(linetype="6.3: Profit over Capital Mexico"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_linetype_manual(values=c("solid"))+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  xlim(90, 200)

p23 <- scenario2 %>%
  ggplot( aes(x=period, y=((ucsla-mean(ucsla))/sd(ucsla)))) +
  geom_line( aes(linetype="6.4: Unit Costs Mexico"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_linetype_manual(values=c("solid"))+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  xlim(90, 200)



grid.arrange(p25, p22, p24, p23, nrow = 2)




pcomp <- scenario1 %>%
  ggplot( aes(x=period, y=((rerla-mean(rerla))/sd(rerla)))) +
  geom_line( aes(linetype="1.1: RER Mexico, Scenario 1"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_linetype_manual(values=c("solid"))+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  coord_cartesian(xlim=c(95, 130))+
  scale_x_continuous(breaks=seq(95, 130, 2))



pcomp1 <- scenario2 %>%
  ggplot( aes(x=period, y=((rerla-mean(rerla))/sd(rerla)))) +
  geom_line( aes(linetype="4.1: RER Mexico, Scenario 2"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_linetype_manual(values=c("solid"))+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  coord_cartesian(xlim=c(95, 130))+
  scale_x_continuous(breaks=seq(95, 130, 2))



grid.arrange(pcomp, pcomp1, nrow = 2)


##### Model RER validation ####

#Last update: 21/06/2021

#Clear all
#rm(list=ls(all=TRUE))

#Set directory
#setwd("C:/Users/Giuliano Toshiro/OneDrive - uniroma1.it/SFC/RER")

#save scenario 1 data in cvs
#write.csv(scenario1,"C:/Users/Giuliano Toshiro/OneDrive - uniroma1.it/SFC/RER/data_model_r.csv", row.names = TRUE)

#download data from WB
Data1<-wb_data(country = c("MX"), indicator = c("NY.GDP.MKTP.CN", "NE.CON.TOTL.CN", "NE.EXP.GNFS.CN","NE.GDI.TOTL.CN"), start_date = 2000, end_date = 2018, return_wide=TRUE, freq = "Y")
names(Data1)[names(Data1) == 'NE.CON.TOTL.CN'] <- 'cons'
names(Data1)[names(Data1) == 'NY.GDP.MKTP.CN'] <- 'gdp'
names(Data1)[names(Data1) == 'NE.EXP.GNFS.CN'] <- 'exp'
names(Data1)[names(Data1) == 'NE.GDI.TOTL.CN'] <- 'inv'
#write.csv(Datawb,"C:/Users/Giuliano Toshiro/OneDrive - uniroma1.it/SFC/RER/data_wb_mx_r.csv", row.names = TRUE)


#Upload data
#Data1<-read.csv("C:/Users/Giuliano Toshiro/OneDrive - uniroma1.it/SFC/RER/data_wb_mx_r.csv")
#Data2<-read.csv("C:/Users/Giuliano Toshiro/OneDrive - uniroma1.it/SFC/RER/data_model_r.csv",header=TRUE, sep=",")
#names(Data1) <- c('gdp','xrla')
#names(Data2) <- c('Y','C')

###########################################

#Model validation

#Install.packages("mFilter") #this command is necessary if mFilter has not been installed in your computer
library(mFilter)

#Create series
Y = as.matrix(Data1$gdp)
Ymod = as.matrix(scenario1$yla[95:113])
C = as.matrix(Data1$cons)
Cmod = as.matrix(scenario1$cla[95:113])
X = as.matrix(Data1$exp)
Xmod = as.matrix(scenario1$xla[95:113])
I = as.matrix(Data1$inv)
Imod = as.matrix(scenario1$ila_d[95:113])

#Create log series
Y_log<-log(Y)
Ymod_log<-log(Ymod)
C_log<-log(C)
Cmod_log<-log(Cmod)
X_log<-log(X)
Xmod_log<-log(Xmod)
I_log<-log(I)
Imod_log<-log(Imod)

### Output auto-correlation ###

#Filter series
Y.hp <- hpfilter((Y_log), freq=6.25, drift=TRUE)       #For actual data
Ymod.hp <- hpfilter((Ymod_log), freq=6.25, drift=TRUE)  #For synthetic data

#compute auto-correlation of cyclical components
acfY=acf(Y.hp$cycle,  plot=F)               #For actual data
acfYmod=acf(Ymod.hp$cycle, plot=F)          #For synthetic data

### Output-consumption cross-correlation ###

#Filter consumption series
C.hp <- hpfilter((C_log), freq=6.25, drift=TRUE)       #For actual data
Cmod.hp <- hpfilter((Cmod_log), freq=6.25, drift=TRUE)  #For synthetic data

#Compute cross-correlation of cyclical components
ccfC=ccf(Y.hp$cycle,C.hp$cycle, plot=F)               #For actual data
ccfCmod=ccf(Ymod.hp$cycle,Cmod.hp$cycle, plot=F)      #For synthetic data

### Output-export cross-correlation ###

#Filter export series
X.hp <- hpfilter((X_log), freq=6.25, drift=TRUE)       #For actual data
Xmod.hp <- hpfilter((Xmod_log), freq=6.25, drift=TRUE)  #For synthetic data

#Compute cross-correlation of cyclical components
ccfX=ccf(Y.hp$cycle,X.hp$cycle, plot=F)               #For actual data
ccfXmod=ccf(Ymod.hp$cycle,Xmod.hp$cycle, plot=F)      #For synthetic data

### Output-investment cross-correlation ###

#Filter export series
I.hp <- hpfilter((I_log), freq=6.25, drift=TRUE)       #For actual data
Imod.hp <- hpfilter((Imod_log), freq=6.25, drift=TRUE)  #For synthetic data

#Compute cross-correlation of cyclical components
ccfI=ccf(Y.hp$cycle,I.hp$cycle, plot=F)               #For actual data
ccfImod=ccf(Ymod.hp$cycle,Imod.hp$cycle, plot=F)      #For synthetic data


#Display, save and plot

png(file="model correlations new.png",width=844,height=622)
layout(mat = matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2))


#Plot and compare auto-correlations
plot(acfY$acf, ylab=" ", xlab="Lag", type="l", lty=1, col=1,lwd=1,ylim=c(-0.5,1),font.main=1,cex.main=0.75,main="A) Output auto-correlation",cex.axis=0.75,cex.lab=0.75)
lines(acfYmod$acf, type="l", lty=3, ylim=c(-0.5,1), col=2)
legend("topright", legend=c("Actual", "Simulated"), lty=c(1,3), col=c(1,2), bty="n", cex=0.75,box.lty=0)

#Plot and compare cross-correlations
plot(ccfC$acf, ylab=" ", xlab="Lag", type="l", lty=1, col=1,lwd=1,ylim=c(-0.5,1),font.main=1,cex.main=0.75,main="B) Output-consumption cross-correlation",cex.axis=0.75,cex.lab=0.75) #,xaxt="n")
lines(ccfCmod$acf, type="l", lty=3, ylim=c(-0.5,1), col=3)
legend("topleft", legend=c("Actual", "Simulated"), lty=c(1,3), col=c(1,3), bty="n", cex=0.75,box.lty=0)

#Plot and compare cross-correlations
plot(ccfX$acf, ylab=" ", xlab="Lag", type="l", lty=1, col=1,lwd=1,ylim=c(-0.5,1),font.main=1,cex.main=0.75,main="C) Output-export cross-correlation",cex.axis=0.75,cex.lab=0.75) #,xaxt="n")
lines(ccfXmod$acf, type="l", lty=3, ylim=c(-0.5,1), col=3)
legend("topleft", legend=c("Actual", "Simulated"), lty=c(1,3), col=c(1,3), bty="n", cex=0.75,box.lty=0)

#Plot and compare cross-correlations
plot(ccfI$acf, ylab=" ", xlab="Lag", type="l", lty=1, col=1,lwd=1,ylim=c(-0.5,1),font.main=1,cex.main=0.75,main="D) Output-Investment cross-correlation",cex.axis=0.75,cex.lab=0.75) #,xaxt="n")
lines(ccfImod$acf, type="l", lty=3, ylim=c(-0.5,1), col=3)
legend("topleft", legend=c("Actual", "Simulated"), lty=c(1,3), col=c(1,3), bty="n", cex=0.75,box.lty=0)


dev.off()

