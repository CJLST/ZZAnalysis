 float kfactor_qqZZ_qcd_dPhi(float GENabsdPhiZZ, int finalState)
{

    // finalState=1 : 4e/4mu/4tau
    // finalState=2 : 2e2mu/2mutau/2e2tau

    float k=0.0;

    if (finalState==1) {        
        k+=1.515838921760*(GENabsdPhiZZ>0.0&&GENabsdPhiZZ<=0.1);
        k+=1.496256665410*(GENabsdPhiZZ>0.1&&GENabsdPhiZZ<=0.2);
        k+=1.495522061910*(GENabsdPhiZZ>0.2&&GENabsdPhiZZ<=0.3);
        k+=1.483273154250*(GENabsdPhiZZ>0.3&&GENabsdPhiZZ<=0.4);
        k+=1.465589701130*(GENabsdPhiZZ>0.4&&GENabsdPhiZZ<=0.5);
        k+=1.491500887510*(GENabsdPhiZZ>0.5&&GENabsdPhiZZ<=0.6);
        k+=1.441183580450*(GENabsdPhiZZ>0.6&&GENabsdPhiZZ<=0.7);
        k+=1.440830603990*(GENabsdPhiZZ>0.7&&GENabsdPhiZZ<=0.8);
        k+=1.414339019120*(GENabsdPhiZZ>0.8&&GENabsdPhiZZ<=0.9);
        k+=1.422534218560*(GENabsdPhiZZ>0.9&&GENabsdPhiZZ<=1.0);
        k+=1.401037066000*(GENabsdPhiZZ>1.0&&GENabsdPhiZZ<=1.1);
        k+=1.408539428810*(GENabsdPhiZZ>1.1&&GENabsdPhiZZ<=1.2);
        k+=1.381247744080*(GENabsdPhiZZ>1.2&&GENabsdPhiZZ<=1.3);
        k+=1.370553357430*(GENabsdPhiZZ>1.3&&GENabsdPhiZZ<=1.4);
        k+=1.347323316000*(GENabsdPhiZZ>1.4&&GENabsdPhiZZ<=1.5);
        k+=1.340113437450*(GENabsdPhiZZ>1.5&&GENabsdPhiZZ<=1.6);
        k+=1.312661036510*(GENabsdPhiZZ>1.6&&GENabsdPhiZZ<=1.7);
        k+=1.290055062010*(GENabsdPhiZZ>1.7&&GENabsdPhiZZ<=1.8);
        k+=1.255322614790*(GENabsdPhiZZ>1.8&&GENabsdPhiZZ<=1.9);
        k+=1.254455642450*(GENabsdPhiZZ>1.9&&GENabsdPhiZZ<=2.0);
        k+=1.224047664420*(GENabsdPhiZZ>2.0&&GENabsdPhiZZ<=2.1);
        k+=1.178816782670*(GENabsdPhiZZ>2.1&&GENabsdPhiZZ<=2.2);
        k+=1.162624827140*(GENabsdPhiZZ>2.2&&GENabsdPhiZZ<=2.3);
        k+=1.105401140940*(GENabsdPhiZZ>2.3&&GENabsdPhiZZ<=2.4);
        k+=1.074749265690*(GENabsdPhiZZ>2.4&&GENabsdPhiZZ<=2.5);
        k+=1.021864599380*(GENabsdPhiZZ>2.5&&GENabsdPhiZZ<=2.6);
        k+=0.946334793286*(GENabsdPhiZZ>2.6&&GENabsdPhiZZ<=2.7);
        k+=0.857458082628*(GENabsdPhiZZ>2.7&&GENabsdPhiZZ<=2.8);
        k+=0.716607670482*(GENabsdPhiZZ>2.8&&GENabsdPhiZZ<=2.9);
        k+=1.132841784840*(GENabsdPhiZZ>2.9&&GENabsdPhiZZ<=3.1416);
    }

    if (finalState==2) {
       k+=1.513834489150*(GENabsdPhiZZ>0.0&&GENabsdPhiZZ<=0.1);
       k+=1.541738780180*(GENabsdPhiZZ>0.1&&GENabsdPhiZZ<=0.2);
       k+=1.497829632510*(GENabsdPhiZZ>0.2&&GENabsdPhiZZ<=0.3);
       k+=1.534956782920*(GENabsdPhiZZ>0.3&&GENabsdPhiZZ<=0.4);
       k+=1.478217033060*(GENabsdPhiZZ>0.4&&GENabsdPhiZZ<=0.5);
       k+=1.504330859290*(GENabsdPhiZZ>0.5&&GENabsdPhiZZ<=0.6);
       k+=1.520626246850*(GENabsdPhiZZ>0.6&&GENabsdPhiZZ<=0.7);
       k+=1.507013090030*(GENabsdPhiZZ>0.7&&GENabsdPhiZZ<=0.8);
       k+=1.494243156250*(GENabsdPhiZZ>0.8&&GENabsdPhiZZ<=0.9);
       k+=1.450536096150*(GENabsdPhiZZ>0.9&&GENabsdPhiZZ<=1.0);
       k+=1.460812521660*(GENabsdPhiZZ>1.0&&GENabsdPhiZZ<=1.1);
       k+=1.471603622200*(GENabsdPhiZZ>1.1&&GENabsdPhiZZ<=1.2);
       k+=1.467700038200*(GENabsdPhiZZ>1.2&&GENabsdPhiZZ<=1.3);
       k+=1.422408690640*(GENabsdPhiZZ>1.3&&GENabsdPhiZZ<=1.4);
       k+=1.397184022730*(GENabsdPhiZZ>1.4&&GENabsdPhiZZ<=1.5);
       k+=1.375593447520*(GENabsdPhiZZ>1.5&&GENabsdPhiZZ<=1.6);
       k+=1.391901318370*(GENabsdPhiZZ>1.6&&GENabsdPhiZZ<=1.7);
       k+=1.368564350560*(GENabsdPhiZZ>1.7&&GENabsdPhiZZ<=1.8);
       k+=1.317884804290*(GENabsdPhiZZ>1.8&&GENabsdPhiZZ<=1.9);
       k+=1.314019950800*(GENabsdPhiZZ>1.9&&GENabsdPhiZZ<=2.0);
       k+=1.274641749910*(GENabsdPhiZZ>2.0&&GENabsdPhiZZ<=2.1);
       k+=1.242346606820*(GENabsdPhiZZ>2.1&&GENabsdPhiZZ<=2.2);
       k+=1.244727403840*(GENabsdPhiZZ>2.2&&GENabsdPhiZZ<=2.3);
       k+=1.146259351670*(GENabsdPhiZZ>2.3&&GENabsdPhiZZ<=2.4);
       k+=1.107804993520*(GENabsdPhiZZ>2.4&&GENabsdPhiZZ<=2.5);
       k+=1.042053646740*(GENabsdPhiZZ>2.5&&GENabsdPhiZZ<=2.6);
       k+=0.973608545141*(GENabsdPhiZZ>2.6&&GENabsdPhiZZ<=2.7);
       k+=0.872169942668*(GENabsdPhiZZ>2.7&&GENabsdPhiZZ<=2.8);
       k+=0.734505279177*(GENabsdPhiZZ>2.8&&GENabsdPhiZZ<=2.9);
       k+=1.163152837230*(GENabsdPhiZZ>2.9&&GENabsdPhiZZ<=3.1416);       
    }
    if (k==0.0) return 1.1; // if something goes wrong return inclusive k-factor
    else return k;

}

float xsec_qqZZ_qcd_M(float GenMassZZ, int finalState, int order){ // Order: 0=LO, 1=NLO, 2=NNLO
  const int nbins=20;
  static const float xsec[2][nbins][4]={
    {
      { 0, 1.095828532, 1.623089848, 2.006355102 },
      { 25, 1.70498411, 2.55680607, 3.00553358 },
      { 50, 0.7214786, 1.07670322, 1.26022261 },
      { 75, 9.93585169, 12.49651341, 12.8890551 },
      { 100, 1.88979055, 2.32301901, 2.44580392 },
      { 125, 0.80689284, 1.06246962, 1.18239288 },
      { 150, 0.420080936, 0.576994967, 0.65408981 },
      { 175, 0.573751417, 0.755001119, 0.83318604 },
      { 200, 0.6426366, 0.84359613, 0.92841113 },
      { 225, 0.460870086, 0.6155138, 0.68305214 },
      { 250, 0.325131608, 0.440635023, 0.493815795 },
      { 275, 0.234474317, 0.321194629, 0.358418201 },
      { 300, 0.173099357, 0.239157872, 0.272415573 },
      { 325, 0.130327192, 0.18182296, 0.208832023 },
      { 350, 0.100194403, 0.140405709, 0.160927729 },
      { 375, 0.078065965, 0.109969337, 0.125995342 },
      { 400, 0.061860538, 0.087368406, 0.099450959 },
      { 425, 0.049470606, 0.070200812, 0.081096816 },
      { 450, 0.040036726, 0.057035286, 0.064837612 },
      { 475, 0.03272652, 0.046808668, 0.053386015 }
    },
    {
      { 0, 5.80657303, 8.38413944, 10.48809451 },
      { 25, 4.62863121, 6.84873304, 8.38692118 },
      { 50, 1.61514729, 2.39428421, 2.85607864 },
      { 75, 19.10610496, 24.00530313, 25.10894849 },
      { 100, 3.86816714, 4.7520914, 5.14762763 },
      { 125, 1.64560534, 2.16841494, 2.38514732 },
      { 150, 0.85582638, 1.17224019, 1.3679862 },
      { 175, 1.1629688, 1.52672919, 1.68549457 },
      { 200, 1.29516749, 1.69836336, 1.87826529 },
      { 225, 0.92621932, 1.23581661, 1.36793012 },
      { 250, 0.65489237, 0.88543635, 0.98456032 },
      { 275, 0.470967807, 0.64402568, 0.73111461 },
      { 300, 0.347254403, 0.480862926, 0.538062005 },
      { 325, 0.261981026, 0.363080611, 0.413543401 },
      { 350, 0.201165289, 0.280556435, 0.323941619 },
      { 375, 0.157331951, 0.220998012, 0.259322746 },
      { 400, 0.124083729, 0.174598837, 0.209681592 },
      { 425, 0.099317367, 0.139860047, 0.166315351 },
      { 450, 0.080436156, 0.114091349, 0.135250739 },
      { 475, 0.065390164, 0.093794404, 0.104784251 }
    }
  };
  int cbin=-1;
  for (int ix=0; ix<nbins-1; ix++){
    if (GenMassZZ>=xsec[finalState-1][ix][0] && GenMassZZ<xsec[finalState-1][ix+1][0]){ cbin=ix; break; }
  }
  if (cbin<0) cbin=nbins-1;
  return xsec[finalState-1][cbin][order+1];
}
float kfactor_qqZZ_qcd_M(float GENmassZZ, int finalState, int order) // Order: 1=NLO, 2=NNLO
{
  // finalState=1 : 4e/4mu/4tau
  // finalState=2 : 2e2mu/2mutau/2e2tau

  float k=1;
  float xsec_lo = xsec_qqZZ_qcd_M(GENmassZZ, finalState, 0);
  float xsec = xsec_qqZZ_qcd_M(GENmassZZ, finalState, order);
  k = xsec/xsec_lo;
  return k;
}


float kfactor_qqZZ_qcd_Pt(float GENpTZZ, int finalState)
{

    // finalState=1 : 4e/4mu/4tau
    // finalState=2 : 2e2mu/2mutau/2e2tau

    float k=0.0;

    if (finalState==1) {
        k+=0.64155491983*(GENpTZZ>0.0&&GENpTZZ<=5.0);
        k+=1.09985240531*(GENpTZZ>5.0&&GENpTZZ<=10.0);
        k+=1.29390628654*(GENpTZZ>10.0&&GENpTZZ<=15.0);
        k+=1.37859998571*(GENpTZZ>15.0&&GENpTZZ<=20.0);
        k+=1.42430263312*(GENpTZZ>20.0&&GENpTZZ<=25.0);
        k+=1.45038493266*(GENpTZZ>25.0&&GENpTZZ<=30.0);
        k+=1.47015377651*(GENpTZZ>30.0&&GENpTZZ<=35.0);
        k+=1.48828685748*(GENpTZZ>35.0&&GENpTZZ<=40.0);
        k+=1.50573440448*(GENpTZZ>40.0&&GENpTZZ<=45.0);
        k+=1.50211655928*(GENpTZZ>45.0&&GENpTZZ<=50.0);
        k+=1.50918720827*(GENpTZZ>50.0&&GENpTZZ<=55.0);
        k+=1.52463089491*(GENpTZZ>55.0&&GENpTZZ<=60.0);
        k+=1.52400838378*(GENpTZZ>60.0&&GENpTZZ<=65.0);
        k+=1.52418067701*(GENpTZZ>65.0&&GENpTZZ<=70.0);
        k+=1.55424382578*(GENpTZZ>70.0&&GENpTZZ<=75.0);
        k+=1.52544284222*(GENpTZZ>75.0&&GENpTZZ<=80.0);
        k+=1.57896384602*(GENpTZZ>80.0&&GENpTZZ<=85.0);
        k+=1.53034682567*(GENpTZZ>85.0&&GENpTZZ<=90.0);
        k+=1.56147329708*(GENpTZZ>90.0&&GENpTZZ<=95.0);
        k+=1.54468169268*(GENpTZZ>95.0&&GENpTZZ<=100.0);
        k+=1.57222952415*(GENpTZZ>100.0);
    }

    if (finalState==2) {
        k+=0.743602533303*(GENpTZZ>0.0&&GENpTZZ<=5.0);
        k+=1.14789453219*(GENpTZZ>5.0&&GENpTZZ<=10.0);
        k+=1.33815867892*(GENpTZZ>10.0&&GENpTZZ<=15.0);
        k+=1.41420044104*(GENpTZZ>15.0&&GENpTZZ<=20.0);
        k+=1.45511318916*(GENpTZZ>20.0&&GENpTZZ<=25.0);
        k+=1.47569225244*(GENpTZZ>25.0&&GENpTZZ<=30.0);
        k+=1.49053003693*(GENpTZZ>30.0&&GENpTZZ<=35.0);
        k+=1.50622827695*(GENpTZZ>35.0&&GENpTZZ<=40.0);
        k+=1.50328889799*(GENpTZZ>40.0&&GENpTZZ<=45.0);
        k+=1.52186945281*(GENpTZZ>45.0&&GENpTZZ<=50.0);
        k+=1.52043468754*(GENpTZZ>50.0&&GENpTZZ<=55.0);
        k+=1.53977869986*(GENpTZZ>55.0&&GENpTZZ<=60.0);
        k+=1.53491994434*(GENpTZZ>60.0&&GENpTZZ<=65.0);
        k+=1.51772882172*(GENpTZZ>65.0&&GENpTZZ<=70.0);
        k+=1.54494489131*(GENpTZZ>70.0&&GENpTZZ<=75.0);
        k+=1.57762411697*(GENpTZZ>75.0&&GENpTZZ<=80.0);
        k+=1.55078339014*(GENpTZZ>80.0&&GENpTZZ<=85.0);
        k+=1.57078191891*(GENpTZZ>85.0&&GENpTZZ<=90.0);
        k+=1.56162666568*(GENpTZZ>90.0&&GENpTZZ<=95.0);
        k+=1.54183774627*(GENpTZZ>95.0&&GENpTZZ<=100.0);
        k+=1.58485762205*(GENpTZZ>100.0);
    }

    if (k==0.0) return 1.1;
    else return k; // if something goes wrong return inclusive k-factor

}
