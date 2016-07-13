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

float xsec_qqZZ_qcd_M(float GenMassZZ, int finalstate, int order){ // Order: 0=LO, 1=NLO, 2=NNLO
  const int nbins=33;
  float xsec[2][nbins][4]={
    {
      { 0, 0.151735412, 0.222726408, 0.293690432 },
      { 15, 1.41980404, 2.11137938, 2.5264739 },
      { 30, 1.00724982, 1.51166097, 1.78968853 },
      { 45, 0.56014662, 0.84401206, 1.00089457 },
      { 60, 0.38335535, 0.56682032, 0.66136386 },
      { 75, 2.42358653, 3.10410551, 3.2618054 },
      { 90, 8.04426827, 10.03302569, 10.29723193 },
      { 105, 1.10387829, 1.35922979, 1.43125373 },
      { 120, 0.6507191, 0.83842407, 0.91275257 },
      { 135, 0.41008289, 0.54721698, 0.61420827 },
      { 150, 0.277400115, 0.37905505, 0.42935168 },
      { 165, 0.21094836, 0.291687546, 0.32936606 },
      { 180, 0.359365418, 0.47142866, 0.51922164 },
      { 195, 0.42556087, 0.55455072, 0.61461623 },
      { 210, 0.36319419, 0.47887024, 0.52313137 },
      { 225, 0.295819962, 0.39414625, 0.43627003 },
      { 240, 0.23916677, 0.32143823, 0.36012781 },
      { 255, 0.194457125, 0.263397475, 0.293796426 },
      { 270, 0.159411523, 0.217591203, 0.244570248 },
      { 285, 0.131620631, 0.180770294, 0.200521622 },
      { 300, 0.109827764, 0.151443805, 0.171604481 },
      { 315, 0.092287447, 0.128034787, 0.147080525 },
      { 330, 0.077948679, 0.108895098, 0.125043034 },
      { 345, 0.066529942, 0.093070947, 0.107159505 },
      { 360, 0.05702712, 0.079941904, 0.09128778 },
      { 375, 0.04907944, 0.069058479, 0.07896666 },
      { 390, 0.042513623, 0.060000127, 0.069100158 },
      { 405, 0.037061048, 0.052371461, 0.059399323 },
      { 420, 0.032355303, 0.045787495, 0.052979827 },
      { 435, 0.028387695, 0.040320993, 0.046097149 },
      { 450, 0.024981894, 0.035524204, 0.040728959 },
      { 465, 0.022130763, 0.031603791, 0.035695107 },
      { 480, 0.025650589, 0.036715959, 0.041799561 }
    },
    {
      { 0, 1.77877173, 2.51973134, 3.27103491 },
      { 15, 5.4789198, 7.9974481, 9.851197 },
      { 30, 2.64739638, 3.9225219, 4.798018 },
      { 45, 1.31051581, 1.95878833, 2.35818671 },
      { 60, 0.83474781, 1.22866702, 1.45265771 },
      { 75, 4.68433656, 5.99954263, 6.35212199 },
      { 90, 15.5067952, 19.3086019, 20.0753438 },
      { 105, 2.26265939, 2.78796761, 3.13020973 },
      { 120, 1.33195548, 1.71376924, 1.81602552 },
      { 135, 0.83413081, 1.11592809, 1.2680224 },
      { 150, 0.56407802, 0.77044142, 0.89728304 },
      { 165, 0.43082914, 0.59329019, 0.67985981 },
      { 180, 0.72880377, 0.95192561, 1.05598974 },
      { 195, 0.85834031, 1.11680623, 1.23640319 },
      { 210, 0.73191143, 0.96486929, 1.06221028 },
      { 225, 0.5948062, 0.79058588, 0.89045297 },
      { 240, 0.48094244, 0.64566121, 0.69954409 },
      { 255, 0.39101624, 0.52942727, 0.59281466 },
      { 270, 0.32112943, 0.43727125, 0.48447369 },
      { 285, 0.264185187, 0.36233303, 0.41631964 },
      { 300, 0.220506788, 0.305050888, 0.33994799 },
      { 315, 0.184834015, 0.256532412, 0.286630441 },
      { 330, 0.157019617, 0.217005367, 0.24685187 },
      { 345, 0.133450361, 0.186079557, 0.215443611 },
      { 360, 0.114589937, 0.159831748, 0.186673113 },
      { 375, 0.098838752, 0.138617897, 0.161329676 },
      { 390, 0.085695022, 0.120077664, 0.144923436 },
      { 405, 0.074362057, 0.10500367, 0.123982312 },
      { 420, 0.06468529, 0.090999441, 0.108940818 },
      { 435, 0.057151926, 0.080758224, 0.096143447 },
      { 450, 0.050213208, 0.071140775, 0.086065149 },
      { 465, 0.044433347, 0.063193852, 0.071480181 },
      { 480, 0.051179765, 0.073551126, 0.08248966 }
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
