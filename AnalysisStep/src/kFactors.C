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


float kfactor_qqZZ_qcd_M(float GENmassZZ, int finalState)
{

    // finalState=1 : 4e/4mu/4tau
    // finalState=2 : 2e2mu/2mutau/2e2tau

    float k=0.0;

    if (finalState==1) {
        k+=1.23613311013*(GENmassZZ>0.0&&GENmassZZ<=25.0);
        k+=1.17550314639*(GENmassZZ>25.0&&GENmassZZ<=50.0);
        k+=1.17044565911*(GENmassZZ>50.0&&GENmassZZ<=75.0);
        k+=1.03141209689*(GENmassZZ>75.0&&GENmassZZ<=100.0);
        k+=1.05285574912*(GENmassZZ>100.0&&GENmassZZ<=125.0);
        k+=1.11287217794*(GENmassZZ>125.0&&GENmassZZ<=150.0);
        k+=1.13361441158*(GENmassZZ>150.0&&GENmassZZ<=175.0);
        k+=1.10355603327*(GENmassZZ>175.0&&GENmassZZ<=200.0);
        k+=1.10053981637*(GENmassZZ>200.0&&GENmassZZ<=225.0);
        k+=1.10972676811*(GENmassZZ>225.0&&GENmassZZ<=250.0);
        k+=1.12069120525*(GENmassZZ>250.0&&GENmassZZ<=275.0);
        k+=1.11589101635*(GENmassZZ>275.0&&GENmassZZ<=300.0);
        k+=1.13906170314*(GENmassZZ>300.0&&GENmassZZ<=325.0);
        k+=1.14854594271*(GENmassZZ>325.0&&GENmassZZ<=350.0);
        k+=1.14616229031*(GENmassZZ>350.0&&GENmassZZ<=375.0);
        k+=1.14573157789*(GENmassZZ>375.0&&GENmassZZ<=400.0);
        k+=1.13829430515*(GENmassZZ>400.0&&GENmassZZ<=425.0);
        k+=1.15521193686*(GENmassZZ>425.0&&GENmassZZ<=450.0);
        k+=1.13679822698*(GENmassZZ>450.0&&GENmassZZ<=475.0);
        k+=1.13223956942*(GENmassZZ>475.0);
    }

    if (finalState==2) {
        k+=1.25094466582*(GENmassZZ>0.0&&GENmassZZ<=25.0);
        k+=1.22459455362*(GENmassZZ>25.0&&GENmassZZ<=50.0);
        k+=1.19287368979*(GENmassZZ>50.0&&GENmassZZ<=75.0);
        k+=1.04597506451*(GENmassZZ>75.0&&GENmassZZ<=100.0);
        k+=1.08323413771*(GENmassZZ>100.0&&GENmassZZ<=125.0);
        k+=1.09994968030*(GENmassZZ>125.0&&GENmassZZ<=150.0);
        k+=1.16698455800*(GENmassZZ>150.0&&GENmassZZ<=175.0);
        k+=1.10399053155*(GENmassZZ>175.0&&GENmassZZ<=200.0);
        k+=1.10592664340*(GENmassZZ>200.0&&GENmassZZ<=225.0);
        k+=1.10690381480*(GENmassZZ>225.0&&GENmassZZ<=250.0);
        k+=1.11194928918*(GENmassZZ>250.0&&GENmassZZ<=275.0);
        k+=1.13522586553*(GENmassZZ>275.0&&GENmassZZ<=300.0);
        k+=1.11895090244*(GENmassZZ>300.0&&GENmassZZ<=325.0);
        k+=1.13898508615*(GENmassZZ>325.0&&GENmassZZ<=350.0);
        k+=1.15463977506*(GENmassZZ>350.0&&GENmassZZ<=375.0);
        k+=1.17341664594*(GENmassZZ>375.0&&GENmassZZ<=400.0);
        k+=1.20093349763*(GENmassZZ>400.0&&GENmassZZ<=425.0);
        k+=1.18915554919*(GENmassZZ>425.0&&GENmassZZ<=450.0);
        k+=1.18546007375*(GENmassZZ>450.0&&GENmassZZ<=475.0);
        k+=1.12864505708*(GENmassZZ>475.0);
    }

    if (k==0.0) return 1.1;
    else return k; // if something goes wrong return inclusive k-factor

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
