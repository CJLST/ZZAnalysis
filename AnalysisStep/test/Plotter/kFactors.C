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

