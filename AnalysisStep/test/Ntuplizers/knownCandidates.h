#include <string>
#include <map>

std::string findLabel(int run, int ls, long event) {
  static std::map<int, std::map<int, std::map<long, std::string> > > knownEventLabels;  

  if (knownEventLabels.size()==0) {
    knownEventLabels[146511][724][504867308]   = "A  "; 
    knownEventLabels[147926][347][368148849]   = "B  "; 
    knownEventLabels[163334][499][286336207]   = "C  "; 
    knownEventLabels[163659][463][344708580]   = "D  "; 
    knownEventLabels[163795][34][30998576]     = "E  "; 
    knownEventLabels[163817][174][155679852]   = "F  "; 
    knownEventLabels[165633][303][394010457]   = "G  "; 
    knownEventLabels[165970][236][275108397]   = "R-A";
    knownEventLabels[166408][724][917379387]   = "H  "; 
    knownEventLabels[166438][79][78213037]     = "I  "; 
    knownEventLabels[166438][768][862270386]   = "R-B";
    knownEventLabels[166512][281][337493970]   = "J  "; 
    knownEventLabels[166950][1373][1491724484] = "K  ";
    knownEventLabels[167281][386][480301165]   = "L  "; 
    knownEventLabels[167282][44][44166176]     = "R-C";
    knownEventLabels[167284][1213][1038911933] = "M  "; 
    knownEventLabels[167675][829][876658967]   = "N  "; 
    knownEventLabels[167807][750][966824024]   = "O  "; 
    knownEventLabels[171106][127][141954801]   = "P  "; 
    knownEventLabels[171369][150][160966858]   = "Q  "; 
    knownEventLabels[172163][128][191231387]   = "R  "; 
    knownEventLabels[172208][75][66033190]     = "S  "; 
    knownEventLabels[172620][242][218903169]   = "R-D";
    knownEventLabels[172799][11][10347106]     = "T  "; 
    knownEventLabels[172802][125][107360878]   = "U  "; 
    knownEventLabels[172822][2004][2554393033] = "V  "; 
    knownEventLabels[172868][689][933807102]   = "W  "; 
    knownEventLabels[172949][840][1188043146]  = "X  "; 
    knownEventLabels[172952][466][559839432]   = "Y  "; 
    knownEventLabels[173243][12][16706390]     = "Z  "; 
    knownEventLabels[173657][85][65557571]     = "AA "; 
    knownEventLabels[173659][270][389185367]   = "AB "; 
    knownEventLabels[173692][2066][2722114329] = "AC "; 
    knownEventLabels[175906][190][227517585]   = "AD "; 
    knownEventLabels[175921][220][297753357]   = "AE "; 
    knownEventLabels[175921][349][495614354]   = "AF "; 
    knownEventLabels[175974][9][7526662]       = "AG "; 
    knownEventLabels[176207][206][256888239]   = "AH "; 
    knownEventLabels[176201][182][261184429]   = "R-E";
    knownEventLabels[176304][299][417897294]   = "AI ";   
    knownEventLabels[176304][300][418052877]   = "AJ "; 
    knownEventLabels[176309][950][1340034258]  = "R-F";
    knownEventLabels[176309][224][257489763]   = "AK ";
    knownEventLabels[176468][128][215855118]   = "AL "; 
    knownEventLabels[176548][231][403771114]   = "AM "; 
    knownEventLabels[176799][24][35688265]     = "AN "; 
    knownEventLabels[176886][631][1057019814]  = "AO "; 
    knownEventLabels[177074][384][588602439]   = "AP "; 
    knownEventLabels[177139][183][290826062]   = "AQ "; 
    knownEventLabels[177222][227][339499459]   = "AR "; 
    knownEventLabels[177730][1248][1995624376] = "R-G";
    knownEventLabels[177782][99][72158025]     = "AS "; 
    knownEventLabels[177790][168][222240677]   = "AT "; 
    knownEventLabels[177790][527][657843813]   = "AU "; 
    knownEventLabels[177875][133][148667118]   = "AV "; 
    knownEventLabels[178100][236][326364918]   = "AW "; 
    knownEventLabels[178116][437][709511403]   = "AX "; 
    knownEventLabels[178421][86][87514902]     = "AY "; 
    knownEventLabels[178421][973][1450980155]  = "AZ "; 
    knownEventLabels[178479][210][298608854]   = "BA "; 
    knownEventLabels[178479][369][589085976]   = "BB ";
    knownEventLabels[178703][137][191352626]   = "BC ";
    knownEventLabels[178731][192][248562036]   = "BD ";
    knownEventLabels[178866][82][140063742]    = "BE ";
    knownEventLabels[178970][67][57399691]     = "BF ";
    knownEventLabels[178970][103][122998167]   = "BG ";
    knownEventLabels[179434][52][86225612]     = "BH ";
    knownEventLabels[179452][1056][1459855927] = "BI ";
    knownEventLabels[179476][30][30532070]     = "BJ ";
    knownEventLabels[179563][177][287281642]   = "BK ";
    knownEventLabels[179563][871][1409064222]  = "BL ";
    knownEventLabels[180076][46][79350642]     = "BM ";
    knownEventLabels[180250][326][591651181]   = "BN ";
    //
//     knownEventLabels[175834][87][67783615]     = "NC ";
//     knownEventLabels[177318][169][270676815]   = "NC ";
//     knownEventLabels[177515][163][286020249]   = "NC ";
//     knownEventLabels[177515][620][983847906]   = "NC ";
  }

  std::string result =  knownEventLabels[run][ls][event];
  
  if (result == "") return "???";
  return result;
}
