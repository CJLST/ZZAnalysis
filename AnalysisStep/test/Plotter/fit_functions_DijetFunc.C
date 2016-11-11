Double_t myfunction(Double_t *x, Double_t *par)
{
   Float_t xx =x[0];
   Double_t f = par[0]*exp(-xx/(par[1]+par[2]*xx));
   return f;
}

Double_t myfunctionErrUp(Double_t *x, Double_t *par)
{
   Float_t xx =x[0];
   Double_t f = par[0]*exp(-xx/(par[1]+par[2]*xx));

   Double_t der[3]; 
   der[0] = exp(-xx/(par[1]+par[2]*xx));
   der[1] = par[0]*exp(-xx/(par[1]+par[2]*xx))/pow(par[1]+par[2]*xx,2);   
   der[2] = par[0]*xx*exp(-xx/(par[1]+par[2]*xx))/pow(par[1]+par[2]*xx,2);
 
   Double_t fsigma = 0.;
   for (int j=0; j<3; j++) {
     for (int i=0; i<3; i++) {
       fsigma += der[i]*par[3+j+3*i]*der[j];
     }   
   }
   return f+sqrt(fabs(fsigma));
}

Double_t myfunctionErrDown(Double_t *x, Double_t *par)
{
   Float_t xx =x[0];
   Double_t f = par[0]*exp(-xx/(par[1]+par[2]*xx));

   Double_t der[3]; 
   der[0] = exp(-xx/(par[1]+par[2]*xx));
   der[1] = par[0]*exp(-xx/(par[1]+par[2]*xx))/pow(par[1]+par[2]*xx,2);   
   der[2] = par[0]*xx*exp(-xx/(par[1]+par[2]*xx))/pow(par[1]+par[2]*xx,2);
 
   Double_t fsigma = 0.;
   for (int j=0; j<3; j++) {
     for (int i=0; i<3; i++) {
       fsigma += der[i]*par[3+j+3*i]*der[j];
     }   
   }
   return f-sqrt(fabs(fsigma));
}

