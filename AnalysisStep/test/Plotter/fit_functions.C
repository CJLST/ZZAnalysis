Double_t myfunction(Double_t *x, Double_t *par)
{
   Float_t xx =x[0];
   Double_t f = par[0]*exp(-par[1]*xx)+par[2]*exp(-par[3]*xx);
   return f;
}

Double_t myfunctionErrUp(Double_t *x, Double_t *par)
{
   Float_t xx =x[0];
   Double_t f = par[0]*exp(-par[1]*xx)+par[2]*exp(-par[3]*xx);

   Double_t der[4]; 
   der[0] = exp(-par[1]*xx);
   der[1] = -par[0]*xx*exp(-par[1]*xx);   
   der[2] = exp(-par[3]*xx);
   der[3] = -par[2]*xx*exp(-par[3]*xx);

   Double_t fsigma = 0.;
   for (int j=0; j<4; j++) {
     for (int i=0; i<4; i++) {
       fsigma += der[i]*par[4+j+4*i]*der[j];
     }   
   }
   return f+sqrt(fabs(fsigma));
}

Double_t myfunctionErrDown(Double_t *x, Double_t *par)
{
   Float_t xx =x[0];
   Double_t f = par[0]*exp(-par[1]*xx)+par[2]*exp(-par[3]*xx);

   Double_t der[4]; 
   der[0] = exp(-par[1]*xx);
   der[1] = -par[0]*xx*exp(-par[1]*xx);   
   der[2] = exp(-par[3]*xx);
   der[3] = -par[2]*xx*exp(-par[3]*xx);

   Double_t fsigma = 0.;
   for (int j=0; j<4; j++) {
     for (int i=0; i<4; i++) {
       fsigma += der[i]*par[4+j+4*i]*der[j];
     }   
   }
   return f-sqrt(fabs(fsigma));
}

