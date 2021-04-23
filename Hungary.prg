cls; format /ld 6,6;
ser=2;

load yyy[69,3] = C:\GAUSS\Hansen\RBCIShif\Pol_lev.dat ;
y=yyy[2:69,ser];
@
load yyy[63,3] = C:\GAUSS\Hansen\RBCIShif\Hungary.dat ;
y=yyy[2:61,ser];
@
sername = "" $+yyy[1,ser];
"Series: " sername;
print "";

call main(y,3,4,7);  @ 0 - break in c only (series have no trend)
                       1 - break in c only (series include trend)
                       2 - break in trend only 
                       3 - break in c and trend @
end;


/*************************************************************************
----PROC MAIN
----FORMAT: call  main(y,model,choice,k)
----INPUT:      y - depend variable
                model - choice for model        
						=1  A
                        =2  B
                        =3  C
                choice - only in ADF test,  
						=1  pre-specified AR lag
                        =2  AIC-chosen AR lag
                        =3  BIC-chosen AR lag
                        =4  downward-t-chosen AR lag
                k - maximum lag for ADF test
----GLOBAL VARIABLES: none
----EXTERNAL PROCEDURES: adf,  phillips
----NB: Constant included in regression
************************************************************************/

/*
****************  Main procedure *******************
*/

proc(0)=main(y,model,choice,k);
   local t,n,final,begin,tstat,x,lag,j,DU, DT,temp1,temp2,temp3,temp4;
   local breakpt1,breakpt2,breakpta,za,zt;
   n=rows(y);
   begin=round(0.15*n);
   final=round(0.85*n);
   temp1=zeros(final-begin+1,1);
   temp2=temp1;
   temp3=temp1;
   temp4=temp1;
   t=begin; 
   do while t<=final;
     DU=zeros(t,1)|ones(n-t,1);
     DT=zeros(t,1)|seqa(1,1,n-t);
     @ adjust regressors for different models @
     if model==0;
        x=ones(n,1)~DU;
     elseif model==1;
        x=ones(n,1)~DU~seqa(1,1,n);
     elseif model==2;
        x=ones(n,1)~seqa(1,1,n)~DT;
     elseif model==3;
        x=ones(n,1)~DU~seqa(1,1,n)~DT;
     endif;

     @ compute ADF for each t  @
     {temp1[t-begin+1],temp2[t-begin+1]}=adf(y,x,k,choice);

   
   t=t+1;endo; 
   
 
   @  ADF test @
   tstat=minc(temp1); 
   lag=temp2[minindc(temp1)];
   breakpta=(minindc(temp1)+begin+1); 
   print "******** Break ADF Test *****************************";
   print "t-statistic           = " tstat;
   print "AR lag                = " lag;
   print "breakpoint            = " breakpta;
   print " ";

   if model==0;
     x = ones(n,1);
     {temp1,temp2}=adf(y,x,k,choice);
     print "******** Standard ADF Test ***************************";
     print "t-statistic           = " temp1;
     print "AR lag                = " temp2;
     print " ";
     print "Tests include: constant only";
     print "";
   else;
     x = ones(n,1)~seqa(1,1,n);
     {temp1,temp2}=adf(y,x,k,choice);
     print "******** Standard ADF Test ***************************";
     print "t-statistic           = " temp1;
     print "AR lag                = " temp2;
     print "";
     print "Tests include: constant and trend";
     print "";
   endif;
  
  
retp;
endp;
@ -------------------------------------------------------------- @


/**********************  PROC ADF  *****************************
**   FORMAT
**          { stat,lag } = adf(y,x)
**   INPUT
**        y - dependent variable
**        x - independent variables
**   OUTPUT
**  stata - ADF statistic
**  lag - the lag length
**   GLOBAL VARIABLES: none
**   EXTERNAL PROCEDURES: estimate
**********************************************************************/

/*
*************** ADF for each breakpoint ********************
*/
proc(2) = adf(y,x,kmax,choice);
   local b,m,e,e1,n,n1,sig2,se,xdy,ydy,j,tstat,dy,temp1,temp2;
   local lag,k,ic,aic,bic;
   @ compute ADF  @
   n=rows(y);
   dy=y[2:n]-y[1:n-1]; @ difference of dependent @
   ic=0;
   k=kmax;
   temp1=zeros(kmax+1,1);
   temp2=zeros(kmax+1,1);
   do while k>=0;
      ydy=dy[k+1:n-1];
      n1=rows(ydy);
      @  set up matrix for lagged dependent + trend dummies @
      xdy=y[k+1:n-1]~x[k+1:n-1,1:cols(x)];
      j=1;
      do while j <= k;
         xdy=xdy~dy[k+1-j:n-1-j];
         j=j+1;
      endo;
      {b,e1,sig2,se}=estimate(ydy,xdy);
      if choice==1;  @ K is pre-specified @
          temp1[k+1]=-1000;   @ set an random negative constant @
          temp2[k+1]=b[1]/se[1];
          break;
      elseif choice==2;  @ K is determined by AIC @
         aic=ln(e1'e1/n1)+2*cols(xdy)/n1;
         ic=aic;
      elseif choice==3;  @ K is determined by BIC @
         bic=ln(e1'e1/n1)+cols(xdy)*ln(n1)/n1;
         ic=bic;
      elseif choice==4; @K is determined by downward t @
         if abs(b[cols(xdy)]/se[cols(xdy)]) >= 1.96 or k==0;
            temp1[k+1]=-1000;    @ set an random negative constant @
            temp2[k+1]=b[1]/se[1];
            break;
         endif;
      endif;
      temp1[k+1]=ic;
      temp2[k+1]=b[1]/se[1];
      k=k-1;
   endo;
   lag=minindc(temp1);
   tstat=temp2[lag];
   retp(tstat,lag-1);
endp;
@ ------------------------------------------------------------ @

/**********************  PROC ESTIMATE  *****************************
**   FORMAT
**          { b,e,sig2,se } = estimate(y,x)
**   INPUT
**        y  - dependent variable
**        x - independent variables
**   OUTPUT
**  b - OLS estimates
**  e - residuals
**  sig2 - variance
**  se - standard error for coefficients
**   GLOBAL VARIABLES: none
**********************************************************************/
/* *****  ols regression ****** */
proc(4) = estimate(y,x);
   local m, b, e, sig2, se;
   m=invpd(moment(x,0));
   b=m*(x'y);
   e=y-x*b;
   sig2=(e'e)/(rows(y)-cols(x));
   se=sqrt(diag(m)*sig2);
   retp(b,e,sig2,se);
endp;
@ ---------------------------------------------------------------- @





/**********************  PROC ESTIMATE  *****************************
**   FORMAT
**          { b,e,sig2,se } = estimate(y,x)
**   INPUT
**        y  - dependent variable
**        x - independent variables
**   OUTPUT
**  b - OLS estimates
**  e - residuals
**  sig2 - variance
**  se - standard error for coefficients
**   GLOBAL VARIABLES: none
**********************************************************************/
/* *****  ols regression ****** */
proc(4) = estimate(y,x);
   local m, b, e, sig2, se;
   m=invpd(moment(x,0));
   b=m*(x'y);
   e=y-x*b;
   sig2=(e'e)/(rows(y)-cols(x));
   se=sqrt(diag(m)*sig2);
   retp(b,e,sig2,se);
endp;
@ ---------------------------------------------------------------- @

