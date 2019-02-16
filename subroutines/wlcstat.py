#!/usr/bin/python

import numpy as np
import os

# Calculate the order parameter for a nematic liquid crystal
# Continued fraction formulation of the partition function

    
def gwlc(r,n,d,k0,kf,kstep,ordmax,ord,reslayer):

# Calculate the end-distribution function for the wormlike chain model
# Perform numerical inverse Fourier transform


# Tolerance for adaptive reduction of eigenvalues

    tol=1e-22

# Maximum k increment (frequency), Kstep is a minimum to enforce manually
    
    kperiod=min(kstep,2*pi/n/10)
    iend=len(R)

    grn=np.zeros(len(R))

    counter=0
    k=k0+kperiod
    while k <= kf:
        counter=counter+1
        
#        [GK,SepRes]=gkwlc(N,K,d,ORDmax,ORD,ResLayer);
        
        
#    % if residual is less the a tolerance, reduce number of eigenvalues for the
#    % next steps
#        if max(K'.*abs(SepRes(:,ORDmax)))<Tol && K>2840
#            ORDmax=ORDmax-1;
#        end    
    
#        if ORDmax==0; 
#            break;
#        end; % stop the code
    
        if counter==1:
            coef=55/24
        elif counter==2:
            coef=-1/6      
        elif counter==3:
            coef=11/8
        else:
            coef=1

#        integ=COEF*K^(d/2)*(besselj(d/2-1,K*R*N)).*(ones(length(R),1)*GK);
    
#        GRN=GRN+integ;
        

    grn=real(kperiod/((2*pi)^(d/2))*(r**(-d/2+1)*power(n,d/2+1))*grn)
    
    return grn

