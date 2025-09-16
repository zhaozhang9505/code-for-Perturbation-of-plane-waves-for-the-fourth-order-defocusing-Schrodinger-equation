function dnvf=DNLSE4_zu(t,nvf,dummy,k,N)
k1=0;
nvf(isnan(nvf));
nf=nvf(1:N);
vf=nvf(N+1:2*N);
n=ifft(nf);
v=ifft(vf);

n1f=1i*k.*nf;
n2f=-k.^2.*nf;
n3f=-1i*k.^3.*nf;
n4f=k.^4.*nf;
n5f=1i*k.^5.*nf;
n1=ifft(n1f);
n2=ifft(n2f);
n3=ifft(n3f);
n4=ifft(n4f);
n5=ifft(n5f);

v1f=1i*k.*vf;
v2f=-k.^2.*vf;
v3f=-1i*k.^3.*vf;
v4f=k.^4.*vf;
v5f=1i*k.^5.*vf;
v1=ifft(v1f);
v2=ifft(v2f);
v3=ifft(v3f);
v4=ifft(v4f);
v5=ifft(v5f);


 dnf= fft(-4*v.*n3-2*n.*v3+6*n2.*v.*n1./n-6*n2.*v1...
     -4*n1.*v2-3*v.*n1.^3./n.^2+3*v1.*n1.^2./n...
     +4*n1.*v.^3+24*n.*v.*n1+12*n.*v1.*v.^2+12*n.^2.*v1);


dvf=fft(4*v.^3.*v1+12*n1.*v.^2+12*n.*n1-4*v3.*v...
    -10*v2.*v1+1/2./n.*n5+24*n.*v.*v1-3./n.*n3.*v.^2 ...
    -3/2./n.^2.*n4.*n1-5/2./n.^2.*n3.*n2+6./n.^3.*n2.^2.*n1...
    +15/4./n.^5.*n1.^5+9./n.^2.*v1.*v.*n1.^2+6./n.^2.*n2.*v.^2.*n1...
    -12./n.*n2.*v.*v1-6./n.*v2.*v.*n1-6./n.*v1.^2.*n1...
    -3./n.^3.*v.^2.*n1.^3+17/4./n.^3.*n3.*n1.^2-21/2./n.^4.*n2.*n1.^3)-5.*n3f;


dnvf=[dnf;dvf];

end