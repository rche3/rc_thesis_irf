function vbias = proc_vb_d(nZtwid,phitwid,omegaUtwid,nk,PX,delta,nq,vbias)

% local variables
% vect[real] EXeta vbias1 vbias2
% vect[rect] Eetaztwid
% rect gamma0 bigD
% int nN iflag

% parameters
% type int nztwid nq nk
% type rect[real] omegaUtwid phitwid PX
% type vector[real] *vbias
% type vector[rect] delta

% dimension of the parameters
% dim vbias1(nk) vbias2(nk) gamma0(nztwid,nztwid) bigD(nk,nk)
% dim EXeta(nk) Eetaztwid(nq+1)

nN=1;
iflag=0;

% just a declaration, not used when iflag=0
Eetaztwid=zeros(nN,nN*nZtwid,nq+1);
EXeta=zeros(nk,nN*nN*nZtwid);
bigD=zeros(nk,nk);
gamma0=zeros(nZtwid,nZtwid);
vbias1=zeros(nk,1);
vbias2=zeros(nk,1);

[bigD,gamma0,vbias1,vbias2,vbias]=proc_vbias(iflag,nZtwid,phitwid,omegaUtwid,nq,nk,nN,PX,Eetaztwid,EXeta,delta,bigD,gamma0,vbias1,vbias2,vbias);
    
end
