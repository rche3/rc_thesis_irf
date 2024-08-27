function vbias = proc_vb_d(nZtwid,phitwid,omegaUtwid,nk,PX,delta,nq)

% Note: This is for special case nN=1 and iflag=0.
% nN, iflag: See documentation on proc_vbias for explanations.
nN=1;
iflag=0;

% delta is passed in as a nZtwid*(nq+1) matrix, which should be translated
% into a nZtwid*nN*(nq+1) matrix Delta before it can be passed into proc_vbias.
Delta=zeros(nZtwid,nN,nq+1);
for i=1:nZtwid
    for j=1:nq+1
        Delta(i,nN,j)=delta(i,j);
    end
end

% Declarations
% delta, bigD, gamma0, vbias1, vbias2: See documentation on proc_vbias for explanations.
Eetaztwid=zeros(nN,nN*nZtwid,nq+1);
EXeta=zeros(nk,nN*nN*nZtwid);
bigD=zeros(nk,nk);
gamma0=zeros(nZtwid,nZtwid);
vbias1=zeros(nk,1);
vbias2=zeros(nk,1);

[bigD,gamma0,vbias1,vbias2,vbias]=proc_vbias(iflag,nZtwid,phitwid,omegaUtwid,nq,nk,nN,PX,Eetaztwid,EXeta,Delta,bigD,gamma0,vbias1,vbias2);
end