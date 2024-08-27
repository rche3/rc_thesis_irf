function vbias = proc_vb_maq(nZtwid,phitwid,omegaUtwid,nk,PX,EetaZtwid,EXeta,nq)

% Note: This is for special case nN=1 and iflag=1.
% nN, iflag: See documentation on proc_vbias for explanations.
nN=1;
iflag=1;

% EetaZtwid is passed in as a (nq+1)*nZtwid matrix, which should be translated
% into a nN*nZtwid*(nq+1) matrix Eetaztwid before it can be passed into proc_vbias.
Eetaztwid=zeros(nN,nZtwid,nq+1);
for i=1:nZtwid
    for j=1:nq+1
        Eetaztwid(nN,i,j)=EetaZtwid(j,i);
    end
end

% Declarations
% delta, bigD, gamma0, vbias1, vbias2: See documentation on proc_vbias for explanations.
delta=zeros(nZtwid,nN,nq+1);
bigD=zeros(nk,nk);
gamma0=zeros(nZtwid,nZtwid);
vbias1=zeros(nk,1);
vbias2=zeros(nk,1);

[bigD,gamma0,vbias1,vbias2,vbias]=proc_vbias(iflag,nZtwid,phitwid,omegaUtwid,nq,nk,nN,PX,Eetaztwid,EXeta,delta,bigD,gamma0,vbias1,vbias2);
end