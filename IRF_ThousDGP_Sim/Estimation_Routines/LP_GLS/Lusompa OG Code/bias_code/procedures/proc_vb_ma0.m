function vbias = proc_vb_ma0(nZtwid,phitwid,omegaUtwid,nk,PX,EetaZtwid0)

% Note: This is for special case nq=0, nN=1, and iflag=1.
% nq is the order of MA of eta. See documentation on proc_vbias for details.
% nN, iflag: See documentation on proc_vbias for explanations.
nq=0;
nN=1;
iflag=1;

% EXeta is an nk*1 vector where its elements are all zero for this special case of nq=0.
EXeta = zeros(nk,1);

% EetaZtwid0 is passed in as a nN*nZtwid matrix, which should be translated
% into a nN*nZtwid*(nq+1) matrix Eetaztwid before it can be passed into proc_vbias.
Eetaztwid=zeros(nN,nZtwid,nq+1);
for i=1:nN
    for j=1:nZtwid
        Eetaztwid(i,j,nq+1)= EetaZtwid0(i,j);
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