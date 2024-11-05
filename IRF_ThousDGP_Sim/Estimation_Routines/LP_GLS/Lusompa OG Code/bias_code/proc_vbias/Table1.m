clear;
clc;

%1.

%2. Declare and set the parameters.

global nz nq nk nN iflag g h i j phi1A phi2A phi1B phi2B phi1C phi2C sigma2_u phi_matrix omegaU phi PX Table1_vb_d Table1_ma0 Eetaztwid0 vbias q_vec delta Eetaz

nz=2;
nk=1;
nN=1;
iflag=0;

phi1A = 1.5;
phi2A = -0.54;
phi1B = 0.55;
phi2B = 0.4;
phi1C = 0.9;
phi2C = 0;
phi_matrix=[phi1A,phi2A;phi1B,phi2B;phi1C,phi2C];

q_vec=0;

sigma2_u=1;
omegaU=zeros(nz,nz);
omegaU(1,1)=sigma2_u;

PX=zeros(nk,nz);
PX(1,1)=1;

%3. Execute the program

Table1_vb_d=zeros(size(phi_matrix,1),size(q_vec,1));
Table1_ma0=zeros(size(phi_matrix,1),size(q_vec,1));

for g=1:size(phi_matrix,1)
    phi = [phi_matrix(g,:);1.0,0.0];
    nq=q_vec(1);

    %dim delta(nq+1) Eetaz(nq+1)
    delta=zeros(nz,1,nq+1);
    Eetaz=zeros(1,nz,nq+1);
    % compute delta
    for h=1:nq+1
        % defining Eetaz in this way sets delta=1
        Eetaz(1,1,h) = sigma2_u;
        delta(1,1,h)=1.0;
    end
    EetaZtwid0 = Eetaz(:,:,1);

    Table1_vb_d(g,1)=proc_vb_d(nz,phi,omegaU,nk,PX,delta,nq);
    Table1_ma0(g,1)=proc_vb_ma0(nz,phi,omegaU,nk,PX,EetaZtwid0);
end

display('Table 1B, col (2), using proc_vb_d');
Table1_vb_d
display('Table 1B, col (2), using proc_vb_ma0');
Table1_ma0
