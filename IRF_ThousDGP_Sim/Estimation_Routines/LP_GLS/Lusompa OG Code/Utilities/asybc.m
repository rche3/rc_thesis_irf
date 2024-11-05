% ASYBC.M
% Lutz Kilian
% University of Michigan
% April 1997
%
% Source: Pope (1990), JTSA

function [bcA]=asybc(A,SIGMA,t,p,K)

T=t-p;
vecSIGMAY=inv(eye((K*p)^2)-kron(A,A))*vec(SIGMA);
SIGMAY=reshape(vecSIGMAY,K*p,K*p); 
I=eye(K*p,K*p);
B=A';

% There are q*p eigenvalues by construction
peigen=eig(A);
sumeig=zeros(K*p,K*p);
for h=1:K*p
	sumeig = sumeig + (peigen(h).*inv(I-peigen(h)*B));
end

bias=SIGMA*(inv(I-B)+B*inv(I-B^2)+sumeig)*inv(SIGMAY);

Abias=-bias/T;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bcstab=9;         % Arbitrary default value > 1
delta=1;          % Adjustment factor

while bcstab >= 1

	% Adjust bias-correction proportionately
	bcA=A-delta*Abias;

	bcmod=abs(eig(bcA));
		if any(bcmod>=1)
			bcstab=1;
		else
			bcstab=0;
        end
	delta=delta-0.01;
	
	if delta <= 0
		bcstab=0;
    end
end