function demo_Leaspy_param

%long-term parameters
t0=75;
% pk=\gamma_k(t0)
% uk= d/dt(\gamma_k(t0))
pk_1=0.5;%A13
uk_1=0.02; % soft
pk_2=0.2;%CDR
uk_2=0.2; % edge

t=60:0.1:90;

%i-subject
tau_i=2;%cognitive reserve in years
alpha_i=0.5;% acceleration greather than 1
wik_1=0.01;
wik_2=0.1;

tij=70:1:75;

gamma_k_1=get_gamma_k(t,t0,pk_1,uk_1); %long term A13
gamma_k_2=get_gamma_k(t,t0,pk_2,uk_2); %long term CDR

eta_k_1=get_eta_k(tij,t0,pk_1,uk_1,tau_i,alpha_i,wik_1); %short term
eta_k_2=get_eta_k(tij,t0,pk_2,uk_2,tau_i,alpha_i,wik_2); %short term

figure(1);clf;
plot(t,gamma_k_1,'b-');hold on;
plot(t,gamma_k_2,'r-');
xlabel('ages');
ylabel('Percentages');
legend('ADAS13','CDRSB');

plot(tij,eta_k_1,'b+');
plot(tij,eta_k_2,'r+');hold off;


end

function gamma_k=get_gamma_k(t,t0,pk,uk)
gamma_k=1./(1+(((1/pk)-1).*exp(-uk*(t-t0)/(pk*(1-pk)))));
end

function eta_k=get_eta_k(tij,t0,pk,uk,tau_i,alpha_i,wik)
eta_k=1./(1+(((1/pk)-1).*exp(-(wik+(uk*alpha_i*(tij-t0-tau_i)))/(pk*(1-pk)))));
end
