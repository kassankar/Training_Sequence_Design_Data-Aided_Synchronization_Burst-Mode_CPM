clc
clear
close all
%% Pulse shape & Variable ini
MAIN          = MainFunctions;
pulse         = 2;      % 1 -> lorentzian pulse
                        % 2 -> GMSK pulse BT = 0.3
                        % 3 -> LRC pulse
                        % 4 -> LREC pulse
L             = 4;      % Pulse length
                        % 1  -> Full response
                        % >1 -> Partial response
os            = 2^4;    % Over sampling frequency
Ts            = 1/os;   % Sampling Time
M             = 2^1;    % M_ary symbols used (2 -> Binary)
h             = 0.5;    % Modulation index
width         = 0.6;    % This variable is used for Lorentzian Pulse only. (Not be used for pulse > 1)
%% -------------- Modulated data ----------------------
snr       = 0:10;
snr_lin   = 10.^(snr/10);
%% frequency pulse
[g_t,q_t]     = MAIN.CREATECPMPULSE(pulse,L,width,os,0); % Function return the CPM pulse and phase.
                                                         % g_t = g(t) is the CPM pulse shape.
                                                         % q_t -> is the phase, integral of g_t.
%-------------------------------------------------------------------------------------------------------------
L0       = 64; % preamble length
bits     = [-1*ones(1,L0/4) ones(1,L0/2) -1*ones(1,L0/4)]; % preamble
%--------------------- A B C calculation -----------------------
%% FIM calculation
Nbits             = L0;
bits_NO           = bits;
A                 = 0;
B                 = 0;
C                 = 0;
t_seq             = 0:Ts:L0-Ts;
dPhi_tauxx        = zeros(Nbits,(Nbits+1)*os+length(g_t)-1);
dPhi_tauxx_1      = zeros(Nbits,length(t_seq));
index_min         = 1;
index_max         = index_min+length(g_t)-1;
for i= 1:(Nbits)
    dPhi_tauxx(i,index_min:index_max) = g_t;
    dPhi_tauxx(i,1:index_min-1)       = 0;
    dPhi_tauxx(i,index_max+1:end)     = 0;
    index_min                         = index_min+os;
    index_max                         = index_max+os;
    dPhi_tauxx_1(i,1:end)             = dPhi_tauxx(i,1:length(t_seq)); %(different between time start_Keep)
    A                                 = A + bits_NO(i)*sum(dPhi_tauxx_1(i,1:end).*t_seq)*Ts;
    B                                 = B + bits_NO(i)*sum(dPhi_tauxx_1(i,1:end))*Ts;
end

for j= 1:(Nbits)
    for i= 1:(Nbits)
        C                                  = C + bits_NO(j)*bits_NO(i)*sum(dPhi_tauxx_1(j,1:end).*dPhi_tauxx_1(i,1:end))*Ts;
    end
end
B = 2*B;
%------------------------------------------------------------
%------------------------ CRB equations -----------------------------------
Tb          = 1;
T0          = L0*Tb;
IN_3        = T0^3/(8*pi^2*h^2*(C*T0^3-(3*A^2+(3*A-B*T0)^2))).* (Tb./snr_lin); %time
IN_1        = (3/(2*pi^2*T0^3))*((C*T0-B^2/4)/(C*T0-B^2/4-(3/(4*T0^2))*(B*T0-4*A)^2)).* (Tb./snr_lin); %freq
IN_2        = (2/T0)*((C*T0^3-3*A^2)/(C*T0^3-3*A^2-(B*T0-3*A)^2)).* (Tb./snr_lin); %phase
%------------------------------------------------------------
%---------------------- plot ----------------------
figure(1)
hold on;
plot(snr,IN_3);
set(gca, 'YScale', 'log');
title('Time offset')
xlabel('E_{b}/N_{0}')
ylabel('CRB (\tau)')
set(gca,'FontName','Arial','FontSize',12);
grid on;
box on;
figure(2)
hold on;
plot(snr,IN_1);
title('Frequency offset')
xlabel('E_{b}/N_{0}')
ylabel('CRB (f_{d})')
set(gca, 'YScale', 'log');
set(gca,'FontName','Arial','FontSize',12);
grid on;
box on;
figure(3)
hold on;
plot(snr,IN_2);
title('Phase offset')
xlabel('E_{b}/N_{0}')
ylabel('CRB (\theta)')
set(gca, 'YScale', 'log');
set(gca,'FontName','Arial','FontSize',12);
grid on;
box on;

















