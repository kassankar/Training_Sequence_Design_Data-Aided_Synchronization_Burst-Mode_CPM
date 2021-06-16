function MAIN = MainFunctions
MAIN.CREATECPMPULSE                    = @CREATECPMPULSE;
end

function [g_t,q_t,t] = CREATECPMPULSE(pulse,pulse_length,w,freq_sa,img)
g_t              = 0;
q_t              = 0;
pulse_length_sa  = pulse_length*freq_sa;
time_sa          = 1/freq_sa;
switch             pulse
    
    case    1
t          = -pulse_length/2:time_sa:pulse_length/2;
t0         = 0;
g_t        = (2*w)./((t-t0).^2+w^2);              % Fix form w --> to w^2 (the error is compensated in the correcting factor)
Cst        = sum(g_t)*time_sa;
nug_t      = (2*pi) / Cst;                        
g_t        = nug_t*g_t;
g_t        = g_t/(4*pi);
q_t        = cumtrapz(g_t)*time_sa;        % So we can obtain the same CPM forme. (0 t<0, 1/2 t>LT).


    case    2
t           = -pulse_length/2:time_sa:pulse_length/2;
% Tb        = pulse_length*time_sa;
BT          = 0.3;
% B         = BT/Tb;
% h         = sqrt(2*pi*B^2/(log(2)))*exp(-t.^2*2*pi^2*B^2/(log(2)));
% h         = h/sum(h);
% g_t       = h;
% K         = 0.5/trapz(t,g_t);
% g_t       = K*g_t;
alpha       = 2*pi*BT/(sqrt(log(2)));
gauss       = qfunc(alpha*(t-0.5)) - qfunc(alpha*(t+0.5));
Cst         = 0.5/(sum(gauss)*time_sa); 
g_t         = Cst*gauss;
q_t         = cumtrapz(g_t)*time_sa;


    case    3
t          = 0:time_sa:pulse_length;
g_t        = (1/(2*pulse_length).*(1- cos(2*pi.*t/(pulse_length))));
K          = 0.5/(sum(g_t)*time_sa);
g_t        = K*g_t;
q_t        = cumtrapz(g_t)*time_sa;  


    case   4       
t                      = 0:time_sa:pulse_length;
g_t                    = 1/(2*pulse_length)*ones(1,length(t));
g_t(1)                 = 0;
% g_t(end)               = 0;
K                      = 0.5/(sum(g_t)*time_sa);
g_t                    = K*g_t;
q_t                    = cumtrapz(g_t)*time_sa;
    case 5

t          = -pulse_length/2:time_sa:pulse_length/2;
t0         = 0;
g_t        = (2*w)./((t-t0).^2+w^2);              % Fix form w --> to w^2 (the error is compensated in the correcting factor)
g_t        = g_t*(1/(4*pi));
Cst        = sum(g_t)*time_sa;
nug_t      = (0.5) / Cst;
g_t        = nug_t*g_t;
q_t        = cumtrapz(g_t)*time_sa; 


    case 6 
t          = -pulse_length/2:time_sa:pulse_length/2;
t0         = 0;
g_t        = (2*w)./((t-t0).^2+w^2);
g_t        = g_t*(1/(4*pi));
Cst        = sum(g_t)*time_sa;
nug_t      = (0.5) / Cst;
g_t        = nug_t*g_t;
q_t        = cumtrapz(g_t)*time_sa;

    case 7
t          = -pulse_length/2:time_sa:pulse_length/2;
t0         = 0;
g_t        = (2*w)./((t-t0).^2+w^2);
g_t        = g_t*(1/(4*pi));



s         = 2.7;
gauss      = exp(-(t.^2)./(2*s^2));




g_t        = g_t.*gauss;
Cst        = sum(g_t)*time_sa;
nug_t      = (0.5) / Cst;
g_t        = nug_t*g_t;
q_t        = cumtrapz(g_t)*time_sa;

    case 8
t          = -pulse_length/2:time_sa:pulse_length/2;
t_1        = (0:pulse_length_sa)*time_sa;
t0         = 0;
g_t        = (2*w)./((t-t0).^2+w^2);
g_t        = g_t*(1/(4*pi));



BT        = 0.3;
s         = 2.7;
alpha       = 2*pi*BT/(sqrt(log(2)));
% gauss       = qfunc(alpha*(t-0.5)) - qfunc(alpha*(t+0.5));
gauss      = exp(-(t.^2)./(2*s^2));
rec        = 1/(2*pulse_length*1).*(1- cos(2*pi.*t_1/(pulse_length*1)));


g_t        = g_t.*gauss;
Cst        = sum(g_t)*time_sa;
nug_t      = (0.5) / Cst;
g_t        = nug_t*g_t;
q_t        = cumtrapz(g_t)*time_sa;
    case 9

t          = 0.0000000001:time_sa:pulse_length/2;
g_t        = 1/(pulse_length)*(sin(2*pi*t/(pulse_length))./(2*pi*t/pulse_length)).*(cos(2*pi*t/pulse_length)./(1-4*(2*t/pulse_length).^2));
g_t        = [flip(g_t) g_t(1) g_t];
K          = 0.5/(sum(g_t)*time_sa);
g_t        = K*g_t;
q_t        = cumtrapz(g_t)*time_sa;

% t           = [flip(t) t];
    case 10
t          = -pulse_length/2:time_sa:pulse_length/2;
t0         = 0;
g_t        = (2*w)./((t-t0).^2+w^2);              % Fix form w --> to w^2 (the error is compensated in the correcting factor)
g_t        = g_t*(1/(4*pi));
s          = 1.37;
gauss      = exp(-(t.^2)./(2*s^2));
gauss      = gauss*(1/(4*pi));
g_t        = gauss.*g_t;
% Cst        = 0.5/(sum(g_t)*time_sa); 
% g_t        = Cst*gauss;
Cst        = sum(g_t)*time_sa;
nug_t      = (0.5) / Cst;
g_t        = nug_t*g_t;
q_t        = cumtrapz(g_t)*time_sa;

case 11
   t          = 0:time_sa:pulse_length;
   g_t        = rcosdesign(1,pulse_length,freq_sa);
   cst        = 0.5/(sum(g_t)*time_sa);
   g_t        = cst*g_t;
   q_t        = cumtrapz(g_t)*time_sa;  


    otherwise
        disp('no pulses for >4')
end



if(img ==1)
set(gca,'FontName','Arial','FontSize',12);
figure (1)
% subplot(1,2,1)
if(pulse==1)
    plot(t,g_t*4*pi)
else
    plot(t,g_t)
end

xlabel("time(t/T_s)")
ylabel("Freqeuncy pulse g(t)")
grid on;
set(gca,'FontName','Arial','FontSize',12)
% subplot(1,2,2)
figure(2)
plot(t,q_t)
xlabel("time(t/T_s)")
ylabel("Phase response q(t)")
grid on;
set(gca,'FontName','Arial','FontSize',12);
else
end


end



