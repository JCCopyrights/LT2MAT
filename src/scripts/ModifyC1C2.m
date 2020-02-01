%% LTSpice 2 Matlab
%   *Put this in the script MIF3 Folder
% Author: JCCopyrights Autumn 2019
% Retrieves waveforms for a WPT system with Class D inverter from LTSpice
f=6.78e6;
w=2*pi*f;
T=1/f;
periods=100;
addpath('../functions');
%file_name='DBS_AUXC2.asc';
file_name='DBS_noAUX.asc';
modify_C2=false;
lwdth=1.0;
L2=1.498e-6;
L1=1.497e-6;
C2=1/(w^2*L2)*1e12; %In pico farads
C1=1/(w^2*L1)*1e12; %In pico farads
if modify_C2==true
    RR=linspace(C2*0.5,C2*1.5,25);
else
    RR=linspace(C1*0.5,C1*1.5,25);
end
RR=round(RR); %avoid too many decimals
wait=waitbar(0,'Initialization');

for cycles=1:1:length(RR)
    if modify_C2==true
        text = sprintf('C2: %f', RR(cycles));
        waitbar(cycles/length(RR),wait,text);
        LTmodify( file_name, 'C2', [num2str(RR(cycles)) 'p']) ;
    else
        text = sprintf('C1: %f', RR(cycles));
        waitbar(cycles/length(RR),wait,text);
        LTmodify( file_name, 'C1', [num2str(RR(cycles)) 'p']) ;
        aux=1/(2*pi*sqrt(RR(cycles)*1e-12*L1));
        LTmodify_param( file_name, 'f', num2str(round(aux))); 
    end
    tic;
    raw=LTautomation(file_name);
    toc;
	sim_length=length(raw.time_vect);
    
	%% Only the last Periods will be integrated
	last=raw.time_vect(sim_length);
	ini=last;
	i=0;
	while last-ini<=periods*T
		ini=raw.time_vect(sim_length-i);
		i=i+1;
	end
	k=sim_length-i;
    
    %% Variable Searching
	% Searchs the varaibles in the raw data.
	% I recommend doing this manually the first time to get the voltage and current references and names fine.
	% Voltages are ALWAYS lower key, Currents keep the component name
    
	LT_t_original=raw.time_vect(k:sim_length);
    LT_t=LT_t_original;
	
    %Input
    LT_Vin=raw.variable_mat(search_var(raw, 'V(vin)'),k:sim_length);
	LT_Iin=-raw.variable_mat(search_var(raw, 'I(V1)'),k:sim_length);
    LT_Pin=LT_Vin.*LT_Iin;
    %Output
	LT_VoutP=raw.variable_mat(search_var(raw, 'V(vout+)'),k:sim_length);
	LT_VoutN=raw.variable_mat(search_var(raw, 'V(vout-)'),k:sim_length);
	LT_Vout=LT_VoutP-LT_VoutN;
	LT_Iout=raw.variable_mat(search_var(raw, 'I(Load)'),k:sim_length);
	LT_Pout=LT_Vout.*LT_Iout;
    %Output of the Inverter-Input of the Link
	LT_Ilinkin=-raw.variable_mat(search_var(raw, 'I(C1)'),k:sim_length);
	LT_Vlinkin=raw.variable_mat(search_var(raw, 'V(vsw)'),k:sim_length);
	LT_Plinkin=LT_Ilinkin.*LT_Vlinkin;
    %Output of the Link-Input of the AUX
	LT_Vlinkout=raw.variable_mat(search_var(raw, 'V(vaux+)'),k:sim_length);
	LT_Ilinkout=-(raw.variable_mat(search_var(raw, 'I(R2)'),k:sim_length));%+raw.variable_mat(search_var(raw, 'I(Cp2)'),k:sim_length));
	LT_Plinkout=LT_Ilinkout.*LT_Vlinkout;
	%Output of the AUX-Input of the Rectifier
	LT_Vrectin=raw.variable_mat(search_var(raw, 'V(vrect)'),k:sim_length);
	LT_Irectin=(raw.variable_mat(search_var(raw, 'I(D1)'),k:sim_length)-raw.variable_mat(search_var(raw, 'I(D3)'),k:sim_length));
	LT_Prectin=LT_Irectin.*LT_Vrectin;
	%Current in the ZVS Network
	LT_Izvs=raw.variable_mat(search_var(raw, 'I(L3)'),k:sim_length);
    %AUX Voltage
    LT_Vaux=LT_Vlinkout-LT_Vrectin;
    
    %% Eliminates all NaN data that appear
	% Some errors during the importing can appear that can ruin the
	% integration with NaN values (@TODO: Maybe a codification problem?)
	NaN_cnt=0;
	reduced_sim_length=length(LT_Pin);
	i=1;
    while(i<=reduced_sim_length-NaN_cnt)
		if isnan( LT_Pin(i))|isnan( LT_Pout(i))|isnan( LT_t(i))|isnan( LT_Prectin(i))|isnan( LT_Plinkin(i))|isnan( LT_Plinkout(i))
			NaN_cnt=NaN_cnt+1;
			 LT_Pin(i)=[];  LT_Pout(i)=[];  LT_t(i)=[]; %Eliminates the element from all vectors
			 LT_Plinkout(i)=[]; LT_Plinkin(i)=[]; LT_Prectin(i)=[];
             LT_Vlinkin(i)=[];  LT_Izvs(i)=[]; LT_Ilinkin(i)=[];
             LT_Irectin(i)=[];  LT_Ilinkout(i)=[]; LT_Vaux(i)=[]; LT_Iin(i)=[];
             LT_Vin(i)=[]; LT_Vlinkout(i)=[]; LT_Vout(i)=[]; LT_Vrectin(i)=[];  
             LT_Iout(i)=[];
        else
			i=i+1; %It doesnt increse if a NaN is found
		end
    end
    %% Power in every stage
	LT_Pinv(cycles)=trapz(LT_t,LT_Pin)/(periods*T);
	LT_Poutv(cycles)=trapz(LT_t,LT_Pout)/(periods*T);
	LT_Plinkoutv(cycles)=trapz(LT_t,LT_Plinkout)/(periods*T);
	LT_Plinkinv(cycles)=trapz(LT_t,LT_Plinkin)/(periods*T);
	LT_Prectinv(cycles)=trapz(LT_t,LT_Prectin)/(periods*T);
	LT_efic_total(cycles)=LT_Poutv(cycles)/LT_Pinv(cycles);
	LT_efic_inverter(cycles)=LT_Plinkinv(cycles)/LT_Pinv(cycles);
	LT_efic_link(cycles)=LT_Plinkoutv(cycles)/LT_Plinkinv(cycles);
	LT_efic_aux(cycles)=LT_Prectinv(cycles)/LT_Plinkoutv(cycles);
	LT_efic_rect(cycles)=LT_Poutv(cycles)/LT_Prectinv(cycles);
    clear('raw');
end
if modify_C2==true
    LTmodify( file_name, 'C2', '{C2}');
else
    LTmodify( file_name, 'C1', '{C1}');
    LTmodify_param( file_name, 'f', '6.78Meg'); 
end

delete(wait)

figure();
hold on
grid on;
plot(RR/C2*100-100,LT_efic_total,'LineWidth',lwdth);
plot(RR/C2*100-100,LT_efic_inverter,'LineWidth',lwdth);
plot(RR/C2*100-100,LT_efic_link,'LineWidth',lwdth);
plot(RR/C2*100-100,LT_efic_aux,'LineWidth',lwdth);
plot(RR/C2*100-100,LT_efic_rect,'LineWidth',lwdth);
legend('\eta total','\eta inverter','\eta link','\eta aux', '\eta rect')
title('\eta')
if modify_C2==true
    xlabel('\DeltaC_2%')
else
    xlabel('\DeltaC_1%')
end


figure();
grid on;
yyaxis left 
plot(RR/C2*100-100,LT_efic_total,'LineWidth',lwdth);
ylabel('\eta')
yyaxis right 
plot(RR/C2*100-100,LT_Poutv,'LineWidth',lwdth);
ylabel('P_o_u_t [W]')
title('\eta vs C_2')

if modify_C2==true
    xlabel('\DeltaC_2%')
else
    xlabel('\DeltaC_1%')
end
delete(wait)

save('../../data/LTsim_ModifyC2NOAUX.mat')


%% Search Var
% Searches variables in the raw data by its name
function index=search_var(raw, variable_name)
	for index=1:1:length(raw.variable_name_list)
		if strfind(raw.variable_name_list{index}, variable_name)%CASE SENSITIVE!!
			return;
        end
    end
    error('Component not found');
    %@TODO: Handle strfind exCEPTIONS
end
