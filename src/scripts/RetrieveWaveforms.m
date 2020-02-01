%% LTSpice 2 Matlab
%   *Put this in the script MIF3 Folder
% Author: JCCopyrights Autumn 2019
% Retrieves waveforms for a WPT system with Class D inverter from LTSpice
% Compare ZVS and no ZVS circuits

f=6.78e6;
T=1/f;
periods=100;
addpath('../functions');

%file_name='DBS.asc';
file_name='DBS_NoZVS.asc'; 
   
with_ZVS=true;

for cycles=1:1:1
    
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
	
	%LT_Izvs=raw.variable_mat(search_var(raw, 'I(L3)'),k:sim_length);
    %AUX Voltage
    %LT_Vaux=LT_Vlinkout-LT_Vrectin;
	
    
    %% Eliminates all NaN data that appear
	% Some errors during the importing can appear that can ruin the
	% integration with NaN values
	NaN_cnt=0;
	reduced_sim_length=length(LT_Pin);
	i=1;
    while(i<=reduced_sim_length-NaN_cnt)
		if isnan( LT_Pin(i))|isnan( LT_Pout(i))|isnan( LT_t(i))|isnan( LT_Prectin(i))|isnan( LT_Plinkin(i))|isnan( LT_Plinkout(i))
			NaN_cnt=NaN_cnt+1;
			 LT_Pin(i)=[];  LT_Pout(i)=[];  LT_t(i)=[]; %Eliminates the element from all vectors
			 LT_Plinkout(i)=[]; LT_Plinkin(i)=[]; LT_Prectin(i)=[];
             LT_Vlinkin(i)=[];   LT_Ilinkin(i)=[];
             LT_Irectin(i)=[];  LT_Ilinkout(i)=[];  LT_Iin(i)=[];
             LT_Vin(i)=[]; LT_Vlinkout(i)=[]; LT_Vout(i)=[]; LT_Vrectin(i)=[];  
             LT_Iout(i)=[];
             %LT_Vaux(i)=[]; LT_Izvs(i)=[]; %ZVS
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
save('../../data/LTsim_delete.mat')


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
