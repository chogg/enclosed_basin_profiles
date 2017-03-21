function []= dimensionless_front_evolution_const_E_rms_err_for_paper()
% A function to plot all data collapsed onto the theoretical line.
%
% Inputs:
%		inputname:
% Outputs:
%		name:  description type units
%		saved data: (does this overwrite a statically named file?)
%		plots:
%
% Standard call:
%
% Written by C. Hogg Date 2013_12_20
% Modified from dimensionless_experimental_results_20131218b
% ...
% modified to calculate root mean square error.

    debugmode=2; % plot
    debugmode=1;

if debugmode==1; dbstop if error;end

% Set parameters
set_standard_variables
exp_data=calculate_derived_experimental_data;

% Load raw data
meta_data_filename='C:\Users\carh5\laboratory_experiments\Scripts\experiment_meta_data\exp_meta_data_current.dfc';
meta_data = DF_dfc_variable_read(meta_data_filename);
height_to_source_from_weir_base=(meta_data.tank_data.total_length-meta_data.tank_data.partition_to_wall-meta_data.tank_data.partition_thickness-meta_data.tank_data.source_to_wall)*sin(meta_data.exps_data{1}.slope_angle*pi/180); % vertical height

% Get the series of experiments at constant inclination and g'
[valid_exp_indices]= get_valid_exp_id_only_vary_Q(exp_data.exps_data_array)
% Get rid of laminar and transitional regime runs
valid_exp_indices([exp_data.exps_data_array(valid_exp_indices).source_flow_rate]<0.05)=[]; %Exclude laminar and transitional flows



[val Q_order]=sort(-[exp_data.exps_data_array(valid_exp_indices).source_flow_rate]);

fig_h1=figure;
cmap=colormap;
hold all

for exp_index=valid_exp_indices(Q_order)
        exp_front_data=analyse_exp_ponded_height(exp_index,'concentration_threshold',0.1,'median_window_size',3,...
        'run_av_window_size',1,'interface_method','first_sample_above_threshold_from_bottom','max_density_method','global_threshold');
    time=exp_front_data.time;
    height=exp_front_data.height;

    Q(exp_index)=exp_data.exps_data_array(exp_index).source_flow_rate/exp_data.tank_data.width/l_per_m3; % in m^2 s^{-1}
    %No longer used E=6.8E-5/1.05E-6*Q-0.02; % From dissertation text using previous ponded region algorithm.
    
    E=62.0*Q(exp_index)-0.022;% from nonlinear_fit_entrainment_1_param_for_paper. see 2014-09-10T17:07:45 phd_notebook 
    
    F_0=Q(exp_index)*meta_data.exps_data{exp_index}.g_prime;
    time_nd{exp_index}=time*E^(2/3)*F_0^(1/3)*cosd(meta_data.exps_data{exp_index}.slope_angle)*(sind(meta_data.exps_data{exp_index}.slope_angle))^(1/3)/height_to_source_from_weir_base;
    
    D_s(exp_index)=height_to_source_from_weir_base+Q(exp_index)*sind(meta_data.exps_data{exp_index}.slope_angle)^(2/3)/(E)^(2/3)/F_0^(1/3);
    height_nd{exp_index}=height/D_s(exp_index);
    
    Re(exp_index)=exp_data.derived_data(exp_index).Re;
    max_Re=max([exp_data.derived_data(valid_exp_indices).Re]);
    colour_spec=cmap(floor(Re(exp_index)/max_Re*length(cmap)),:);
    
    figure(fig_h1)
    plot(time_nd{exp_index},height_nd{exp_index},'Color',colour_spec)
    
    clear traverse_ind density_profile density_profiles max_density min_density normalised_density height_ind ponding_interface height_from_weir_a height_from_weir_b height_from_weir first_traverse_time_delay time
    
end

% Tidy
xlabel('dimensionless time, $\tau$','interpreter','none')
ylabel('dimensionless height, $\eta$','interpreter','none')
max_x=2;
xlim([0 max_x]);
ylim([0 1]);

% Plot theory

figure(fig_h1)
theory_t = 0:0.01:max_x;
theory_h = 1-(-lambertw(-exp(-theory_t-1)));
plot(theory_t,theory_h,'k--','LineWidth',3)

% show color bar at specified location
caxis([0,Q(valid_exp_indices(Q_order(1)))]*1000)
cbar_hand=colorbar;
pos=get(cbar_hand,'position');
set(cbar_hand,'position',pos+[-0.01 0  -0.02 0])

pos2=get(gca,'position');
pos2(3)=pos(1)-pos2(1)-0.05;
set(gca,'position',pos2())
add_text_on_plot('$Q_s$ [m$^{2}$ s$^{-1}$]',0.85,0.5);
add_text_on_plot('$\times 10^{-3}$',0.89,0.95);

% Calculate standard deviation
obs_error_all=[];
for exp_index=valid_exp_indices(Q_order)
    theory_height{exp_index}=1-(-lambertw(-exp(-time_nd{exp_index}-1)));
    obs_error{exp_index}=height_nd{exp_index}-theory_height{exp_index};
    obs_error_all=[obs_error_all,obs_error{exp_index}];
end
root_mean_sq_error=sqrt(mean(obs_error_all.^2))
%scatter(1,(root_mean_sq_error));

% Print for thesis
if debugmode==2;
    filename=['dimensionless_front_evolution_variable_E_rms_err_for_paper_output_',datestr(now,30)];
    printCH_for_thesis(gcf,filename,'C:\Users\carh5\matlab\paperHDHI_plot_scripts\','BottomMargin',20,'RightMargin',15,'LeftMargin',10)
end

%% Debug. Break point here.
if debugmode==1|2; keyboard;end

end
