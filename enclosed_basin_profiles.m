function []= fillingbox_evolution()
% A function to do the Germeles (JFM 1975) calculation of a profile of a
% constant volume basin receiving buoyancy from an axisymmetric source. The
% function outputs plots of the evolving profile in the lake basin.
%
% Written by C. Hogg Date 2013_01_17
% Modified to fix bug spotted by megan on line 47 on 20170320

% Initiate time step sizes
total_time=5;
delta_tau=0.01;
total_time_steps=ceil(total_time/delta_tau);
% Set time steps to plot
index_snapshots=ceil([3 total_time/10/delta_tau 2*total_time/ ...
	10/delta_tau total_time/delta_tau])

% The first layer is calculated asymptotically
initial_plume_elements=1;
total_plume_elements=initial_plume_elements;
d_a=0;
zeta(1)=1;
f(1)=1;
q_2(1)=1;
m_2(1)=1;

% Do first time step to avoid crashing from having only one layer
total_plume_elements=total_plume_elements+1;
top_thickness=sqrt(q_2(end))*delta_tau;
top_d_a=d_a(end)+f(end)/(q_2(end)^(1/2));
zeta(1)=zeta(1)-top_thickness;
zeta(2)=1;
d_a(2)=top_d_a;

for ii=2:total_time_steps
    %% Calculate plume evolution
		
		% Calculate	plume fluxes at the top of the ponded region
		q_2(1)=zeta(1)^2;
    m_2(1)=zeta(1)^2;
    f(1)=1;

		% iterate through the layers of the plume
    for i=2:total_plume_elements
        delta_zeta=zeta(i)-zeta(i-1);
        % calculate terms for the Runge-Kutta scheme
        zeta_halfstep=delta_zeta/2;
        q_2_halfstep=q_2(i-1) + zeta_halfstep*(2*m_2(i-1)^(1/2));
        m_2_halfstep = m_2(i-1) + ...
						zeta_halfstep *(2*q_2(i-1)^(1/2)*f(i-1));
        f_halfstep = f(i-1) + zeta_halfstep * ...
				(-q_2(i-1)^(1/2)*(d_a(i)-d_a(i-1))/(delta_zeta));
        % calculate the next values in the iteration
        q_2(i) = q_2(i-1) + delta_zeta*(2*m_2_halfstep^(1/2));
        m_2(i) = m_2(i-1) + delta_zeta  ...
						*(2*q_2_halfstep^(1/2)*f_halfstep);
        f(i)   = f(i-1) + delta_zeta ...
		*(-q_2_halfstep^(1/2)*(d_a(i)-d_a(i-1))/(delta_zeta));
        
    end
    

		
    %%  Find new stratification
    % Calculate parameters for new layer
		total_plume_elements    =total_plume_elements+1;
		top_thickness           =sqrt(q_2(end))*delta_tau;
		top_d_a                 =d_a(end)+f(end)/q_2(end)^(1/2);
    
	% Plot snapshot
     if sum(ii==index_snapshots)>0
           f_han=figure;
           plot(sqrt(q_2),zeta(1:total_plume_elements-1));
           set(gca,'YDir','reverse');
           hold on
           plot(sqrt(m_2),zeta(1:total_plume_elements-1),'r');
           plot(f,zeta(1:total_plume_elements-1),'g');
           plot(d_a,zeta(1:total_plume_elements-1),'c');
           legend({'volume flux','momentum','buoyancy flux'...
					  ,'density profile'});
           ylabel('dimensionless height, \zeta')
           title(['t=',num2str(ii*delta_tau)])
           ylim([0 1])
     end


    %% Set new layer depths
    if ii ~= total_time_steps % Not on final iteration
		zeta(1:total_plume_elements-2)= ...
			zeta(1:total_plume_elements-2)- ...
		  sqrt(q_2(1:total_plume_elements-2))*delta_tau;
		zeta(total_plume_elements-1)=1-top_thickness;
		zeta(total_plume_elements)=1;
		d_a=[d_a,top_d_a];
    end
end