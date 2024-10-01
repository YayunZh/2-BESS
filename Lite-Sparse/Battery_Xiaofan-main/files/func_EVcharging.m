function TrajCont = func_EVcharging(OutVar, TrajENV, EV_demand, T_delay)
   TrajENV.EV_demand = EV_demand;  % Add ths to make parfor work
   TrajENV.T_delay = T_delay;      % Add ths to make parfor work
%     OutVar_trad = func_dc_energyflow_two_layers(ENV_trad.Avg_Conv.e_lim_mat, Bat_Info_MC(:,k2), CL, ENV_trad);
%     var_energy_flow(k0,k1,k2,:,:) = OutVar.Var_energy_flow;
%     avg_energy_flow(k0,k1,k2,:,:) = OutVar.Avg_energy_flow;
% 
%     var_energy_process(k0,k1,k2) = OutVar.Var_energy_process;
%     total_energy_process(k0,k1,k2) = OutVar.Total_energy_process;
%     avg_energy_process(k0,k1,k2) = OutVar.Avg_energy_process;
% 
%     bat_output_energy_individual(k0,k1,k2,:) = OutVar.BatOutputEnergy;
%     var_conv_energy_rating(k0,k1,k2,:) = OutVar.Var_conv_energy_rating;
%     avg_conv_energy_rating(k0,k1,k2,:) = ENV.Avg_Conv.e_lim_vec;
    TrajCont = InitializeTraj();
    bat_energy_rating = OutVar.Bat_energy_rating;
    max_output_energy = abs(OutVar.Maximum_output_energy);          
    bat_uratio_energy = abs(OutVar.Maximum_output_energy)/sum(OutVar.Bat_energy_rating);
    TrajCont.Ebess.value(1) = sum(bat_energy_rating);
    interval_ind = 1;
    while (TrajCont.Pgrid_avg.abs_time(end) <= 24)
        % Averaged Grid Power
        T_cycle = Discharge_time(TrajENV.EV_demand(interval_ind), TrajENV.Pgrid_fit, TrajCont.Pgrid_avg.abs_time(end));
        TrajCont.Pgrid_avg.value(end + 1) = TrajENV.EV_demand(interval_ind)/T_cycle;
        TrajCont.Pgrid_avg = UpdateTimeIndividual(TrajCont.Pgrid_avg, T_cycle);

        T_discharge_max = max_output_energy/(TrajENV.Pcharge_level3 - TrajCont.Pgrid_avg.value(end)); 
        if (T_discharge_max*TrajENV.Pcharge_level3 < TrajENV.EV_demand(interval_ind))  % EV demand is large, require grid to charge EV
        % Phase 1: BESS and grid charge EV together
            T_discharge_ph1 = T_discharge_max;

            TrajCont.Pbessout.value(end + 1) = TrajENV.Pcharge_level3 - TrajCont.Pgrid_avg.value(end);
            TrajCont.Ebess.value(end + 1) = sum(bat_energy_rating) - max_output_energy;
            TrajCont.Pevin.value(end + 1) = TrajENV.Pcharge_level3;
            
            TrajCont.Ebessout.value(end + 1) = max_output_energy;

            TrajCont = UpdateTime(TrajCont, T_discharge_ph1);
        % Phase 2: grid charge EV
            T_discharge_ph2 = (TrajENV.EV_demand(interval_ind)-TrajENV.Pcharge_level3*T_discharge_max)/TrajCont.Pgrid_avg.value(end);

            TrajCont.Pbessout.value(end + 1) = 0;
            TrajCont.Ebess.value(end + 1) = TrajCont.Ebess.value(end);   % Battery energy stays the same
            TrajCont.Pevin.value(end + 1) = TrajCont.Pgrid_avg.value(end);
            
            TrajCont.TotalChargeTime.value(end + 1) = T_discharge_ph1 + T_discharge_ph2;
            TrajCont.L3ChargeTime.value(end + 1) = T_discharge_ph1;
            TrajCont.CurtChargeTime.value(end + 1) = T_discharge_ph2;

            TrajCont = UpdateTime(TrajCont, T_discharge_ph2);

        % Phase 3: grid charge BESS
            T_charge_ph3 = (sum(bat_energy_rating) - TrajCont.Ebess.value(end))/TrajCont.Pgrid_avg.value(end);

            TrajCont.Pbessout.value(end + 1) = - TrajCont.Pgrid_avg.value(end);
            TrajCont.Ebess.value(end + 1) = sum(bat_energy_rating);     % Battery energy back to full
            TrajCont.Pevin.value(end + 1) = 0;

            TrajCont = UpdateTime(TrajCont, T_charge_ph3);
            
        % Utilization Ratio of Battery Energy
            TrajCont.uratio_energy.value(end + 1) = max_output_energy/sum(bat_energy_rating);
        % Error Check
            TrajCont.Time_errorcheck(end+1,:) = [T_discharge_ph1 + T_discharge_ph2 + T_charge_ph3, T_cycle, 3];

        else   % EV demand is small, bess is enough to charge EV
        % Phase 1: BESS and grid charge EV together
            T_discharge_ph1 = TrajENV.EV_demand(interval_ind)/TrajENV.Pcharge_level3;

            TrajCont.Pbessout.value(end + 1) = TrajENV.Pcharge_level3 - TrajCont.Pgrid_avg.value(end);
            TrajCont.Ebess.value(end + 1) = sum(bat_energy_rating) - T_discharge_ph1*(TrajENV.Pcharge_level3 - TrajCont.Pgrid_avg.value(end));
            TrajCont.Pevin.value(end + 1) = TrajENV.Pcharge_level3;
            
            TrajCont.Ebessout.value(end + 1) = T_discharge_ph1*(TrajENV.Pcharge_level3 - TrajCont.Pgrid_avg.value(end));
            TrajCont.TotalChargeTime.value(end + 1) = T_discharge_ph1;
            TrajCont.L3ChargeTime.value(end + 1) = T_discharge_ph1;
            TrajCont.CurtChargeTime.value(end + 1) = 0;

            TrajCont = UpdateTime(TrajCont, T_discharge_ph1);
        % Phase 2: grid charge BESS
            T_charge_ph2 = (sum(bat_energy_rating) - TrajCont.Ebess.value(end))/TrajCont.Pgrid_avg.value(end);  % compensate the energy back

            TrajCont.Pbessout.value(end + 1) = - TrajCont.Pgrid_avg.value(end);
            TrajCont.Ebess.value(end + 1) = sum(bat_energy_rating);
            TrajCont.Pevin.value(end + 1) = 0;

            TrajCont = UpdateTime(TrajCont, T_charge_ph2); 
        % Utilization Ratio of Battery Energy
            TrajCont.uratio_energy.value(end + 1) = T_discharge_ph1*(TrajENV.Pcharge_level3 - TrajCont.Pgrid_avg.value(end))/sum(bat_energy_rating);
        % Error Check
            TrajCont.Time_errorcheck(end+1,:) = [T_discharge_ph1 + T_charge_ph2, T_cycle, 2];
        end    

        % Waiting for the next EV to come
        TrajCont.Pbessout.value(end + 1) = 0;
        TrajCont.Ebess.value(end + 1) = TrajCont.Ebess.value(end);
        TrajCont.Pevin.value(end + 1) = 0;  
        TrajCont = UpdateTime(TrajCont, TrajENV.T_delay(interval_ind));
        % Update the averaged power from grid
        TrajCont.Pgrid_avg.value(end + 1) = 0;
        TrajCont.Pgrid_avg = UpdateTimeIndividual(TrajCont.Pgrid_avg, TrajENV.T_delay(interval_ind));
        % Update the Demand of EVs and Tdelay
        TrajCont.Evdemand.value(end + 1) = TrajENV.EV_demand(interval_ind);
        TrajCont.Tdelay.value(end + 1) = TrajENV.T_delay(interval_ind);
        
        % Update time
        TrajCont.uratio_energy = UpdateTimeIndividual(TrajCont.uratio_energy, T_cycle + TrajENV.T_delay(interval_ind));
        TrajCont.Evdemand = UpdateTimeIndividual(TrajCont.Evdemand, T_cycle + TrajENV.T_delay(interval_ind)); % Update the EV demand time 
        TrajCont.Tdelay = UpdateTimeIndividual(TrajCont.Tdelay, T_cycle + TrajENV.T_delay(interval_ind));   % Update the T_Delay time
        
        interval_ind = interval_ind + 1;
    end
end

function T_discharge = Discharge_time(EV_demand, Pgrid_fit, starting_time)
    dt = 0.1;
    timer = starting_time;
    Grid_energy = 0;
    while(Grid_energy <= EV_demand)
        Grid_energy = Pgrid_fit(timer)*dt + Grid_energy;
        timer = timer + dt;
    end
    T_discharge = timer - starting_time;
end

function TrajCont = UpdateTime(TrajCont, delta_time)
    TrajCont.Pbessout.delta_time(end + 1) = delta_time;
    TrajCont.Ebess.delta_time(end + 1) = delta_time;
    TrajCont.Pevin.delta_time(end + 1) = delta_time;

    TrajCont.Pbessout.abs_time(end + 1) = TrajCont.Pbessout.abs_time(end) + TrajCont.Pbessout.delta_time(end);
    TrajCont.Ebess.abs_time(end + 1) = TrajCont.Ebess.abs_time(end) + TrajCont.Ebess.delta_time(end);
    TrajCont.Pevin.abs_time(end + 1) = TrajCont.Pevin.abs_time(end) + TrajCont.Pevin.delta_time(end);
end

function x = UpdateTimeIndividual(x, delta_time)
    x.delta_time(end + 1) = delta_time;
    x.abs_time(end + 1) = x.abs_time(end) + delta_time;
end


function TrajCont = InitializeTraj()
    
    TrajCont.Time_errorcheck(1,1:3) = 0;
    
    TrajCont.Pgrid_avg.value = 0;
    TrajCont.Pbessout.value = 0;
    TrajCont.Ebess.value = 0;
    TrajCont.Pevin.value = 0;
    TrajCont.uratio_energy.value = 0;
    TrajCont.Evdemand.value = 0;
    TrajCont.Tdelay.value = 0;
    TrajCont.TotalChargeTime.value = 0;
    TrajCont.L3ChargeTime.value = 0;
    TrajCont.CurtChargeTime.value = 0;
    TrajCont.Ebessout.value = 0;

    TrajCont.Pgrid_avg.delta_time = 0;
    TrajCont.Pbessout.delta_time = 0;
    TrajCont.Ebess.delta_time = 0;
    TrajCont.Pevin.delta_time = 0;
    TrajCont.uratio_energy.delta_time = 0;
    TrajCont.Evdemand.delta_time = 0;
    TrajCont.Tdelay.delta_time = 0;
    
    TrajCont.Pgrid_avg.abs_time = 0;
    TrajCont.Pbessout.abs_time = 0;
    TrajCont.Ebess.abs_time = 0;
    TrajCont.Pevin.abs_time = 0;
    TrajCont.uratio_energy.abs_time = 0;
    TrajCont.Evdemand.abs_time = 0;
    TrajCont.Tdelay.abs_time = 0;
end