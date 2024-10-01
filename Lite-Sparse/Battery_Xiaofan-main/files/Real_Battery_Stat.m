function RealBatStat = Real_Battery_Stat(Bat_info_tab)

[~,ind_1] = sort(Bat_info_tab(:,1));
Bat_power_info_tab = Bat_info_tab(ind_1,:);
[~,ind_2] = sort(Bat_info_tab(:,2));
Bat_energy_info_tab = Bat_info_tab(ind_2,:);
if(~isequal(Bat_power_info_tab, Bat_energy_info_tab))
    fprintf('Wrong Assumptions');
end
Bat_pdf_info_tab = [0.5*Bat_power_info_tab(1,1), 0.5*Bat_power_info_tab(1,2), 0;Bat_power_info_tab];
Bat_cdf_info_tab = [];
for i = 1:size(Bat_pdf_info_tab,1)
    Bat_cdf_info_tab(i,:) = Bat_pdf_info_tab(i,:);
    Bat_cdf_info_tab(i,3) = sum(Bat_pdf_info_tab(1:i,3));
end
cdf_power_func = fit(Bat_cdf_info_tab(:,1), Bat_cdf_info_tab(:,3), 'linearinterp');
inv_cdf_power_func = fit(Bat_cdf_info_tab(:,3), Bat_cdf_info_tab(:,1), 'linearinterp');
cdf_energy_func = fit(Bat_cdf_info_tab(:,2), Bat_cdf_info_tab(:,3), 'linearinterp');
inv_cdf_energy_func = fit(Bat_cdf_info_tab(:,3), Bat_cdf_info_tab(:,2), 'linearinterp');

RealBatStat.real_value.power_dist.cdf = cdf_power_func;
RealBatStat.real_value.power_dist.invcdf = inv_cdf_power_func;
RealBatStat.real_value.energy_dist.cdf = cdf_energy_func;
RealBatStat.real_value.energy_dist.invcdf = inv_cdf_energy_func;

Bat_cdf_info_tab_power_base = (Bat_cdf_info_tab(1,1) + Bat_cdf_info_tab(end,1))/2;   % base power value for normalization
Bat_cdf_info_tab_energy_base = (Bat_cdf_info_tab(1,2) + Bat_cdf_info_tab(end,2))/2;  % base energyvalue for normalization
cdf_power_nom_func = fit(Bat_cdf_info_tab(:,1)/Bat_cdf_info_tab_power_base, Bat_cdf_info_tab(:,3), 'linearinterp');
inv_cdf_power_nom_func = fit(Bat_cdf_info_tab(:,3), Bat_cdf_info_tab(:,1)/Bat_cdf_info_tab_power_base, 'linearinterp');
cdf_energy_nom_func = fit(Bat_cdf_info_tab(:,2)/Bat_cdf_info_tab_energy_base, Bat_cdf_info_tab(:,3), 'linearinterp');
inv_cdf_energy_nom_func = fit(Bat_cdf_info_tab(:,3), Bat_cdf_info_tab(:,2)/Bat_cdf_info_tab_energy_base, 'linearinterp');

RealBatStat.norm_value.power_dist.cdf = cdf_power_nom_func;
RealBatStat.norm_value.power_dist.invcdf = inv_cdf_power_nom_func;
RealBatStat.norm_value.energy_dist.cdf = cdf_energy_nom_func;
RealBatStat.norm_value.energy_dist.invcdf = inv_cdf_energy_nom_func;
% for function testing
% power_sample = inv_cdf_power_func(rand(1,50000));
% histogram(power_sample,100);
% energy_sample = inv_cdf_energy_func(rand(1,50000));
% histogram(energy_sample,100);

end
