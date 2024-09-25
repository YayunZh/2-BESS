function func_figure5(heter,Ipeak)

V_s_peak = 216;  
L1 = 0.1e-3;   
L2 = 0.1e-3;    
C = 5e-6;      
f = 60;          
omega = 2 * pi * f;  

Z_L1 = 1j * omega * L1;
Z_L2 = 1j * omega * L2;
Z_C = -1j / (omega * C);

phi_values = linspace(2*pi*0.001, 2*pi*0.999, 100);


figure;
hold on;

label = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];

for i = 1:10
   
        P_g_values = zeros(1, 100);
        Q_g_values = zeros(1, 100);
        I_s_peak = Ipeak(1,i);

    for k = 1:length(phi_values)
        phi = phi_values(k);
        
        
        V_s_rms = V_s_peak / sqrt(2); 
        I_s_rms = I_s_peak / sqrt(2); 
        V_s = V_s_rms * exp(1j * 0);  
        I_s = I_s_rms * exp(1j * phi); 

        V_C = V_s - I_s * Z_L1;

        I_C = V_C / Z_C;

        I_g = I_s - I_C;

        V_g = V_C - I_g * Z_L2;

        P_g_values(k) = real(V_g * conj(I_g));
        Q_g_values(k) = imag(V_g * conj(I_g));
    end

      scatter(Q_g_values, P_g_values,  5, 'filled','DisplayName', sprintf('%.1f A', label(i)));
end

set(gcf, 'Color', 'w');

percentages = [5, 10, 15, 20];  
title({'MAC-2BESS';['Battery Capacity Variation ', num2str(percentages(heter)), '%']});

xlabel('Q (VAR)');
ylabel('P (W)');

xlim([-6000, 6000]);
ylim([-6000, 6000]);
axis equal;
legend;
grid on;
hold off;

end
