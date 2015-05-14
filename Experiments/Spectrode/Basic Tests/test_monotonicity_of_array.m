%Test monotonicity of array
%% set parameters of geometric model
a = 0.5;
b = 1.5;
K = 10;
[t,w]  = geometric_model(a,b,K);
gamma = 1/6;

%% test ODE method
epsi_array = 10.^(-(1:3))';
gamma_array = 2.^(linspace(-2,2,4))';
K_err = zeros(length(epsi_array),length(gamma_array));
l_err = inf(length(epsi_array),length(gamma_array));
u_err = inf(length(epsi_array),length(gamma_array));

for i=1:length(epsi_array)
    epsi = epsi_array(i);
    for j=1:length(gamma_array)
        gamma = gamma_array(j);
        f = compute_esd_ode(t, w, gamma,epsi);
        figure, plot(f, f-sort(f))
        str = sprintf('epsi= %d; gamma = %d;', epsi, gamma);
        title(str);
    end
end