%% Compressor P-Q Map
load cmp_PQ.mat

figure('Name','Compressor P-Q Map')
hold on
grid on
box on

for i = 1:size(N_table,2)
    plot(Wcr(i,:),pr);
end

plot(out.W_cmp.Data,out.pout_cmp.Data./out.pin_cmp.Data,'r');
plot(out.W_cmp.Data(end),out.pout_cmp.Data(end)./out.pin_cmp.Data(end),'rx');

xlim([0 0.1]);
xlabel("Corrected mass flow rate [kg/s]");
ylabel("Pressure ratio");