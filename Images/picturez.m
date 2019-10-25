y = [900;1060;1200;1080;800;560;580;1000;900;1060;1200;1080;800;560;580;1000;900;1060;1200;1080;800;560;580;1000];
figure(1)
plot(y,'-b','LineWidth',3)
grid on
% axis([0 24 24 28]);
title 'Hourly Load Curve of Grid Station'
xlabel 'Time (Hour of Day)'
ylabel 'Hourly Peak Demand (MW)'
legend('(MW)','Location','best');

    ysolar = 200*[0.0008
        0.0004
        0.0053
        0.0035
        0.0032
        0.0345
        0.0909
        0.1314
        0.2034
        0.3605
        0.6993
        0.8269
        0.8168
        0.8069
        0.7280
        0.4920
        0.2317
        0.1155
        0.0841
        0.0672
        0.0062
        0.0024
        0.0081
        0.0097
        ];

    figure(2)
plot(ysolar,'-b','LineWidth',3)
grid on
axis([0 24 0 200]);
title 'Power Curve of Solar'
xlabel 'Time (Hour of Day)'
ylabel 'Hourly Power (MW)'
legend('(MW)','Location','best');
    ybattery = 50*[0.6692
        0.1904
        0.3689
        0.4607
        0.9816
        0.1564
        0.8555
        0.6448
        0.3763
        0.1909
        0.4283
        0.4820
        0.1206
        0.5895
        0.2262
        0.3846
        0.5830
        0.2518
        0.2904
        0.6171
        0.2653
        0.8244
        0.9827
        0.7302
        ];
    
    figure(3)
plot(ybattery,'-b','LineWidth',3)
grid on
title 'Battery Storage'
xlabel 'Time (Hour of Day)'
ylabel 'Hourly Power (MW)'
legend('(MW)','Location','best');
axis([0 24 0 50]);
    ywind = 100*[0.3395
        0.9516
        0.9203
        0.0527
        0.7379
        0.2691
        0.4228
        0.5479
        0.9427
        0.4177
        0.9831
        0.3015
        0.7011
        0.6663
        0.5391
        0.6981
        0.6665
        0.1781
        0.1280
        0.9991
        0.1711
        0.0326
        0.5612
        0.8819
        ];
figure(4)
plot(ywind,'-b','LineWidth',3)
grid on
axis([0 24 0 100]);
title 'Power Curve of Wind'
xlabel 'Time (Hour of Day)'
ylabel 'Hourly Power (MW)'
legend('(MW)','Location','best');
y1 = y-ysolar-ybattery-ywind;
figure(5)
plot(y1,'-b','LineWidth',3)
grid on
axis([0 24 0 1500]);
title 'Load Curve for Generators'
xlabel 'Time (Hour of Day)'
ylabel 'Hourly Power (MW)'
legend('(MW)','Location','best');
