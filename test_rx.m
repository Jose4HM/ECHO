clc;
clear;

b = bladeRF('*:serial=8c14dec0ebc947b88424691d0dbf744f');

b.rx.channel = 'RX2';
b.rx.frequency = 5.8e9;
b.rx.samplerate = 1e6;
b.rx.bandwidth = 1e6;
b.rx.agc = 'manual';
b.rx.gain = -5;

b.rx.start();

[rx,~,count,ov] = b.receive(1024);

b.rx.stop();
 
fprintf("count = %d\n",count);
fprintf("overrun = %d\n",ov);

fprintf("min I = %.4f\n",min(real(rx)));
fprintf("max I = %.4f\n",max(real(rx)));

fprintf("min Q = %.4f\n",min(imag(rx)));
fprintf("max Q = %.4f\n",max(imag(rx)));

figure
plot(real(rx))

figure
plot(imag(rx))

figure
scatter(real(rx),imag(rx),5,'filled')
axis equal
grid on