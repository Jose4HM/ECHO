function [z_out, f_off, vueltas, deriva_pp] = quitar_offset_portadora(z, Fs, fc_track)
% La correccion por defecto elimina solo el offset lineal TX/RX. Un
% seguimiento paso-bajo de fase a 0.06 Hz tenia una constante temporal de
% varios segundos y generaba una joroba artificial en los bordes, que el
% detector confundia con respiracion. Se conserva fc_track solo para
% diagnosticos offline explicitos.
if nargin < 3 || isempty(fc_track), fc_track = 0; end

N = numel(z);
n = (0:N-1).';

% Pico grueso cerca de DC.
Nf = 2^nextpow2(8*N);
Zf = fft(z, Nf);
f = (0:Nf-1).'/Nf*Fs;
f(f > Fs/2) = f(f > Fs/2)-Fs;
mag = abs(Zf);
mag(abs(f) > 2) = 0;
[~, k] = max(mag);
f0 = f(k);

z1 = z .* exp(-1j*2*pi*f0*n/Fs);

% Refinamiento mediante pendiente lineal de fase.
psi = unwrap(angle(z1));
u = linspace(-1,1,N).';
p = polyfit(u, psi, 1);
% du/dn = 2/(N-1). Convertir pendiente por u a rad/muestra.
rad_per_sample = p(1) * 2/max(N-1,1);
f1 = rad_per_sample*Fs/(2*pi);

z2 = z1 .* exp(-1j*2*pi*f1*n/Fs);
f_off = f0 + f1;
vueltas = f_off*N/Fs;

if fc_track <= 0
    z_out = z2;
    deriva_pp = 0;
    return;
end

psi = unwrap(angle(z2));
Wn = fc_track/(Fs/2);
[bb, aa] = butter(2, Wn, 'low');
psi_lento = filtfilt(bb, aa, psi);
z_out = z2 .* exp(-1j*psi_lento);
deriva_pp = max(psi_lento)-min(psi_lento);
end
