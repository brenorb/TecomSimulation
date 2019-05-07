%%%%%%%%%%%%%%%%%%%%%%%%
%Simulação 2 de TECOM
%Modulação SSB de sinais de voz
%1.1 Obtenha duas amostras de voz em formato .wav

[y1,Fs1,nbits1]=wavread('BRB');
[y2,Fs2,nbits2]=wavread('simulacao_de_tecom');

%igualando os tamanhos dos vetores
if length(y1)>length(y2)
    diferenca=length(y1)-length(y2);
    coluna=zeros(diferenca,1);
    y2=[y2;coluna];
else
    if length(y1)<length(y2)
        diferenca=abs(length(y1)-length(y2));
        coluna=zeros(diferenca,1);
        y1=[y1;coluna];
    end
end

L=length(y1);
t1=0:1/Fs1:(L-1)/Fs1;
t1=t1';
if Fs1~=Fs2
    t2=0:1/Fs2:(L-1)/Fs2;
    t2=t2';
end

%1.2 Faça a multiplexação dos dois sinais utilizando modulação SSB.
% Um dos sinais pode ser transmitido em banda base (sem modulação),
% enquanto o outro pode ser transmitido com modulação USB com 
% frequência de portadora de 4kHz. Explique os passos desta multiplexação.

%fazendo a transformada de Hilbert e pegando a parte USB
y1_mais=hilbert(y1)/2;
y1_menos=conj(y1_mais);
y1h=imag(y1_mais);
fc=4e3;

y1_USB=y1_mais.*exp(1i*2*pi*fc.*t1)+y1_menos.*exp(-1i*2*pi*fc.*t1);
mux_y1_USB_y2=y1_USB+y2;

%1.3 Mostre o espectro dos sinais originais, assim como o do sinal multiplexado.

NFFT = 2^nextpow2(L);
Ys1 = fft(y1,NFFT)/L;
f = Fs1/2*linspace(0,1,NFFT/2+1);

plot(f,2*abs(Ys1(1:NFFT/2+1)));
title('Espectro do sinal y1');
xlabel('Frequência (Hz)');
ylabel('Módulo da Amplitude');
print -djpeg -r300 Espectroy1

NFFT = 2^nextpow2(L);
Ys2 = fft(y2,NFFT)/L;
f = Fs2/2*linspace(0,1,NFFT/2+1);

 
plot(f,2*abs(Ys2(1:NFFT/2+1))) 
title('Espectro do sinal y2')
xlabel('Frequência (Hz)')
ylabel('Módulo da Amplitude')
print -djpeg -r300 Espectroy2

NFFT = 2^nextpow2(L);
Ysmux = fft(mux_y1_USB_y2,NFFT)/L;
f = Fs1/2*linspace(0,1,NFFT/2+1);

 
plot(f,2*abs(Ysmux(1:NFFT/2+1))) 
title('Espectro do sinal multiplexado')
xlabel('Frequência (Hz)')
ylabel('Módulo da Amplitude')
print -djpeg -r300 EspectroYmux

%1.4 Faça a demultiplexação do sinal, ou seja, recupere os sinais originais a partir do
% sinal multiplexado. Expliquem o processo. Como vocês avaliam a qualidade de
% voz?

ft = Fs1/2*linspace(0,1,NFFT);

Ys2demux=Ysmux;
for n=1:length(Ysmux)
    if (ft(n)>3.9e3)&&(ft(n)<(max(ft)-3.9e3))
        Ys2demux(n)=0;
    end
end
y2demux=ifft(Ys2demux,NFFT);
convy2=(max(y2)-min(y2))/(max(y2demux)-min(y2demux));
y2normdemux=y2demux*convy2;
wavwrite(y2normdemux,Fs2,nbits1,'simulacao_demux');

Ys1USBdemux=Ysmux-Ys2demux;
y1USBdemux=ifft(Ys1USBdemux,NFFT);
td=0:1/Fs1:(length(y1USBdemux)-1)/Fs1;
td=td';
y1demoddemux=y1USBdemux.*cos(2*pi*fc.*td);
Ys1demoddemux=fft(y1demoddemux,NFFT)/L;

W=Ys1demoddemux;
for n=1:length(Ys1demoddemux)
    if (ft(n)<3.9e3)||(ft(n)>(max(ft)-3.9e3))
        W(n)=0;
    end
end
y1demux=ifft(W,NFFT);
convy1=(max(y1)-min(y1))/(max(y1demux)-min(y1demux));
y1normdemux=y1demux*convy1;
wavwrite(y1normdemux,Fs2,nbits1,'BRB_demux');


%1.5 Mostre o espectro dos sinais demultiplexados.
% % 

NFFT = 2^nextpow2(L);
Ys1demux = fft(y1demux,NFFT)/L;
f = Fs1/2*linspace(0,1,NFFT/2+1);

 
plot(f,2*abs(Ys1demux(1:NFFT/2+1))) 
title('Espectro do sinal demultiplexado Y1')
xlabel('Frequência (Hz)')
ylabel('Módulo da Amplitude')
print -djpeg -r300 EspectroDemuxY1
NFFT = 2^nextpow2(L);
Ys2demux = fft(y2demux,NFFT)/L;
f = Fs1/2*linspace(0,1,NFFT/2+1);

 
plot(f,2*abs(Ys2demux(1:NFFT/2+1))) 
title('Espectro do sinal demultiplexado Y2')
xlabel('Frequência (Hz)')
ylabel('Módulo da Amplitude')
print -djpeg -r300 EspectroDemuxY2

% 2 Modulação de sinais de imagem
% 2.1 Importe uma imagem qualquer em formato bitmap. Mostre seu espectro,
% supondo que as amostras (pixels) são transmitidas a intervalos de 1ms.

img=imread('tsunami.bmp');
[x,y,z]=size(img);
imgcol=reshape(img,1,[],1);
Fsimg=1/1e-3;
Limg=x*y*z;

NFFT = 2^nextpow2(Limg);
Ysimg = fft(imgcol,NFFT)/Limg;
f = Fsimg/2*linspace(0,1,NFFT/2+1);

 
plot(f,2*abs(Ysimg(1:NFFT/2+1))) 
title('Espectro do sinal da imagem')
xlabel('Frequência (Hz)')
ylabel('Módulo da Amplitude')
print -djpeg -r300 Ymagem

% 2.2 Com base no espectro obtido, estime a largura de banda essencial B do
% sinal (99,9999% da energia).

Pot=abs(sum(Ysimg.^2));
Pot_aprox=0;
n=0;
while (Pot_aprox/Pot) < 0.999999
    n=n+1;
    Pot_aprox=Pot_aprox+Ysimg(n).^2;
end
B=f(n+1);

% 2.3 Faça a modulação do sinal em uma portadora de frequência 8B, usando 
% os seguintes esquemas. Descreva o que foi feito.
% 2.3.1 AM-SSB-USB
%ajustando o tamanho para a FFT
W=zeros(1,NFFT);
W(1:length(imgcol))=imgcol(1:length(imgcol));
imgcol=W;
Limgcol= length(imgcol);

timg=0:1/Fsimg:(NFFT-1)/Fsimg;

img_mais=hilbert(imgcol)/2; %Transformada de Hilbert
img_menos=conj(img_mais);   %
imgh=imag(img_mais);        %

%modulação
img_USB=img_mais.*exp(1i*2*pi*8*B.*timg)+img_menos.*exp(-1i*2*pi*8*B.*timg);

% 2.3.2 AM-VSB (use o filtro que achar mais adequado)

% Gerando espectro DSB
fc=8*B; %Frequencia da portadora

Fsm=2*(fc+Fsimg/2); 

Tsm=1/Fsm;

t=0:Tsm:(Limgcol)/Fsm-Tsm;
imgdsb=imgcol.*cos(2*pi*fc*t); %Modulando o sinal

Limgdsb=length(imgdsb);

NFFT = 2^nextpow2(Limgdsb); %N elementos FFT DSB
espectro_dsb = fftshift(fft(imgdsb,NFFT)/Limgdsb); 


%Espectro DSB com 0 deslocado para o meio do espectro para a visualização
% % da parte negativa
% %espectro_dsb=fft(sinal_dsb,NFFT)/tam_sinal_dsb;
% 
% freq_amostra=NFFT/Fs_modulado;  
% % Dado em amostra/Hz, nos da o inverso da frequencia por amostra. 
% % Para facilitar o trabalho com os gráficos, pois a funcao FFT nos da a
% % amplitude pelos indices, e não pelas frequencias.
% % Apenas multiplicando freq_amostra pelo indice, obtemos a frequencia de
% % tal. Lembrando de deslocar para o meio somando FFT/2
% f_dsb = Fs_modulado/2*linspace(-1,1,NFFT); %Frequencias da FFT DSB
% abs_espectro_dsb = abs(espectro_dsb);
% 
% figure(21);
% plot(f_dsb,abs_espectro_dsb);
% ylabel('Magnitude')
% xlabel('Frequency (Hz)')
% title('Espectro do sinal DSB')
% 
% 
% %%%%%% Gerando SSB-USB
% tam_espectro_dsb=length(espectro_dsb);
% n_amostras_filtro_usb=round(B99*freq_amostra);  %Número de amostras em que o filtro irá atuar.
% 
% ind_fc_pos=NFFT/2+1+round(fc*freq_amostra);     %Indice da freqencia da portadora, no primeiro quadrante
% ind_fc_neg=NFFT/2+1-round(fc*freq_amostra);     %Indice da frequencia da portadora, no segundo quadrante
% 
% ind_sig_pos=ind_fc_pos+round(B99*freq_amostra); %Indice da portadora até a largura de banda do sinal original
% ind_sig_neg=ind_fc_neg-round(B99*freq_amostra);
% 
% auxp=espectro_dsb(ind_fc_pos:ind_sig_pos);  %Valor auxiliar para definir os valores que o sinal USB irá assumir do DSB, guiado pelo número de amostras do filtro, indicando esta faixa
% auxn=espectro_dsb(ind_sig_neg:ind_fc_neg);
% espectro_usb=zeros(1,tam_espectro_dsb);         %Alocando memória para o espectro USB
% espectro_usb(ind_fc_pos:ind_sig_pos)=auxp;  %Agora o espectro é todo 0, menos aonde está a banda do sinal USB.
% espectro_usb(ind_sig_neg:ind_fc_neg)=auxn;
% 
% figure(22);
% abs_espectro_usb=abs(espectro_usb);
% plot(f_dsb,abs_espectro_usb);
% ylabel('Magnitude')
% xlabel('Frequency (Hz)')
% title('Espectro do sinal SSB-USB')
% 
% %%%%%%% Gerando SSB-VSB
% 
% filtro_vestigio_pos=linspace(0,0.5,n_amostras_filtro_usb/5);   
% %Filtro vestigial.Tomando a origem o começo do sinal USB, será uma reta
% %com inclinação de 45º começando em -B99/5 a 0. De 0 a B99 o filtro é 
% %ideal, ou seja, =1.
% aux2p=filtro_vestigio_pos.*espectro_dsb(ind_fc_pos-length(filtro_vestigio_pos)+1:ind_fc_pos); %Multiplicacao do filtro com um trecho da banda inferior vestigial
% filtro_vestigio_neg=linspace(0.5,0,n_amostras_filtro_usb/5);    %Filtro vestigial para a parte negativa.
% aux2n=filtro_vestigio_neg.*espectro_dsb(ind_fc_neg:ind_fc_neg+length(filtro_vestigio_neg)-1);
% espectro_vsb=espectro_usb; %Alocando memória para o espectro vsb
% espectro_vsb(ind_fc_pos-length(filtro_vestigio_pos)+1:ind_fc_pos)=aux2p;   
% espectro_vsb(ind_fc_neg:ind_fc_neg+length(filtro_vestigio_neg)-1)=aux2n; 
% abs_espectro_vsb=abs(espectro_vsb);
% 
% 
% 
% figure(23);
% plot(f_dsb,abs_espectro_vsb);
% ylabel('Magnitude')
% xlabel('Frequency (Hz)')
% title('Espectro do sinal SSB-VSB')
% 
% % 2.5 %%%%%%%%%%%%%%%%
% 
% %%%%% Demodulando o sinal USB
% 
% sinal_usb=ifft(espectro_usb,NFFT);
% 
% sinal_usb_x_portadora=sinal_usb.*cos(2*pi*fc*t);
% NFFT = 2^nextpow2(length(sinal_usb_x_portadora));
% espectro_sinal_usb_x_portadora = fft(sinal_usb_x_portadora,NFFT);
% 
% figure(24)
% plot(f_dsb,abs(espectro_sinal_usb_x_portadora));
% 
% espectro_sinal_usb_demodulado=zeros(1,length(t));        %Alocando memoria para o sinal demodulado
% espectro_imagem_usb=espectro_sinal_usb_x_portadora(NFFT/2+1:NFFT/2+1+round(B99*freq_amostra));
% f_img=freq_amostra*length(espectro_imagem_usb)/Fs*linspace(0,1,length(espectro_imagem_usb)); 
% 
% figure(25);
% plot(f_img,abs(espectro_imagem_usb));
% ylabel('Magnitude')
% xlabel('Frequency (Hz)')
% title('Espectro do sinal demodulado - USB')
% 
% ganho_usb=1000000*1.1;
% imagem_usb=ifft(espectro_imagem_usb,NFFT)*ganho_usb;
% 
% figure(26);
% plot(1:length(img_line),abs(img_line));
% 
% figure(27);
% plot(1:length(imagem_usb),abs(imagem_usb));
% 
% imagem_usb=uint8(imagem_usb);
% imagem_usb_parc=imagem_usb(1:length(img_line_parc));
% imagem=reshape(imagem_usb_parc,256,256,3);
% 
% imwrite(imagem,'imagem_usb.bmp','bmp');
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %%%%%%%%Demodulando VSB
% sinal_vsb=ifft(espectro_vsb,NFFT);
% 
% sinal_vsb_x_portadora=sinal_vsb.*cos(2*pi*fc*t);
% 
% NFFT = 2^nextpow2(length(sinal_vsb_x_portadora));
% espectro_sinal_vsb_x_portadora = fft(sinal_vsb_x_portadora,NFFT);
% 
% figure(30)
% plot(f_dsb,abs(espectro_sinal_vsb_x_portadora));
% 
% espectro_sinal_vsb_demodulado=zeros(1,length(t));        %Alocando memoria para o sinal demodulado
% espectro_imagem_vsb=espectro_sinal_vsb_x_portadora(NFFT/2+1:NFFT/2+1+round(B99*freq_amostra));
% f_img=freq_amostra*length(espectro_imagem_vsb)/Fs*linspace(0,1,length(espectro_imagem_vsb));
% 
% filtro_vestigio_demod=linspace(2/3,1,length(espectro_imagem_vsb)/5);
% espectro_imagem_vsb(1:round(length(espectro_imagem_vsb)/5))=filtro_vestigio_demod.*espectro_imagem_vsb(1:round(length(espectro_imagem_vsb)/5));
% 
% figure(31);
% plot(f_img,abs(espectro_imagem_vsb));
% ylabel('Magnitude')
% xlabel('Frequency (Hz)')
% title('Espectro do sinal demodulado - VSB')
% 
% ganho_vsb=1000000*2;
% imagem_vsb=ifft(espectro_imagem_vsb,NFFT)*ganho_vsb;
% 
% figure(32);
% plot(1:length(img_line),img_line);
% 
% figure(33);
% plot(1:length(imagem_vsb),abs(imagem_vsb));
% 
% imagem_vsb=uint8(imagem_vsb);
% imagem_vsb_parc=imagem_vsb(1:length(img_line_parc));
% imagem=reshape(imagem_vsb_parc,256,256,3);
% 
% imwrite(imagem,'imagem_vsb.bmp','bmp');
% 
% 

% 2.4 Mostre o espectro dos sinais modulados.

NFFT = 2^nextpow2(length(img_USB));
Ysimg_USB = fft(img_USB,NFFT)/length(img_USB);
f = Fsimg/2*linspace(0,1,NFFT/2+1);

 
plot(f,2*abs(Ysimg_USB(1:NFFT/2+1))) 
title('Espectro do sinal da imagem modulado USB')
xlabel('Frequência (Hz)')
ylabel('Módulo da Amplitude')
print -djpeg -r300 ImgUSB

% 2.5 Faça a demodulação destes sinais e compare a imagem demodulada nos
% dois casos com a imagem original. Mostre-as no relatório.

%do sinal SSB

% 2.6 Mostre também o espectro das mensagens demoduladas nos dois casos.