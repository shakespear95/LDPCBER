clc;
clear all;
m=5;n=10;
R=(n-m)/n
frame=1;
Eb_N0=[0 0.5 1 1.5];

H=getH(m,n);
[G,valid]=H2G(H); 
while valid==0                  
H = getH(m,n);               
[G,valid]=H2G(H); 
end


  for i=1:length(Eb_N0)
    %sigma_2=1/(2*(10^(Eb_N0(i)/10))*R); 
    sigma=sqrt(1./(2*10^(Eb_N0(i)/10)*R));
    %N0 = 1/(exp(Eb_N0(i)*log(10)/10));
    %sigma=sqrt(N0/2)
    ber0(i)=0;
    %sigma_2
    %i
    for num=1:frame  
        num
        %x= rand(1,n)>0.5;
        %x=round(rand(1,m));
        x = (sign(randn(1,size(G,1)))+1)/2; % random bits
        y = mod(x*G,2);                     % coding 
        bpskmod = 2*y-1;                          %BPSK modulation
        h = 1/sqrt(2)*abs([randn(1,n) + j*randn(1,n)]); % Rayleigh channel
        z =h.*bpskmod + sigma*randn(size(bpskmod));%经Rayleigh信道的信息序列，也是准备译码的序列
        %z=bpskmod + sigma*randn(1,size(G,2));   % AWGN transmission
        f1=1./(1+exp(-2*z*0.8862/(sigma^2)));         % likelihoods
        f0=1-f1;
        [z_hat, success, k] = ldpc_decode(z,f0,f1,H);
        x_hat = z_hat(size(G,2)+1-size(G,1):size(G,2));
        x_hat = x_hat';
        err_max=find(x~=x_hat);             %寻找错误信息位
        num_eer=length(err_max)           %求出错误信息位位数
        ber(i)=num_eer/n;               %计算比特误码率BER
        ber0(i)=ber(i)+ber0(i);
        %ber(i);
        %ber0(i);
        %k;
        
    end %for num
    ber0(i)=ber0(i)/frame
    
end %for i
semilogy(Eb_N0,ber0,'b-o');
hold on

clc; 
clear all;

N = 10;
M= 5;
method = 1;

% Eliminate length-4 cycle
noCycle = 1;

% Number of 1s per column for LDPC matrix
onePerCol = 3;

% LDPC matrix reorder strategy (0 = First; 1 = Mincol; 2 = Minprod)
strategy = 2;

% EbN0 in dB
EbN0 = [0 0.5 1 1.5];

% Number of iteration;
iter = 1;

% Number of frame (N bits per frame)
frame = 1;

% Make the LDPC matrix
H = makeLdpc(M, N, 1, 1, onePerCol);

for i = 1:length(EbN0)
   
   ber1(i) = 0;
   ber2(i) = 0;
   ber3(i) = 0; 
   
   
   % Make random data (0/1)
   dSource = round(rand(M, frame));
   
   for j = 1:frame
      fprintf('Frame : %d\n', j);
      
      % Encoding message
      [c, newH] = makeParityChk(dSource(:, j), H, strategy);
      u = [c; dSource(:, j)];

      % BPSK modulation
      bpskMod = 2*u - 1;
      % Additional white gaussian noise
      N0 = 1/(exp(EbN0(i)*log(10)/10));
      tx = bpskMod + sqrt(N0/2)*randn(size(bpskMod));

      % Decoding (select decoding method)
      %vhat = decodeProbDomain(tx, newH,N0, iter);
      vhat1 = decodeLogDomain(tx, newH, N0, iter);
      vhat2 = decodeLogDomainSimple(tx, newH, iter);
      vhat3 = decodeBitFlip(tx, newH, iter);
   
      % Get bit error rate (for brevity, BER calculation includes parity bits)
      [num1, rat1] = biterr(vhat1', u);
      ber1(i) = (ber1(i) + rat1);
      
      [num2, rat2] = biterr(vhat2', u);
      ber2(i) = (ber2(i) + rat2);
      
       [num3, rat3] = biterr(vhat3', u);
      ber3(i) = (ber3(i) + rat3);
      
      
   end % for j
   
   % Get average of BER
   ber1(i) = ber1(i)/frame;
  
   ber2(i) = ber2(i)/frame;
   
   ber3(i) = ber3(i)/frame;
   
 
end % for i
 semilogy(EbN0, ber1, 'o-');
    hold on;
    semilogy(EbN0, ber2, 'b--');
    hold on;
    semilogy(EbN0, ber3, 'b-o')
xlabel('SNR(dB)');
ylabel('BER');
grid on;
hold off;

%semilogy(EbN0, ber4, 'r-*')
grid on;
hold off;