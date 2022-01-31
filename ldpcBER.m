% Bit error rate of BPSK modulated LDPC codes under AWGN channel
%
%
% Copyright Bagawan S. Nugroho, 2007 
% http://bsnugroho.googlepages.com

clc; 
clear all;

% LDPC matrix size, rate must be 1/2
% Warning: encoding - decoding can be very long for large LDPC matrix!
M = 256;
N = 512;

% Method for creating LDPC matrix (0 = Evencol; 1 = Evenboth)
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
iter = 10;

% Number of frame (N bits per frame)
frame = 10;

% Make the LDPC matrix
H = makeLdpc(M, N, 1, 1, onePerCol);

for i = 1:length(EbN0)
   
   ber1(i) = 0;
   ber2(i) = 0;
   ber3(i) = 0;
   ber4(i) = 0;
   
   
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
      vhat1 = decodeProbDomain(tx, newH,N0, iter);
      vhat2 = decodeLogDomain(tx, newH, N0, iter);
      vhat3 = decodeLogDomainSimple(tx, newH, iter);
      vhat4 = decodeBitFlip(tx, newH, iter);
   
      % Get bit error rate (for brevity, BER calculation includes parity bits)
      [num1, rat1] = biterr(vhat1', u);
      ber1(i) = (ber1(i) + rat1);
      
      [num2, rat2] = biterr(vhat2', u);
      ber2(i) = (ber2(i) + rat2);
      
       [num3, rat3] = biterr(vhat3', u);
      ber3(i) = (ber3(i) + rat3);
      
      [num4, rat4] = biterr(vhat4', u);
      ber4(i) = (ber4(i) + rat4);
      
   end % for j
   
    %Get average of BER
   ber1(i) = ber1(i)/frame;
   ber2(i) = ber2(i)/frame;
   ber3(i) = ber3(i)/frame;
   ber4(i) = ber4(i)/frame;
   
end % for i

% Plot the result
semilogy(EbN0, ber1, 'o-');
hold on;
semilogy(EbN0, ber2, 'g--');
hold on;
semilogy(EbN0, ber3, 'b-o')
hold on;
semilogy(EbN0, ber4, 'r-*')

xlabel('SNR(dB)');
ylabel('BER');
grid on;
grid on;
hold off;
