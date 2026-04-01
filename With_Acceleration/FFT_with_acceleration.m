%%Using a Fast Fourier Transform Algorithm FFT to compute a Discrete Fourier
%%Transform DFT
data=readmatrix("rpm_ramp_2sec_to_2M_with_noise.csv","NumHeaderLines",1);

N=length(data); %length(data); %length of signal/number of samples
f_s = 200000; %sampling frequency of sensor in Hz, I think ours is smaller than that right?
T = 1/f_s; %sampling period of sensor

f_c = 1000; %sampling frequence for RPM calculation
c_s = f_s/f_c; %window(rows of data) over which each value of RPM is calculated

%original signal plot
Xn = data(:,2);
plot(Xn, '.-');
title("Signal");
xlabel("t");
ylabel("Xn(t)");

%One-sided, strength of each frequency iterated through every RPM
%calculation sampling window
count=0; % for indexing the two arrays below
RPM_list=[];
time_list=[];

for i = 1:c_s:(N-c_s) %step by the sampling window for RPM calculation
    segment = Xn(i:i+c_s-1);

    Y=fft(segment); % do the fft
    count=count+1;
    F2 = abs(Y/c_s);
    F1 = F2(1:(c_s/2+1)); %y
    f_axis = f_s*(0:(c_s/2))/c_s; %x
    %F1 not scaled yet so cannot use it to determine Power but we don't
    %need that for our purpose
    [~, idx] = max(F1); % column at which max peak magnitude occurs
    peak_freq = f_axis(idx); % find the frequency corresponding to the peak magnitude
    peak_rpm = peak_freq*60/40; % number of turbine blades is 40
    RPM_list(count) = peak_rpm;
    time_list(count) = (i + (c_s/2) - 1)/f_s;
    % fprintf("peak RPM %f\n", peak_rpm); % Convert frequency to RPM
end

%plot of RPM Signal over time
figure
plot(time_list, RPM_list)
xlabel('Time (s)')
ylabel('RPM')
title('RPM Signal over time')
axis tight