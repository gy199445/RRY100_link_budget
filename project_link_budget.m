% GEO stationary satellite link budget calculation
%constants
Re = 6371E3; %earth r
k = 1.38064852E-23; %boltzman's constant
h = 35786E3; %geostationary orbit height
% center frequency
%f = 26.644E9; %26.644 GHz
f=8E9;
% speed of light
c = 3E8;
% antennas
groundAntennaGain = 0.6*(pi*3.5*f/c)^2; %5-meter dish
groundAntennaGain = 10*log10(groundAntennaGain);
groundAntennaHeight = 5/1000;
satelliteAntennaGain = 10*log10(0.7*(pi*0.5*f/c)^2); %assume 1.5m dish
PtxSat = 10*log10(300);
PtxGround = 10*log10(1000); %1000W
%tx rate
Rb = 30*10^6;%10Mbps
%amplifiers
%locations
%*********Plug in your settings here************
%sub satellite points (in deg.)
satLongitude = -30;
satLatitude = 0;
%target ground station location
targetLongitude = -82;
targetLatitude = -54;
%ground station 0.01% rain rate
targetRainRate = 145;
%polarization angle
polarizationAngle = pi/4; %assume circular polarization
%coding gain
codingGain = 5;
%code rate
codeRate = 0.5;
%roll off factor of pulse shaping
rollOffFactor = 0.5;
%calculate bandwidth
infoBitRate = 10*10^6; %in bps
rawBitRate = infoBitRate/codeRate;
MCS = 4;%for QPSK (2 bits/symbol)
BW = rawBitRate/log2(MCS) * (1+rollOffFactor);
%*********Plug in your settings here***********
%free space loss
%link distance in meters
DFs = ...
    sqrt(1+0.42*(1-cos(targetLatitude*pi/180)*cos(abs(targetLongitude-satLongitude)*pi/180)))*h;
LFs = 10*log10((4*pi*DFs*f/c)^2);
%rain atennuation
%k_H =  0.1724; k_v = 0.1669 ; alpha_H = 0.9884;alpha_v = 0.9421; %(26GHz)
%k_H =  0.0001321; k_v =  0.0001464; alpha_H = 1.1209;alpha_v = 1.0085; %(2.5GHz)
k_H = 0.004115;alpha_H = 1.3905;k_v = 0.003450;alpha_v = 1.3797;%(8ghz)
[elevationAngle,~] = ...
    calcElevationAngle([satLatitude;satLongitude],[targetLatitude;targetLongitude],Re,h);
rainAttenuation = ...
    calcRainAtt(targetRainRate,[k_H;k_v],[alpha_H;alpha_v],groundAntennaHeight,elevationAngle,polarizationAngle,targetLatitude*pi/180,f/(10^9));
%% noise temperature at receivers
%DOWNLINK
%antenna temperature
rainAttenuationLinear = 10^(rainAttenuation/10);
TRain = 275*(1-1/rainAttenuationLinear);
Tsky = 9/(rainAttenuationLinear);
TpFeeder = 290;%ground station antenna feeder phy. temp.
LFeeder = 1.1220; %loss = 0.5dB
TFeederEff = TpFeeder*(LFeeder-1);
noiseTempAMP = noiseTempCalc([1 2 2],[15 30 40]);
noiseTemp = (TRain + Tsky + TFeederEff)/LFeeder + noiseTempAMP;
noiseTempUplink = (290+LFeeder)/LFeeder + noiseTempAMP;
N0 = noiseTemp*k;
N0Uplink = noiseTempUplink*k;
N = N0 * BW;
NUplink = N0Uplink * BW;
%% carrier power
%DOWNLINK
PcDownLink = sum([satelliteAntennaGain PtxSat -LFs -rainAttenuation groundAntennaGain -0.5 + codingGain]);
%UPLINK
PcUpLink = sum([satelliteAntennaGain PtxGround -LFs -rainAttenuation groundAntennaGain -0.5 +codingGain]);
%% C/N (Es/N0)
C_N = PcDownLink - 10*log10(N);
C_N_Uplink = PcUpLink - 10*log10(N);
%Eb/N0
EbN0 = C_N  + 10*log10(infoBitRate/BW);
EbN0_Uplink = C_N_Uplink + 10*log10(infoBitRate/BW);
%% save result to xls:
xlsName = 'link_budget';
saveData = {'DOWNLINK',' ','UPLINK',' ';
    'Sat. EIRP(dBW)',satelliteAntennaGain+PtxSat,'Grd. EIRP(dBW)',groundAntennaGain+PtxGround;
    'Link distance(km)',DFs/1000,'Link distance(km)',DFs/1000;
    'Free space loss(dB)',LFs,'Free space loss(dB)',LFs;
    'Rain rate(mm/h)',targetRainRate,'Rain rate(mm/h)',targetRainRate;
    'Rain att.(dB)',rainAttenuation,'Rain att.(dB)',rainAttenuation;
    'Rx antenna gain(dBi)',groundAntennaGain,'Rx antenna gain(dBi)',satelliteAntennaGain;
    'Rx feeder loss(dB)',0.5,'Rx feeder loss(dB)',0.5;
    'Elevation(deg)',elevationAngle*180/pi,'Elevation(deg)',elevationAngle*180/pi;
    'Carrier Power(dBW)',PcDownLink,'Carrier Power(dBW)',PcUpLink;
    'Eff. noise temp.(K)',noiseTemp,'Eff. noise temp.(K)',noiseTempUplink;
    'N0(W/Hz)',N0,'N0(W/Hz)',N0Uplink;
    'Coding gain(dB)',5,'Coding gain(dB)',5;
    'BandWidth(MHz)',BW/(10^6),'BandWidth(MHz)',BW/(10^6);
    'Eb/N0(dB)',EbN0,'Eb/N0(dB)',EbN0_Uplink;
    'BER(%)',berawgn(EbN0,'qam',4),'BER(%)',berawgn(EbN0_Uplink,'qam',4)};
xlswrite('london',saveData)
function [elevation,azimuth] = calcElevationAngle(sate,target,Re,h)
    %compute the antenna elevation angle of ground station
    %sate = [satellite_latitude; satellite_longitude] in degree
    %target = [ground_latitude; ground_longtitude] in degree
    %elevation = elevation angle in rad
    %Re = earth's radius
    %h = satellite height
    %lecture 04 page 21
    %central angle
    sate = sate*pi/180;
    target = target*pi/180;
    delta = acos( sin(sate(1)) * sin(target(1)) ...
        + cos(sate(1)) * cos(target(1)) * cos(sate(2) - target(2)));
    %elevation and azimuth
    elevation = atan((cos(delta)-(Re/(Re+h))) / sin(delta));
    azimuth = sign(sin(sate(2) - target(2))) ...
        * acos( (sin(sate(1)) - sin(target(1))*cos(delta)) / (sin(delta)*cos(target(1))));
end
function [Ap] = calcRainAtt...
    (RR,k,alpha,hs,elevAngle,polarAngle,targetLatitude,f)
    %compute rain attenuation
    %follow ITU_P618 page 7
    %RR 0.01% rain rate mm/h
    %k [k_h k_v]
    %alpha [alpha_h, alpha_v] %from ITU-R P.838-3
    %hs ground station height above mean sea level (km)
    %elevAngle = ground station elevation angle in rad
    %polarAngle in rad
    %targetLatitude target latitude in rad
    %f: frequency in GHZ
    %Re effective radius of earth 8500km
    %rain height (ITU-R P.839)
    Re = 8500;
    hr = calcRainHeight(targetLatitude) + 0.36;
    elevAngle_deg = elevAngle*180/pi;
    %slant-path
    if(elevAngle_deg>=5)
        Ls = (hr-hs)/sin(elevAngle);
    else
        Ls = 2*(hr-hs)/(sqrt(sin(elevAngle)^2+(2*(hr-hs))/Re)+sin(elevAngle));
    end
    %horizontal projection
    Lg = Ls*cos(elevAngle);
    %specific attenuation
    %k and gamma_r is from ITU-R P838
    K = (k(1)+k(2)+(k(1)-k(2))*cos(elevAngle)^2*cos(2*polarAngle))/2;
    ALPHA = (k(1)*alpha(1)+k(2)*alpha(2)+...
        (k(1)*alpha(1)-k(2)*alpha(2))*cos(elevAngle)^2*...
    cos(2*polarAngle))/2/K;
    gamma_r = K*(RR)^ALPHA;
    %horizontal reduction factor
    r_001 = 1/(1+0.78*sqrt(Lg*gamma_r/f)-0.38*(1-exp(-2*Lg)));
    %vertical adjustment factor
    zeta = atan((hr-hs)/(Lg*r_001));%in rad
    if zeta>elevAngle
        LR = Lg*r_001/cos(elevAngle);
    else
        LR = (hr-hs)/sin(elevAngle);
    end
    if abs(targetLatitude)<36*pi/180
        chi = 36*pi/180 - abs(targetLatitude);%rad
    else
        chi = 0;
    end
    chi_deg = chi*180/pi;
    v_001 = 1/(1 + sqrt(sin(elevAngle)) * (31 * (1-exp(-(elevAngle_deg/(1+chi_deg)))) * sqrt(LR*gamma_r) / (f^2) -0.45 ));
    LE = LR*v_001;
    A_001 = gamma_r*LE;
    %scale to 0.001
    targetLatitude_deg = targetLatitude * 180/pi;
    if abs(targetLatitude_deg)<36 && elevAngle_deg>=25
        beta = -0.005*(abs(targetLatitude_deg) - 36);
    else
        beta = -0.005*(abs(targetLatitude_deg) - 36)+1.8-4.25*sin(elevAngle);
    end
    p = 0.001;
    Ap = A_001*(p/0.01)^(- (0.655+ 0.033*log(p) - 0.045*log(A_001)-beta*(1-p)*sin(elevAngle)) );
end
function [noiseTemp] = noiseTempCalc(noiseFigures,gains)
    %compute the effective inuput noise temperature of
    %3 cascaded amplifiers
    %noiseFigure [1*3] noise figure in dB
    %gains [1*3] AMP gains in dB
    %noiseTemp effective input noise temp in K
    Tref = 290;
    gains = 10.^(gains/10);%convert to linear
    EffInputT = Tref*(10.^(noiseFigures/10) - 1);%in K
    noiseTemp = EffInputT(1) + ...
    EffInputT(2)/(gains(1)) + ...
        EffInputT(3)/(gains(1)*gains(2));
end
function [h0] = calcRainHeight(targetLatitude)
%h0   height of isotherm
%targetLatitude: latitude in deg
    if targetLatitude>23
        h0 = 5-0.075*(targetLatitude-23);
    elseif targetLatitude>0 && targetLatitude<23
        h0 = 5;
    elseif targetLatitude<=0 && targetLatitude>=-21
        h0 = 5;
    elseif targetLatitude<-21 && targetLatitude>=-71
        h0=5+0.1*(targetLatitude+21);
    else
        h0=0;
    end
end