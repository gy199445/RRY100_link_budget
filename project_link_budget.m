% GEO stationary satellite link budget calculation
%constants
Re = 6371E3; %earth r
k = 1.38064852E-32; %boltzman's constant
h = 35786E3; %geostationary orbit height
% center frequency
f = 26.644E9; %26.644 GHz
% speed of light
c = 3E8;
% antennas
groundAntennaGain = 0.6*(pi*5*f/c)^2; %5-meter dish
groundAntennaGain = 10*log10(groundAntennaGain);
satelliteAntennaGain = 10*log10(0.7*(pi*3.5*f/c)^2); %assume 3.5m dish
PtxSat = 30; %1000W
PtxGround = 10*log10(5000); %5000W
%amplifiers
%locations
%*********Plug in your settings here************
%sub satellite points (in deg.)
satLongitude = 45;
satLatitude = 0;
%target ground station location
targetLongitude = -56;
targetLatitude = -67;
%ground station 0.001% rain rate
targetRainRate = 250; %zone P
%polarization angle
polarizationAngle = pi/2; %assume vertical polarization
%coding gain
codingGain = 5;
%*********Plug in your settings here***********
%free space loss
%link distance in meters
DFs = ...
    sqrt(1+0.42*(1-cos(targetLatitude*pi/180)*cos(abs(targetLongitude-satLongitude)*pi/180)))*h;
LFs = 10*log10((4*pi*DFs*f/c)^2);
%rain atennuation
k_H = 0.0751; k_v = 0.0691; alpha_H = 1.099;alpha_v = 1.065;
[elevationAngle,~] = ...
    calcElevationAngle([satLatitude;satLongitude],[targetLatitude;targetLongitude],Re,h);
rainAttenuation = calcRainAtt([k_H;k_v],[alpha_H;alpha_v],targetRainRate,elevationAngle,polarizationAngle,DFs);

%display result:
fprintf('Free space loss: %.4f dB\n Rain attenuation %4f dB\n',LFs,rainAttenuation)
%functions to be used
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
function [rainAttenuation] = calcRainAtt...
    (k,alpha,RR,elevAngle,polarAngle,linkDistance)
    %compute rain attenuation
    %follow the calculation of 'introduction to RF propagation' page 243
    %k = [k_H;k_v]
    %alpha = [alpha_H;alpha_v]
    %RR = 99.999% rain rate (mm/h)
    %elevAngle = ground station elevation angle in rad
    %polarAngle = polarization angle in rad
    %sate = [satellite_latitude; satellite_longitude]
    %target = [ground_latitude; ground_longtitude]
    %linkDistance = link distance in meter
    linkDistance = linkDistance/1000;%required
    K = (k(1)+k(2)+(k(1)-k(2))*cos(elevAngle)^2*cos(2*polarAngle))/2;
    ALPHA = (k'*alpha + (k(1)*alpha(1) -...
        k(2)*alpha(2))*cos(elevAngle)^2 * cos(2*polarAngle))/(2*K);
    d0 = 35*exp(-0.015*RR);
    r = 1/(1+linkDistance/d0);
    %the rain attenuation at 99.99% (0.01)
    rainAttenuation = K*(RR^ALPHA)*linkDistance*r;
end