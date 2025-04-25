function [s, ds, dds] = generatePulse(md, tau_0, tau, normalize)
%function s = generatePulse(md, tau_0, tau, normalize)
%
% Generates a pulse (RC, RRC) specified by md at position tau_0.
%The pulse is defined over the vector tau;
% Also generates the analytical derivatives (First and Second order).
%normalize is optional:
% 1: normalize to unit maximum amplitude
% 2: normalize to unit energy (discrete time)
% 3: normalize to unit energy (continuous time, acc. to tau vector)
% Note: use rcos_win.m for frequency domain RC pulse generation
%
% Example: 
%
% md.type = 'RRC';
% md.Tp = 1e-9;
% md.beta = 0.6;
% tau = (-5*md.Tp:1e-12:5*md.Tp)';
% tau_0 = 1.23456e-12;
% [s, ds, dds] = generatePulse(md, tau_0, tau, 0);
% figure; plot(tau, s); figure; plot(tau, ds); figure; plot(tau, dds);


% Paul Meissner, Dec. 2012, update Feb. 2016 Stefan Hinteregger

%% Initialize and preliminaries
if( ~(strcmp(md.type, 'RC') || strcmp(md.type, 'RRC')) )
  error('Only RC and RRC pulse supported!')
end

tau = tau(:);

if length(tau) > 1
  tau_res = tau(2) - tau(1);
else
  tau_res = 0;
end
%% RC
if strcmp(md.type, 'RC')
  t = tau-tau_0;
  i=find(t==0);
  % RC pulse
  a = sin(pi*t/md.Tp);
  b = cos(pi*md.beta*t/md.Tp);
  c = pi*t/md.Tp;
  d = 1-4*md.beta^2*t.^2/md.Tp^2;
  
  s = a.*b./(c.*d);
  s(abs(abs(t) - 0) < tau_res/1e3) = 1;
  s(abs(abs(t) - md.Tp/2/md.beta) < tau_res/1e3) = sinc(1/2/md.beta)*pi/4;

  if nargout>1
    adot = pi/md.Tp*cos(pi*t/md.Tp);
    bdot = -pi*md.beta/md.Tp*sin(pi*md.beta*t/md.Tp);
    cdot = pi/md.Tp;
    ddot = -8*md.beta^2*t/md.Tp^2;
    ds = ((adot.*b+bdot.*a).*c.*d - a.*b.*(cdot.*d+ddot.*c))./(c.*d).^2;
    ds(i) = 0;
  end
  
  if nargout>2
    adotdot = -pi^2/md.Tp^2*sin(pi*t/md.Tp);
    bdotdot = -pi^2*md.beta^2/md.Tp^2*cos(pi*md.beta*t/md.Tp);
    cdotdot = 0;
    ddotdot = -8*md.beta^2/md.Tp^2;
    
    e = adot.*b.*c.*d + a.*bdot.*c.*d - a.*b.*cdot.*d - a.*b.*c.*ddot;
    f = (c.*d).^2;
    
    edot = adotdot.*b.*c.*d + adot.*bdot.*c.*d + adot.*b.*cdot.*d + adot.*b.*c.*ddot ...
      + adot.*bdot.*c.*d + a.*bdotdot.*c.*d + a.*bdot.*cdot.*d + a.*bdot.*c.*ddot ...
      - adot.*b.*cdot.*d - a.*bdot.*cdot.*d - a.*b.*cdotdot.*d - a.*b.*cdot.*ddot ...
      - adot.*b.*c.*ddot - a.*bdot.*c.*ddot - a.*b.*cdot.*ddot - a.*b.*c.*ddotdot;
    fdot = 2*c.*d.*(cdot.*d + ddot.*c);
    
    dds = (edot.*f - fdot.*e)./f.^2;
    idx = find(abs(t) < tau_res/1e3);
    for i = 1:length(idx)
      dds(idx(i)) = (dds(idx(i)-1) + dds(idx(i)+1))/2;
    end
  end
%% RRC
elseif strcmp(md.type, 'RRC')
  t = tau-tau_0;
  i=find(t==0);
%   t(i)= 1;
  tau_res = tau(2)-tau(1);
  
  a = 4*md.beta*cos(pi*(1+md.beta)*t/md.Tp);
  b = sin(pi*(1-md.beta)*t/md.Tp);
  c = md.Tp./t;
  d = 1-16*md.beta^2*t.^2/md.Tp^2;
  
  s = 1/(pi*sqrt(md.Tp)) * (a+b.*c)./d;
  s(i) = 1/sqrt(md.Tp)*(1 - md.beta + 4*md.beta/pi);
  s(abs(abs(t) - md.Tp/4/md.beta) < tau_res/1e3) = ...
    md.beta/sqrt(2*md.Tp)*( (1+2/pi)*sin(pi/4/md.beta) + (1-2/pi)*cos(pi/4/md.beta) );
  if nargout > 1
    adot = -4*pi*(1+md.beta)*md.beta/md.Tp*sin(pi*(1+md.beta)*t/md.Tp);
    bdot = pi*(1-md.beta)/md.Tp*cos(pi*(1-md.beta)*t/md.Tp);
    cdot = -md.Tp./t.^2;
    ddot = -32*md.beta^2*t/md.Tp^2;
    
    ds = 1/(pi*sqrt(md.Tp)) * ((d.*adot-a.*ddot)./d.^2 + (d.*(bdot.*c+cdot.*b)-b.*c.*ddot)./d.^2);
    ds(i) = 0;
    idx = abs(abs(t)-md.Tp/2) < tau_res/1e3;
    for i = 1:sum(idx)
      ds(find(idx,i)) = (ds(find(idx,i)-1) + ds(find(idx,i) + 1))/2;
    end
  end
  
  if nargout > 2
    adotdot = -4*pi^2*md.beta*(1+md.beta)^2/md.Tp^2*cos(pi*(1+md.beta)*t/md.Tp);
    bdotdot = -pi^2*(1-md.beta)^2/md.Tp^2*sin(pi*(1-md.beta)*t/md.Tp);
    cdotdot = 2*md.Tp./t.^3;
    ddotdot = -32*md.beta^2/md.Tp^2;
    
    e = adot.*d-a.*ddot;
    f = bdot.*c.*d+b.*cdot.*d-b.*c.*ddot;
    g = d.^2;

    edot = adotdot.*d-a.*ddotdot;
    fdot = bdotdot.*c.*d + bdot.*cdot.*d + bdot.*c.*ddot ...
      + bdot.*cdot.*d + b.*cdotdot.*d + b.*cdot.*ddot ...
      - bdot.*c.*ddot - b.*cdot.*ddot - b.*c.*ddotdot;
    gdot = 2*d.*ddot;
    
    dds = 1/(pi*sqrt(md.Tp))*((edot.*g-gdot.*e)./g.^2 + (fdot.*g-gdot.*f)./g.^2);
    idx = find(abs(t) < tau_res/1e3);
    for i = 1:length(idx)
      dds(idx(i)) = (dds(idx(i)-1) + dds(idx(i)+1))/2;
    end
    
  end

end

%% If desired, normalize
if(nargin > 3)
  if(normalize == 0)
    if strcmp(md.type, 'RC')
      normFactor = sqrt(0.875*md.Tp);
    else
      normFactor = 1;
    end
  elseif(normalize == 1)
    normFactor = max(abs(s));
  elseif(normalize == 2)
    normFactor = sqrt(sigEnergy(s));
  elseif(normalize == 3)
    normFactor = sqrt(sigEnergy(s, tau));
  end
  
  s = s./normFactor;
  if nargout>1
    ds = ds./normFactor;
  end
  if nargout >2
    dds = dds./normFactor;
  end
  
end

%% compute 



