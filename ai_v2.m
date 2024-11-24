%% Function ai_v2
function [ai, Pi, Tfoot, Ti, Tmax, typetxt] = ai_v2(p,samplerate)
% Originally written as a script by ADH 12.6.01
% reconfigured as a function 19/04/19
% Bug fix suggested by Richard Scott implemented 11.04.24
% Calculates AIx according to Kelly method and
% classifies waveform according to Murgo et al.
% Uses zero 4th derivative to determine first shoulder
% AI = (max pressure - inflection pressure)/pulse pressure (Ps-Pi)/(Ps-Pd) %
% References
% Kelly R, Hayward C, Avolio A, O'Rourke M. Noninvasive determination
% of age-related changes in the human arterial pulse. Circulation 1989; 80(6): 1652-9.
% Murgo JP, Westerhof N, Giolma JP, Altobelli SA. Aortic input impedance in normal man:
% relationship to pressure wave forms. Circulation 1980; 62(1): 105-16.
% Nichols & O'Rourke. McDonald's Blood Flow in Arteries 1998.
% M. Karamanoglu. Diagnostic Applanation Tonometry 1996.
%%
samples=1:length(p);
dpdt=fsg71(p);
d2p=fsg71(dpdt);
d3p=fsg71(d2p);
d4p=fsg71(d3p);
% d5p=fsg71(d3p);       % not used currently
%ndp=dpdt./max(dpdt);
%nd2p=d2p./max(d2p);
%nd3p=d3p./max(d3p);
nd4p=1e-6+fix(10*(d4p./max(d4p))); % this filters out minor crossings (<10%)
%np=(p-min(p))./(max(p)-min(p));
d4ptable=[samples;nd4p].'; % create column oriented table from data
zerod4p=round(mminterp(d4ptable,2,0)); % establish zero crossings
zcross=zerod4p(:,1);
%zeroslope=d5p(zcross); % establish where zeros correspond to negative slope d5p - not used so far

%Define reference points with regard to 4th derivative (Kelly et al., 1989)
[pmax, tmax]=max(p);
tfoot=zcross(1);
ti=zcross(3);

% Allow for type B & C otherwise Pi may = Ps.  The choice of 5 samples (0.025ms) is arbitrary
if tmax-ti<5
    ti=zcross(4);
end

pfoot = p(tfoot); % as defined by O'Rourke p216 and Kelly 1989
Pi=p(ti); % as defined by O'Rourke p216 and Kelly 1989
PsPd=pmax-pfoot; % as defined by O'Rourke p216
PsPi=pmax-Pi;
Tfoot=tfoot/samplerate;     % calculate in s
Ti = ti/samplerate;
Tmax=tmax/samplerate;
ai=round(PsPi/PsPd*100);

% Identify Type C
if tmax<ti
    ai=-ai;
end

% Name as type (A, B or C) based on criteria in Murgo et al., Circ 1980;
% 62: 105-116. Murgo et al. aren't very clear about type B but by exclusion
% ai<=12 for type B.
if tmax>=ti
    if ai>12
        typetxt=('Type A');
    elseif ai<=12
        typetxt=('Type B');
    end
elseif tmax<ti
    typetxt=('Type C');
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS
%%
function y=mminterp(tab,col,val)
% MMINTERP 1-D Table Search by Linear Interpolation.
% Y=MMINTERP(TAB,COL,VAL) linearly interpolates the table
% TAB searching for the scalar value VAL in the column COL.
% All crossings are found and TAB(:,COL) need not be monotonic.
% Each crossing is returned as a separate row in Y and Y has as
% many columns as TAB. Naturally, the column COL of Y contains
% the value VAL. If VAL is not found in the table, Y=[].
% Based on code by D.C. Hanselman,
% University of Maine, Orono ME,  04469
% 1/26/94
% Copyright (c) 1996 by Prentice-Hall, Inc.

[rt,ct]=size(tab);
if length(val)>1, error('VAL must be a scalar.'), end
if col>ct||col<1,  error('Chosen column outside table width.'), end
if rt<2, error('Table too small or not oriented in columns.'), end
above=tab(:,col)>val;
below=tab(:,col)<val;
equal=tab(:,col)==val;
if all(above==0)||all(below==0) % handle simplest case
    y=tab(equal,:);
    return
end
pslope=find(below(1:rt-1)&above(2:rt)); %indices where slope is pos
nslope=find(below(2:rt)&above(1:rt-1)); %indices where slope is neg

ib=sort([pslope;nslope+1]);	% put indices below in order
ia=sort([nslope;pslope+1]);	% put indices above in order
ie=find(equal);				% indices where equal to val

[tmp,ix]=sort([ib;ie]);		% find where equals fit in result
ieq=ix>length(ib);			% True where equals values fit
ry=length(tmp);				% # of rows in result y

y=zeros(ry,ct);				% poke data into a zero matrix

alpha=(val-tab(ib,col))./(tab(ia,col)-tab(ib,col));
alpha=alpha(:,ones(1,ct));
y(~ieq,:)=alpha.*tab(ia,:)+(1-alpha).*tab(ib,:);	% interpolated values

y(ieq,:)=tab(ie,:);			% equal values
y(:,col)=val*ones(ry,1);	% remove roundoff error
end
%%
function dx=fsg71(x)
%% dx=fsg71(x)
% 7 point SavGol filter, 1st derivative
% input x
% output dx
% corrected for time shift
%% 2nd order polynomial
C=[.107143,.071429,.035714];
%% 3rd order polynomial
%C=[-.087302,.265873,.230159];
B=zeros(1,3);
for i=1:3
    B(i)=C(i);
end
B(4)=0.;
for i=5:7
    B(i)=-C(8-i);
end
A=[1,0];
s=size(x,2);
dx=filter(B,A,x);
dx=[0,0,0,dx(7:s),0,0,0];
end
