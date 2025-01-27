% 
% PlotHydrogenMolecularOrbital is a function designed to plot a hydrogen
% orbital. The orbital is defined by the quantum numbers, which in tern
% determine the number of nodes and harmonic frequencies observed in the
% wavefunction. This program is intended both to display the functions and
% to allow the user to mix and match wave functions (such as for hybrid
% orbitals).
% 
% Radial nodes are determined by Laguerre polynomials, and angular nodes by 
% Legendre polynomials. Details of the calculation are in the radial wave 
% function and angular wave function sections.
% 
% -------------------------------------------------------------------------
% Some notes on the calculation
%  
% As the orbital increases in size, increase cons.spatialLen to see it all
% Increase cons.resolution to make the plot prettier 
%     remember this makes an N^3 increase in computional time
% 
% -------------------------------------------------------------------------
% notes on the plots
% 
% the functions PlotDensityIsosurface, PlotCrossSectionIsosurface, and
%     PlotWaveFnSignIsosurface are a good start to visualizing the data
% PlotCrossSectionIsosurface has a color bar. 
% 

%% main function

% this function constructs an atomic orbital based on the spherical hamonic
% simple interpretaion of the hydrogen atom.

function [WaveFn, r] = Hydrogen_2p
%--------------------------------------------------------------------------
% define constants and orbital shape

% quantum numbers
n = 2;           % principal quantum number
l = 1;           % angular momentum quantum number (subshell type: s, p, d, f, etc.)

% constants: shape and cutoffs of 3D plots (needs fiddling to get a good plot)
cons.resolution = 201;             % number of calculations per XYZ dimension   
cons.spatialLen = 5;             % length of physical space to view over (angstroms)

% constants: scale size of plot (don't play with)
cons.a0         = 0.53;            % Bohr radius (meters)

%--------------------------------------------------------------------------
% generate XYZ space, convert to spherical radius, theta, and phi

[xSpace, ySpace, zSpace] = make3Dspace(cons);
[theta, phi, r] = convert2polar(xSpace, ySpace, zSpace);

%--------------------------------------------------------------------------
% plot figure

% An index for the different orbital orientations with the same n and l
k = 1;

WaveFn = psiCalculation(n,l,k,r,theta,phi,cons);
cons.orbitName = orbitNameIs(n,l,k);
%PlotDensityIsosurface(xSpace, ySpace, zSpace, WaveFn, cons);
PlotCrossSectionIsosurface(xSpace, ySpace, zSpace, WaveFn, cons);
%PlotWaveFnSignIsosurface(xSpace, ySpace, zSpace, WaveFn, cons);

% lighting, aspect ratio, plot options
daspect([1 1 1]);
view(3); 
axis vis3d;
camlight right;
camlight headlight;
lighting gouraud;
xlabel('x','FontName','Arial','FontSize',18);
ylabel('y','FontName','Arial','FontSize',18);
zlabel('Angstroms','FontName','Arial','FontSize',18);
rotate3d on;

% dimensions of plot
plotScale = cons.spatialLen;
xlim([-plotScale plotScale]);
ylim([-plotScale plotScale]);
zlim([-plotScale plotScale]);
end



%==========================================================================
%%                      calculate wave function
%==========================================================================

% the wave function is a 3D vector of size [resolution, resolution, resolution]
% calculate it by multiplying each element of the radial and angular components
function [WaveFn] = psiCalculation(n,l,k,r,theta,phi,cons)

% WaveFn = Radial part * Angular part
WaveFn = radialCalculation(n,l,r,cons) .* angularCalculation(l,k,theta,phi);

% correction 1: remove NaN at atomic nucleus
center = ceil(cons.resolution/2);
WaveFn(center, center, center) = 0;

% correction 2: flip wave function properly when m < 0
if (k < 0)
    WaveFn = permute(WaveFn, [2 1 3]);
end

WaveFn = -WaveFn;
end

%% radial wave function

% the radial component is composed of two parts: an exponential term and a
% polynomial term. The exponential term adds attraction to the nucleus of
% the atom. The polynomial term adds electron shell harmonics to create
% nodal spheres

function [RadialFn] = radialCalculation(n,l,r,cons)

% import contants
a0 = cons.a0;

% scaling factors
coeff1 = sqrt((2/(a0*n))^3 * factorial(n-l-1)/(2*n*factorial(n+l)));
coeff2 = 1/factorial(n + l);

% Part 1:  exponential component (attraction to nucleus)
nuclearComponent = (2*r/(a0*n)) .^ l .* exp(-r/(a0*n));

% Part 2:  polynomial component (generates radial nodes)
radialNodeComponent = LaguerrePolynomial(n-l-1, 2*l+1, 2*r/(a0*n));

% combine components to calculate radial wave function
RadialFn = coeff1 * coeff2 * nuclearComponent .* radialNodeComponent;

end

% use Laguerre polynomials to introduce nodal spheres into radial function
function [NodalFn] = LaguerrePolynomial(n,m,r)

% initiate polynomial function
NodalFn = zeros(size(r));

% add n coefficient terms to the polynomial
for i = 0 : n
    polynomialCoeff = factorial(m+n) * nchoosek(m+n,n-i) / factorial(i);
    NodalFn = NodalFn + polynomialCoeff * (-r).^i;
end

end

%% angular wave function

% the angular component is a cosine function with Legendre polynomials used
% to add nodal planes.

function [AngularFn] = angularCalculation(l,k,theta,phi)

if (k == 2)
    phi = phi + pi/4;
end

% normalization and scaling factors
normFac = abs(sign(k)*(2^0.5) + (sign(abs(k)) - 1)*2);
scalFac = sqrt((2*l+1) / (4*pi) * factorial(l-abs(k)) / factorial(l+abs(k)));

% add nodes to spherical harmonics functions
SphFn1 = scalFac * LegendrePolynomial(l,k,cos(theta))  .* exp(1i*k*phi);
SphFn2 = scalFac * LegendrePolynomial(l,-k,cos(theta)) .* exp(1i*-k*phi);

AngularFn = (SphFn1 + SphFn2) / normFac;

end

% use Legendre polynomials to introduce nodal planes into angular function
function [NodalFn] = LegendrePolynomial(l,m,x)

% initiate polynomial function
NodalFn = zeros(size(x));
numCoeffs = floor(1/2*l - 1/2*abs(m));

% add n coefficient terms to the polynomial
for n = 0 : numCoeffs
    polynomialCoeff = (-1)^n * nchoosek(l-2*n,abs(m)) * nchoosek(l,n) * nchoosek(2*l-2*n,l);
    exponent = l - 2*n - abs(m);
    NodalFn = NodalFn + polynomialCoeff * x.^exponent;
end
NodalFn = (-1)^m * (1-x.^2).^(abs(m)/2) .* (factorial(abs(m))/2^l*NodalFn);

end

%==========================================================================
%%               generate 3D space (X Y Z theta phi r)
%==========================================================================

% generate 3D cartesian space based on inputs
function [xSpace, ySpace, zSpace] = make3Dspace(cons)

% import dimensions and scaling factors
resolution = cons.resolution;
spatialLen = cons.spatialLen;

% generate XYZ space using meshgrid
xRange = linspace(-spatialLen, spatialLen, resolution);
yRange = linspace(-spatialLen, spatialLen, resolution);
zRange = linspace(-spatialLen, spatialLen, resolution);
[xSpace, ySpace, zSpace] = meshgrid(xRange, yRange, zRange);

end

% most of this function is equivalent to the "cart2sph" function.
% Matlab interprets angles a little differently than a classic calculus
% textbook, as in the engineering-style of polar coordinatles.
% phi and theta are swapped, and the zero-angle is at a -pi/2 angle. 
% Consequently we manually convert the points here.
function [theta, phi, r] = convert2polar(x, y, z)

r     = sqrt(x.^2+y.^2+z.^2);
theta = acos(z./r);     % this is "phi" in the "cart2sph" function
phi   = atan2(y,x);     % this is "theta" in the "cart2sph" function

end

%==========================================================================
%%                      plot wave functions
%==========================================================================

% convert n l m into a useful name
function [orbitName] = orbitNameIs(n,l,k)
subShells = {'s', 'p', 'd', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm'};

orbitName = [num2str(n) subShells{l+1}];
end

% plot a solid isosurface of the orbital
function PlotDensityIsosurface(xSpace, ySpace, zSpace, WaveFn, cons)

% import dimensions and scaling factors
orbitName  = cons.orbitName;

% square the wavefunction to make an electron density map
WaveFn2 = abs(WaveFn) .^ 2;

% maximum value
vmax = max(max(max(WaveFn2)));

% patch and interpret isosurface for better coloration
cmap = jet(21);

isovalue = 0.01 * vmax;
isoColor   = cmap(3, :);
p = patch(isosurface(xSpace,ySpace,zSpace,WaveFn2,isovalue));
isonormals(xSpace,ySpace,zSpace,WaveFn2,p);
set(p,'FaceColor',isoColor,'EdgeColor','none');

caxis([0 vmax]);
title([orbitName ' orbital'],'FontName','Arial','FontSize',18);

end

% plot a cross section of a solid isosurface of the orbital
function PlotCrossSectionIsosurface(xSpace, ySpace, zSpace, WaveFn, cons)

% import dimensions and scaling factors
orbitName   = cons.orbitName;

% slice off half the isosurface, change XYZ accordingly
halfLength = 0;
clipRange = [nan,nan,halfLength,nan,nan,nan];
%clipRange = [halfLength,nan,nan,nan,nan,nan];
[x,y,z,v] = subvolume(xSpace,ySpace,zSpace,WaveFn,clipRange);

% maximum value
%absv = abs(v);
vmax = max(max(max(v)));

% patch and interpret isosurface for better coloration
cmap = jet(21);

isovalue = 0.1 * vmax;
%p1 = patch(isosurface(x,y,z,absv,isovalue),'FaceColor',cmap(3,:),'EdgeColor','none');
%isonormals(x,y,z,absv,p1);
%patch(isocaps(x,y,z,absv,isovalue),'FaceColor','interp','EdgeColor','none');

p1 = patch(isosurface(x,y,z,v,isovalue),'FaceColor',cmap(13,:),'EdgeColor','none');
isonormals(x,y,z,v,p1);
patch(isocaps(x,y,z,v,isovalue,'above'),'FaceColor','interp','EdgeColor','none');

p2 = patch(isosurface(x,y,z,v,-isovalue),'FaceColor',cmap(9,:),'EdgeColor','none');
isonormals(x,y,z,v,p2);
patch(isocaps(x,y,z,v,-isovalue,'below'),'FaceColor','interp','EdgeColor','none');

% nodal plane
p0 = patch(isosurface(x,y,z,v,0),'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');
isonormals(x,y,z,v,p0);

% lighting, aspect ratio, plot options
caxis([-vmax vmax] * 0.6);
title([orbitName ' orbital cross section'],'FontName','Arial','FontSize',18);
colormap jet;
colorbar;

end

% plot a solid isosurface of the orbital colored by the sign of the
% wavefunction.
function PlotWaveFnSignIsosurface(xSpace, ySpace, zSpace, WaveFn, cons)

% import dimensions and scaling factors
orbitName  = cons.orbitName;
isoColor   = sign(WaveFn);

% square the wavefunction to make an electron density map
WaveFn2 = abs(WaveFn) .^ 2;

% maximum value
vmax = max(max(max(WaveFn2)));

% make isosurface
isovalue = 0.01 * vmax;
isosurface(xSpace,ySpace,zSpace,WaveFn,isovalue,isoColor)

caxis([-1 1]);
title([orbitName ' orbital'],'FontName','Arial','FontSize',18);

end


