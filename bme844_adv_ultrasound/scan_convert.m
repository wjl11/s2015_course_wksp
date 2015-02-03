function [out,ax,l,outside] = scan_convert(iin,mode,min_phi,span_phi,apex,flag,fs,varargin)
%
% SCAN_CONVERT
%
% BEWARE OF UNITS !!!! BEWARE OF UNITS !!!! 2DBEWARE OF UNITS !!!!
% converts A-line data and geometric information to a scan-converted image
%
% EXAMPLES:
%
% All coordinates and default spacing (.2mm)
% [out,ax,la,outside]=scan_convert(in,mode,min_phi,span_phi,apex,interp)
%
% User Defined Coordinates ans Spacing
% [out,ax,la,outside]=scan_convert(in,mode,min_phi,span_phi,apex,interp,[ax_min ax_max ax_inc lat_min lat_max lat_inc]);
% 
% Type this...
% [out,ax,la,outside]=scan_convert(bimage,'sector',min(blat),max(blat)-min(blat),fh.rfbm.ApexVerticalCm,2,[0 0 2e-5 0 0 1e-4]); 
%
% INPUT:
%  iin = image in
%  mode = scan type ('sector')
%  min_phi = minimum angle (degrees)
%  span_phi = angle span (degrees)
%  apex = distance to radial center (cm)
%  fs = sampling frequency; 40 MHz is a default value
%  flag = 1 (nearest neighbor) 2 (bilinear interpolation)
%  vargin = [ ax_min ax_max ax_inc lat_min lat_max lat_inc ] (meters)
%
% OUTPUT:
%  out = image
%  ax = axial coordinates
%  lat = lateral coordinates
%  outside = ???

% sampling frequency default value 
if (nargin < 7) || isempty(fs)
    fs =40e6;
end

% initialize dimensions
q = [0 0 0 0 0 0];

% insert user defined dimensions
if (nargin>7) 
  q = varargin{7-6};

  % allows code to be compatible with earlier versions
  if (length(q)==5)
    q(4:6)=q(3:5);
    q(3) = q(6);
  end;

end;

% calculate scan conversion map
[idx,i00,dr,dth,l,ax]=scmap(size(iin),'sector',min_phi,span_phi,apex,1,fs,q);

% bottom left index
i10 = i00+1;

% top right index
i01 = i00+size(iin,1)+1;

% bottom right index
i11 = i01+1;

% initialize output array
out = ones([length(ax) length(l) size(iin,3)])*min(iin(:));

% for each frame...
for(i=1:size(iin,3))
  % grab single frame
  in = iin(:,:,i);

  % extend image by 1 row and 1 column (allows interpolation)
  in(end+1,:)=in(end,:);
  in(:,end+1)=in(:,end);

  % initialize output frame
  imageout = out(:,:,i);

  % nearest neighbor
  imageout(idx) = in(i00);

  % bilinear interpolation
  if(flag==2)
    imageout(idx) = (in(i10)-in(i00)).*dr+(in(i01)-in(i00)).*dth+...
                    (in(i11)+in(i00)-in(i10)-in(i01)).*dth.*dr+...
		    imageout(idx);
  end;

  % assign output frame
  out(:,:,i) = imageout;

end;

outside = 0;
