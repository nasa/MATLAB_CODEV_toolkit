function   [raygridx,raygridy,centx,centy,rms_spot_size,size_x,size_y] = cvspot(f,nrd,s,w,z,center)
%CVSPOT  returns x,y spot diagram data
%
%   function   [raygridx,raygridy,rmsx,rmsy] = cvspot(f,nrd,s,w,z,center);
%
%   INPUTS: f : can be a vector of f numbers (default = 1)
%           w, z: (defaults =1)
%           nrd: number of rays across diameter to be traced.
%           center: choose 1,2, or 3.  
%               1 = no centering,
%               2 = centering on zero, 
%               3 = centering on zero and spacing of fs
%               4 = no display
%           surface = surface at which data is returned (def = image)
%
%   OUTPUT: If there are no output arguments, the output is a scatter 
%               plot showing the spot diagram
%           raygridx: grid of x coordinates at specified suface
%           raygridy: grid of y coordinates at specified suface
%           centx: x spot centroid (length = # of fs input) 
%           centy: y spot centroid (length = # of fs input)
%           rms_spot_size
%
%   See also CVPMA, CVPSF
%

if nargin <1, f = 1; end
if nargin <2, nrd = 25; end
if nargin <3, s = cvnum; end
if nargin <4, w = 1; end
if nargin <5, z = 1; end
if nargin <6, center = 1; end
colors = {'r.','y.','g.','b.','m.'};


for jj=1:length(w)
    for ii=1:length(f)
        pause(0.3);
        [data,mask(:,:,ii,jj)] = cvrgrid(nrd,ii,s,jj,z);
        raygridx(:,:,ii,jj) = data(:,:,1);
        raygridy(:,:,ii,jj) = data(:,:,2);

        raygridx(:,:,ii,jj) = raygridx(:,:,ii,jj).*mask(:,:,ii,jj);
        raygridy(:,:,ii,jj) = raygridy(:,:,ii,jj).*mask(:,:,ii,jj);

        rx = raygridx(:,:,ii,jj); ry = raygridy(:,:,ii,jj);
        centx(ii,jj) = mean(rx(rx~=0));
        if length(find(rx(:)~=0))==0, disp(['no data x ' num2str(ii)]);end
        centy(ii,jj) = mean(ry(ry~=0));
        if length(find(ry(:)~=0))==0, disp(['no data y ' num2str(ii)]);end

        %     rmsopd = norm(opd(mask)-mean(opd(mask)))./sqrt(length(find(mask)));

        rms_spot_size(ii,jj) = sqrt(mean( (rx(rx~=0)-mean(rx(rx~=0))).^2 + (ry(ry~=0)-mean(ry(ry~=0))).^2 ));
        size_x(ii,jj) = sqrt(mean( (rx(rx~=0)-mean(rx(rx~=0))).^2  ));
        size_y(ii,jj) = sqrt(mean( (ry(ry~=0)-mean(ry(ry~=0))).^2 ));
    end
end

if nargout<1 %plot data if no output
    if center==1; center = 'one'; % No centering of fs to zero or anything else
    elseif center==2; center = 'two'; % Remove mean, center on zero (use only for one f point)
    elseif center==3; center = 'three';
    elseif center==4; center = 'four';
    end

    switch center
        case 'one'
            figure; axis equal
            for ii=1:length(f);
            for jj=1:length(w);
                rx = raygridx(:,:,ii,jj);
                ry = raygridy(:,:,ii,jj);
                if length(find(rx))<length(find(ry)); len1=rx;else len1=ry; end
                hold on
                plot(rx(find(len1)),ry(find(len1)),colors{jj},centx(ii,jj),centy(ii,jj),'kh');
                clear rx ry len1
            end
            end
        case 'two'
            for ii=1:length(f);
            for jj=1:length(w);
                rx = raygridx(:,:,ii,jj);
                ry = raygridy(:,:,ii,jj);
                if length(find(rx))<length(find(ry)); len1=rx;else len1=ry; end
                figure(10); hold on
                scatter(rx(find(len1))-mean(rx(find(len1))),ry(find(len1))-mean(ry(find(len1))));
                clear rx ry len1
            end
            end
        case 'three'
            temp = raygridx(:,:,1);
            offset = 20*std(temp(find(temp))); clear temp
            for ii=1:length(f);
            for jj=1:length(w);
                rx = raygridx(:,:,ii,jj);
                ry = raygridy(:,:,ii,jj);
                if length(find(rx))<length(find(ry)); len1=rx;else len1=ry; end
                roffy = offset*ii;
                figure(10); hold on; axis equal
                scatter(rx(find(len1))-mean(rx(find(len1))),...
                    ry(find(len1))-mean(ry(find(len1)))+ roffy);
                clear rx ry roffy len1
            end
            end
            clear roffy
        case 'four'
    end
    hold off
end

if nargout==1
    raygridx = [raygridx(:) raygridy(:)];
end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     