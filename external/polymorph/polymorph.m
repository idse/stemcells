function P=polymorph(P1,P2,numsteps,numpts,exact)
%POLYMORPH Morphs two polygons
%   P=POLYMORPH(P1,P2) returns an array for the morphing of polygon P1 to
%   polygon P2. Both polygons must be two dimensional but can have
%   different number of points. 
%
%   POLYMORPH(P1,P2,NUMSTEPS) morphs the polygons P1 into polygon P2 within
%   a specified NUMSTEPS number of steps. (Default is 10)
%
%   POLYMORPH(P1,P2,NUMSTEPS,NUMPTS) interpolates the two polygons so they
%   have at least NUMPTS number of points. All data points in original
%   input polygons are contained in the output polygons, while also
%   ensuring interpolated polygons have the same number of points. The
%   algorithm uses the least common multiple of #P1 and #P2. 
%   If NUMPTS<LCM then LCM is used. (Default is 1, implying use LCM)
%
%   POLYMORPH(P1,P2,NUMSTEPS,NUMPTS,EXACT) with non-zero value for EXACT
%   if you want to have exactly NUMPTS number of points. This is beneficial
%   for a couple scenarios:
%      - Splicing multiple sets of polygons together that may have
%        different number of points between each polygon. P1, P2, P3,...
%      - Forced downsampling for larger polygons where the overall shape is
%        important, but not the inclusion of every original data point.
%        This will save computation time.
%   Note: Due to forcing a new constraint on number of output points, 
%   some of the original data points may not appear in the output
%   interpolated polygons when EXACT is non-zero. (Default is 0)
%
%   Algorithm will automatically try and solve for a morphing that produces
%   no self-intersections through time, and yields a minimum total
%   distance-squared for the path of each point. If no such path is found,
%   the best self-intersecting morph is given.
%
%   With more NUMPTS a higher chance of discovering a non-intersecting
%   morph is found, however the intermediate shapes might not be preserved
%   - that is, from a triangle to a triangle an intermediate shape might
%   be a quadrilateral.
%
%   Note: For the output, the last point for each polygon is identical to
%   the first point, for closure. If EXACT=1 there will be NUMPTS rows of
%   points given, but NUMPTS-1 unique points in the polygons.
%
%   Output is an array where
%   P(:,:,1)          P1 polygon (with interpolated points)
%   P(:,:,t)          an intermediate polygon
%   P(:,:,NUMSTEPS)   P2 polygon (with interpolated points)
%
%
%   EXAMPLE:
%
%   numsteps=20; numpts=10;
%   P1=rand(3,2);
%   x2=rand(8,1); y2=rand(8,1); k2 = convhull(x2,y2); P2=[x2(k2) y2(k2)];
%   P=polymorph(P1,P2,numsteps,numpts);
%   for i=1:size(P,3)
%      hold on, plot(P(:,1,i),P(:,2,i)), axis([0 1 0 1]), pause(0.1)
%   end
%

%   Mike Sheppard
%   Last Modified 14-Jul-2011


%Check to see if required program is available
%CAN COMMENT OUT AFTER DOWNLOADED
if ~exist('intersections.m','file')
    str1='Requires program #11837 from File Exchange';
    str2='"Fast and Robust Curve Intersections"';
    str3='Available to download here:';
    str4='http://www.mathworks.com/matlabcentral/fileexchange/11837';
    errorstr=strcat(str1,' \n ',str2,' \n ',str3,' \n ',str4);
    error('polymorph:RequiredProgram', errorstr);
end

%Error Checking
if nargin < 2
    error('polymorph:TooFewInputs',...
        'Requires at least two input arguments.');
end
if nargin<3
    numsteps=10;
end
if nargin<4
    numpts=1; %Default to LCM since 1<=LCM for all cases
end
if nargin<5
    exact=0;
end
if (size(P1,2)~=2 || size(P2,2)~=2)
    error('polymorph:InvalidInput',...
        'Polygons must be nx2 matrix, but may have different number of rows for each.');
end
if numsteps<1 || numsteps~=round(numsteps) || numpts<1 || numpts~=round(numpts)
    error('polymorph:InvalidInput',...
        'NUMSTEPS and NUMPTS must be positive integers.');
end

%Make exact boolean
exact=(exact~=0);

%Remove any rows where either coordinate is not finite
k1=unique([~isfinite(P1(:,1)) ~isfinite(P1(:,2))]); P1(k1,:)=[];
k2=unique([~isfinite(P2(:,1)) ~isfinite(P2(:,2))]); P2(k2,:)=[];
if (size(P1,1)<3 || size(P2,1)<3)
    error('polymorph:IncorrectInput',...
        'Requires both polygons to be at least 3 points.');
end

%Delete possible closure end point
[P1,P2]=delendpt(P1,P2);

%Number of points to use
%Use NUMPTS if EXACT=1, else use LCM method
if ~exact
    lc=lcm(size(P1,1),size(P2,1));  %Least Common Multiple
    numpts=1+(lc*ceil(numpts/lc));  %Use at least NUMPTS number of points
    %+1 due to extra point of last point=first point to make the polygon closed
end

%Add repeated end point to interpolate closed polygon
[P1,P2]=addendpt(P1,P2);

%Interpolate both polygons with NUMPTS number of points
numvec=[size(P1,1) size(P2,1)];
num1v=linspace(1,numvec(1),numpts);
num2v=linspace(1,numvec(2),numpts);
P1=[interp1(1:numvec(1),P1(:,1),num1v); interp1(1:numvec(1),P1(:,2),num1v)]';
P2=[interp1(1:numvec(2),P2(:,1),num2v); interp1(1:numvec(2),P2(:,2),num2v)]';

%Delete closure end point for algorithm to cycle through
[P1,P2]=delendpt(P1,P2);

%Find matching points so that total distance-squared traveled among all
%points is a minimum
%'rotate' the indices of P2 until best match between P1 and P2 is obtained
%Go both 'clockwise' and 'counterclockwise' direction
%Since either polygon may be in either orientation, CW and CCW are not
%used; but normal-index and "flipped" index ([1 2 3] vs [3 2 1])
minsumd=[Inf Inf]; flipindx=[NaN NaN]; indx=[NaN NaN];
%Notation: [Best with Intersections Allowed; Best with no Intersections]
X=P1; Y1=P2; Y2=flipud(P2);
for k=1:numpts
    tempY1=Y1([k:end 1:k-1],:); D1=sum(sum((X-tempY1).^2));
    tempY2=Y2([k:end 1:k-1],:); D2=sum(sum((X-tempY2).^2));
    if any(D1<minsumd) || any(D2<minsumd) %New shorter distance
        [D,I]=min([D1 D2]);
        minsumd(1)=D; indx(1)=k; flipindx(1)=I-1; %D1 is no flip, D2 is flip
        %Now check to see if either path yields no intersections
        if D1<minsumd(2)
            P=allpoly(P1,tempY1,numpts,numsteps); ind1=checkinter(P,numsteps);
            if ~ind1
                minsumd(2)=D1; indx(2)=k; flipindx(2)=0;
            end
        end
        if D2<minsumd(2)
            P=allpoly(P1,tempY2,numpts,numsteps); ind2=checkinter(P,numsteps);
            if ~ind2
                minsumd(2)=D2; indx(2)=k; flipindx(2)=1;
            end
        end
        
    end
end


%Use no intersections if possible, else use best scenario with
%intersections
if ~isfinite(minsumd(2)) 
%If INF appears for no intersections, use best scenario with intersections
    flipindx=flipindx(1); indx=indx(1);
else
%A morphing exists with no intersections
    flipindx=flipindx(2); indx=indx(2);
end

%Re-order to match points to P1 to yield minimized distance
if flipindx==0
    P2=Y1([indx:end 1:indx-1],:);
else
    P2=Y2([indx:end 1:indx-1],:);
end

P=allpoly(P1,P2,numpts,numsteps);

return


function ind=checkinter(P,numsteps)
%Uses INTERSECTIONS available on File Exchange:
%http://www.mathworks.com/matlabcentral/fileexchange/11837
ind=0; k=1;
s = warning('off', 'all'); %temporarily disable any warnings
%In theory could check the paths that each point takes, to see if any paths
%cross. However, there are degenerate cases - such as a square morphing
%into another square that is parallel. This yields overlapping paths among
%most of the points from beginning location to ending location - 
%but for each time step the polygon does not actually self-intersect.
%
%The procedure chosen was to check each successive polygon
%within the morphing to see if it self-intersects; stopping at the first
%non-trivial self-intersection.
%
while (k<numsteps && ~ind)
    try
        [X0,Y0] = intersections(P(:,1,k),P(:,2,k));
        %Exclude all points in polygon
        M=[X0 Y0]; M=setdiff(M,P(:,:,k),'rows');
        %Exclude if equal to last point (which is listed twice on purpose)
        kdel=sum((M-repmat(P(end,:,k),size(M,1),1)).^2,2);
        M(kdel<100*eps,:)=[];  %If not self-intersecting M should be empty 0x2
        ind=size(M,1);
    catch errorcatch %Just in case error in program
        ind=0; %assume no intersections found
    end
    k=k+1;
end
%IND>0 there is a crossing at some point
s = warning('on', 'all');
return


function P=allpoly(P1,P2,numpts,numsteps)
%CREATE ALL INTERMEDIATE POLYGONS FOR OUTPUT
[P1,P2]=addendpt(P1,P2); %Close the final polygons
t=linspace(0,1,numsteps);
t=reshape(t,[1 1 numsteps]);
t=repmat(t,numpts,2);
a=repmat(P1,[1,1,numsteps]);
b=repmat(P2,[1,1,numsteps]);
P=a.*(1-t)+b.*t;
return


function [P1,P2]=addendpt(P1,P2)
%ADDS LAST END POINT THE SAME AS FIRST POINT
%MAKES OPEN POLYGON CLOSED
if ~isequal(P1(1,:),P1(end,:))
    P1(end+1,:)=P1(1,:);
end
if ~isequal(P2(1,:),P2(end,:))
    P2(end+1,:)=P2(1,:);
end
return


function [P1,P2]=delendpt(P1,P2)
%DELETES LAST END POINT IF IDENTICAL TO FIRST POINT
%MAKES CLOSED POLYGON OPEN
if isequal(P1(1,:),P1(end,:))
    P1(end,:)=[];
end
if isequal(P2(1,:),P2(end,:))
    P2(end,:)=[];
end
return
