function [x,y,z,avgr,iIt] = psphere1(n,isDemo)
% [x,y,z,r,count] = psphere1(N,isDemo)
%
% Distributes N points "equally" about a unit sphere.
%
% N The number of points to distribute
% x,y,z Each is 1 x N vector
% r The smallest linear distance between two neighboring 
% points. If the function is run several times for the 
% same N, r should not change by more than the convergence 
% criteria, whose default is +-0.01 on a unit sphere.
% count The number of iterations it took to converge.
% isDemo 0 or 1. If 1, the points are displayed over a unit
% sphere for each iteration. isDemo defaults to 0 if
% only N is entered.
%
% Distributes N points about a unit sphere so that the straight line
% distance between neighboring points is roughly the same. The
% actual criteria for stopping is slightly different. The difference
% between a given point and every other point is stored in a matrix.
% The iterations stop once the maximum difference between any element
% in successive distance matrices is less than 0.01. An absolute 
% criteria was chosen due to self-distances being 0, and programming 
% around this for relative convergence seemed too much work for too 
% little reward.
%
% The algorithm first generates N random points. Then a repulsive 
% force vector, based on 1/r^2, is calculated for each point due to 
% all of the other points. The resultant force vector for each point
% is normalized, and then each point is displaced a distance S = 1 
% in the direction of the force. Any value of S between 0.0 and 1.
% 0 seems to work with most values between 0.2 and 1 taking an average 
% of 20 to 30 iterations to converge. If S is too high, too much 
% "energy" is being added to the system and it won't converge. If S is
% too low, the points won't be evenly distributed even though the
% convergence criteria is met. Generally speaking, the larger N is
% the larger S can be without causing instabilities. After 
% displacement, the point is projected back down onto the unit sphere. 
% When the system nears convergence, the displacement vector for a 
% given point is nearly in the same direction as the radius vector for 
% that point due to the points being equally distributed. A good check 
% to make sure the code is working is to specify N = 4 and then check 
% that the resulting points form a regular tetrahedron (or whatever 
% it's called). How you would do this I don't know (check angles 
% maybe). That's why I supplied the demo option so you could look at 
% it in progress.
%
% Jason Bowman
% jbowman90@hotmail.com
% Last Revised June 2000
%
% 160729PJ: Enforced the first point to lie at [1,0,0]
%
%
% You may freely distribute this code. I only ask that you give credit 
% where credit is due.
%
   
   if nargin == 1,
      isDemo = 0;
   end;% if

   isSeed = 1;
   diffTol = 1e-10;% 0.01
   iIt = 0;
%   gold = (1+sqrt(5))/2;% golden ratio

   % Since rand produces number from 0 to 1, subtract off -0.5 so that
   % the points are centered about (0,0,0).
   x = rand(1,n) - 0.5;
   y = rand(1,n) - 0.5;
   z = rand(1,n) - 0.5;
   
   % Make the matrix R matrices for comparison.
   rm_new = ones(n);
   rm_old = zeros(n);

   % Scale the coordinates so that their distance from the origin is 1.
   r = sqrt(x.^2 + y.^2 + z.^2);

   x = x./r;
   y = y./r;
   z = z./r;

   if n == 12,
     x(1:4) = 0; 
     y(1:2) = abs(y(1:2)); y(3:4) = -abs(y(3:4));
     z([1 3]) = abs(z([1 3])); z([2 4]) = -abs(z([2 4]));
     z(5:8) = 0;
     x(5:6) = abs(x(5:6)); x(7:8) = -abs(x(7:8));
     y([5 7]) = abs(y([5 7])); y([6 8]) = -abs(y([6 8]));
     y(9:12) = 0;   
     z(9:10) = abs(z(9:10)); z(11:12) = -abs(z(11:12));
     x([9 11]) = abs(x([9 11])); x([10 12]) = -abs(x([10 12]));
   elseif n == 20,
     x(1) = 1/sqrt(3); y(1) = 1/sqrt(3); z(1) = 1/sqrt(3);
     x(2) = 1/sqrt(3); y(2) = 1/sqrt(3); z(2) = -1/sqrt(3);
     x(3) = 1/sqrt(3); y(3) = -1/sqrt(3); z(3) = 1/sqrt(3);
     x(4) = 1/sqrt(3); y(4) = -1/sqrt(3); z(4) = -1/sqrt(3);
     x(5) = -1/sqrt(3); y(5) = 1/sqrt(3); z(5) = 1/sqrt(3);
     x(6) = -1/sqrt(3); y(6) = 1/sqrt(3); z(6) = -1/sqrt(3);
     x(7) = -1/sqrt(3); y(7) = -1/sqrt(3); z(7) = 1/sqrt(3);
     x(8) = -1/sqrt(3); y(8) = -1/sqrt(3); z(8) = -1/sqrt(3);
     x(9:12) = 0; 
     y(9:10) = abs(y(9:10)); y(11:12) = -abs(y(11:12));
     z([9 11]) = abs(z([9 11])); z([10 12]) = -abs(z([10 12]));
     z(13:16) = 0;
     x(13:14) = abs(x(13:14)); x(15:16) = -abs(x(15:16));
     y([13 15]) = abs(y([13 15])); y([14 16]) = -abs(y([14 16]));
     y(17:20) = 0;   
     z(17:18) = abs(z(17:18)); z(19:20) = -abs(z(19:20));
     x([17 19]) = abs(x([17 19])); x([18 20]) = -abs(x([18 20]));
   elseif n == 32,
     x(1) = 1; y(1) = 0; z(1) = 0;
     x(2) = -1; y(2) = 0; z(2) = 0;
     if isSeed,
       x(3) = cos(pi/4); y(3) = sin(pi/4); z(3) = 0;
       x(4) = cos(2*pi/4); y(4) = sin(2*pi/4); z(4) = 0;
       x(5) = cos(3*pi/4); y(5) = sin(3*pi/4); z(5) = 0;
       x(6) = cos(5*pi/4); y(6) = sin(5*pi/4); z(6) = 0;
       x(7) = cos(6*pi/4); y(7) = sin(6*pi/4); z(7) = 0;
       x(8) = cos(7*pi/4); y(8) = sin(7*pi/4); z(8) = 0;
     end;% if
     z(9:20) = abs(z(9:20));
     z(21:32) = -abs(z(21:32));
   elseif n == 48,
     x(1) = 1; y(1) = 0; z(1) = 0;
     x(2) = -1; y(2) = 0; z(2) = 0;
     z(1:12) = 0;
     z(13:30) = abs(z(13:30));
     z(31:48) = -abs(z(31:48));   
   elseif n == 64,
     x(1) = 1; y(1) = 0; z(1) = 0;
     x(2) = -1; y(2) = 0; z(2) = 0;
     z(1:20) = 0;
     z(21:42) = abs(z(21:42));
     z(43:64) = -abs(z(43:64));   
   elseif n == 100,
     x(1) = 1; y(1) = 0; z(1) = 0;
     x(2) = -1; y(2) = 0; z(2) = 0;
     z(1:30) = 0;
     z(31:65) = abs(z(31:65));
     z(66:100) = -abs(z(66:100));   
   end;%if

   not_done = 1;

   s = 1;

   % Turns off the divide by 0 warning
   warning off;

   while (not_done),
      iIt = iIt + 1;
      for i = 1:n,

        if n == 12,
          x(1:4) = 0; 
          y(1:2) = abs(y(1:2)); y(3:4) = -abs(y(3:4));
          z([1 3]) = abs(z([1 3])); z([2 4]) = -abs(z([2 4]));
          z(5:8) = 0;
          x(5:6) = abs(x(5:6)); x(7:8) = -abs(x(7:8));
          y([5 7]) = abs(y([5 7])); y([6 8]) = -abs(y([6 8]));
          y(9:12) = 0;   
          z(9:10) = abs(z(9:10)); z(11:12) = -abs(z(11:12));
          x([9 11]) = abs(x([9 11])); x([10 12]) = -abs(x([10 12]));
        elseif n == 20,
          x(1) = 1/sqrt(3); y(1) = 1/sqrt(3); z(1) = 1/sqrt(3);
          x(2) = 1/sqrt(3); y(2) = 1/sqrt(3); z(2) = -1/sqrt(3);
          x(3) = 1/sqrt(3); y(3) = -1/sqrt(3); z(3) = 1/sqrt(3);
          x(4) = 1/sqrt(3); y(4) = -1/sqrt(3); z(4) = -1/sqrt(3);
          x(5) = -1/sqrt(3); y(5) = 1/sqrt(3); z(5) = 1/sqrt(3);
          x(6) = -1/sqrt(3); y(6) = 1/sqrt(3); z(6) = -1/sqrt(3);
          x(7) = -1/sqrt(3); y(7) = -1/sqrt(3); z(7) = 1/sqrt(3);
          x(8) = -1/sqrt(3); y(8) = -1/sqrt(3); z(8) = -1/sqrt(3);
        elseif n == 32,
          x(1) = 1; y(1) = 0; z(1) = 0;
          x(2) = -1; y(2) = 0; z(2) = 0;
          z(3) = 0;
          if (isSeed && i < 3),
            z(2) = 0;
            z(4) = 0;
            z(6) = 0;
            z(8) = 0;
          end;% if
        elseif n == 48,
          x(1) = 1; y(1) = 0; z(1) = 0;
          x(2) = -1; y(2) = 0; z(2) = 0;
          if (isSeed && i < 5),
            z(3) = 0;
            z(4) = 0;
            z(5) = 0;
            z(6) = 0;
          end;% if
        elseif n == 64,
          x(1) = 1; y(1) = 0; z(1) = 0;
          x(2) = -1; y(2) = 0; z(2) = 0;
        elseif n == 100,
          x(1) = 1; y(1) = 0; z(1) = 0;
          x(2) = -1; y(2) = 0; z(2) = 0;
        end;% if
         
         % Calculate the i,j,k vectors for the direction of the repulsive forces.
         ii = x(i) - x;
         jj = y(i) - y;
         kk = z(i) - z;

         rm_new(i,:) = sqrt(ii.^2 + jj.^2 + kk.^2);

         ii = ii./rm_new(i,:);
         jj = jj./rm_new(i,:);
         kk = kk./rm_new(i,:);

         % Take care of the self terms.
         ii(i) = 0;
         jj(i) = 0;
         kk(i) = 0;

         % Use a 1/r^2 repulsive force, but add 0.01 to the denominator to
         % avoid a 0 * Inf below. The self term automatically disappears since
         % the ii,jj,kk vectors were set to zero for self terms.
         f = 1./(0.01 + rm_new(i,:).^2);

         % Sum the forces.
         fi = sum(f.*ii);
         fj = sum(f.*jj);
         fk = sum(f.*kk);

         % Find magnitude
         fn = sqrt(fi.^2 + fj.^2 + fk.^2);

         % Find the unit direction of repulsion.
         fi = fi/fn;
         fj = fj/fn;
         fk = fk/fn;

         % Step a distance s in the direciton of repulsion
         x(i) = x(i) + s.*fi;
         y(i) = y(i) + s.*fj;
         z(i) = z(i) + s.*fk;

         % Scale the coordinates back down to the unit sphere.
         r = sqrt(x(i).^2 + y(i).^2 + z(i).^2);

         x(i) = x(i)/r;
         y(i) = y(i)/r;
         z(i) = z(i)/r;

      end;% for

      if isDemo,

         figure(1);
         cla;
         axis equal;

         [xs,ys,zs] = sphere(20);

         h = surf(xs,ys,zs);
         set(h,'Facecolor',[1 0 0]);
         l = light;
         lighting phong;

         for m = 1:length(x),
             plot3(x(m)*1.01,y(m)*1.01,z(m)*1.01,'*');
             hold on;
         end;

         axis([-1.2 1.2 -1.2 1.2 -1.2 1.2]);

         drawnow;
         pause(0.1);
      end;% if

      % Check convergence
      diff = abs(rm_new - rm_old);

      not_done = any(diff(:) > diffTol);

      rm_old = rm_new;

   end;% while

   % Find the smallest distance between neighboring points. To do this
   % exclude the self terms which are 0.
   tmp = rm_new(:);

   indices = find(tmp~=0);
   
   avgr = min(tmp(indices));

   % Turn back on the default warning state.
   warning backtrace;
   