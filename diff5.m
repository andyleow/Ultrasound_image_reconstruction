function dy = diff5(y)
% Differentiation using a five-point stencil:
% LS quadratic fitting on 5 successive points
dy = y;
y = y(:);
dy(1) = [-54 13 40 27 -26]*y(1:5)/70;
dy(2) = [-34 3 20 17 -6]*y(1:5)/70;
dy(end-1) = [6 -17 -20 -3 34]*y(end-4:end)/70;
dy(end) = [26 -27 -40 -13 54]*y(end-4:end)/70;
dy(3:end-2) = conv(y,[2 1 0 -1 -2],'valid')/10;
end