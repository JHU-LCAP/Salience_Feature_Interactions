function W_all = draw_weights(weights,Wtype,experiment,showsbj)
% normalize weights for visualization
nrm = 'max'; valself = 2; valcross = 4;
switch Wtype
    case 'GT'
        cols = [0.4660 0.6740 0.1880];
    case 'beh'
        cols = [0 0.4470 0.7410];
    case 'pupil'
        cols = [0.9290 0.6940 0.1250];
    case 'rand'
        cols = [0.2627 0.1255 0.4588];
end

% weights 6x6xsubjects(or folds if applicable)
W = mean(weights,3); % average weights over subjects or flods (if applicable)
W = combine_weights(W); % combine the 6x6 matrix to 3x3, to represent liudness, pitch, and timbre
W = normalize_weights(W,nrm);

% convert weight into visually discrete levels
Lradius = ceil(W(1,1)* valself);
Pradius = ceil(W(2,2)* valself);
Tradius = ceil(W(3,3)* valself);

LPthick = floor(W(1,2)* valcross);
PLthick = floor(W(2,1)* valcross);
LTthick = floor(W(1,3)* valcross);
TLthick = floor(W(3,1)* valcross);
PTthick = floor(W(2,3)* valcross);
TPthick = floor(W(3,2)* valcross);


sizes.Lradius = Lradius;
sizes.Pradius = Pradius;
sizes.Tradius = Tradius;
sizes.LPthick = LPthick;
sizes.PLthick = PLthick;
sizes.LTthick = LTthick;
sizes.TLthick = TLthick;
sizes.PTthick = PTthick;
sizes.TPthick = TPthick;

% Plot weights
% figure
hold on
Px_center = 0; Py_center = 0;
Tx_center = 30; Ty_center = 0;
Lx_center = 15; Ly_center = 26;

% Self-interactions
draw_circle(Px_center,Py_center,Pradius,'P',cols)
draw_circle(Lx_center,Ly_center,Lradius,'L',cols)
draw_circle(Tx_center,Ty_center,Tradius,'T',cols)

% Cross-interactions
draw_arrow(Lx_center, Ly_center, Tx_center, Ty_center, Lradius, Tradius, LTthick);
draw_arrow(Tx_center, Ty_center, Lx_center, Ly_center, Tradius, Lradius, TLthick);
draw_arrow(Px_center, Py_center, Tx_center, Ty_center, Pradius,Tradius,PTthick);
draw_arrow(Tx_center, Ty_center, Px_center, Py_center, Tradius,Pradius,TPthick);
draw_arrow(Px_center, Py_center, Lx_center, Ly_center, Pradius, Lradius, PLthick);
draw_arrow(Lx_center, Ly_center, Px_center, Py_center, Lradius, Pradius, LPthick);

sgtitle(string(experiment)+' '+string(Wtype)+' weights')
axis off

if showsbj
    figure
    Wsbj = weights;
    for sbj = 1:min(size(Wsbj,3),16)
        subplot(4,4,sbj); hold on; axis equal
        W = Wsbj(:,:,sbj);
        W = combine_weights(W);
        W = normalize_weights(W,nrm);

        % convert weight into visually discrete levels
        Lradius = ceil(W(1,1)* valself);
        Pradius = ceil(W(2,2)* valself);
        Tradius = ceil(W(3,3)* valself);

        LPthick = floor(W(1,2)* valcross);
        PLthick = floor(W(2,1)* valcross);
        LTthick = floor(W(1,3)* valcross);
        TLthick = floor(W(3,1)* valcross);
        PTthick = floor(W(2,3)* valcross);
        TPthick = floor(W(3,2)* valcross);

        % Plot weights
        Px_center = 0; Py_center = 0;
        Tx_center = 30; Ty_center = 0;
        Lx_center = 15; Ly_center = 26;

        
        % Self-interactions
        draw_circle(Px_center,Py_center,Pradius,'P',cols)
        draw_circle(Lx_center,Ly_center,Lradius,'L',cols)
        draw_circle(Tx_center,Ty_center,Tradius,'T',cols)

        % Cross-interactions
        draw_arrow(Lx_center, Ly_center, Tx_center, Ty_center, Lradius, Tradius, LTthick);
        draw_arrow(Tx_center, Ty_center, Lx_center, Ly_center, Tradius, Lradius, TLthick);
        draw_arrow(Px_center, Py_center, Tx_center, Ty_center, Pradius,Tradius,PTthick);
        draw_arrow(Tx_center, Ty_center, Px_center, Py_center, Tradius,Pradius,TPthick);
        draw_arrow(Px_center, Py_center, Lx_center, Ly_center, Pradius, Lradius, PLthick);
        draw_arrow(Lx_center, Ly_center, Px_center, Py_center, Lradius, Pradius, LPthick);

        axis off
    end
    sgtitle(string(experiment)+' '+string(Wtype))

end
    

end

function Wcomb = combine_weights(W)
diagindx = {1, 2:3, 4:6};
Wcomb = zeros(3,3);
for i = 1:3
    Wcomb(i,i) = mean((W(diagindx{i},diagindx{i})),"all");
end
Wcomb(1,2) = mean(W(1,2:3),"all");
Wcomb(2,1) = mean(W(2:3,1),"all");
Wcomb(1,3) = mean(W(1,4:6),"all");
Wcomb(3,1) = mean(W(4:6,1),"all");
Wcomb(2,3) = mean(W(2:3,4:6),"all");
Wcomb(3,2) = mean(W(4:6,2:3),"all");
end

function Wnorm = normalize_weights(W,nrm)
Wnorm = [];
switch nrm
    case 'max',
        Wnorm = W./max(W(:));
    case 'sum',
        Wnorm = W./sum(W(:));
end
end

function draw_circle(x,y,r,txt,col)
theta = linspace(0, 2*pi, 100);
xval = r * cos(theta) + x;
yval = r * sin(theta) + y;
fill(xval, yval,col, 'EdgeColor', 'none');
text(x, y, txt, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 10);
end

function draw_arrow(x1, y1, x2, y2, r1, r2, w,col)
if ~exist('col','var'), col = 'k'; end
if w,
    center1 = [x1, y1];
    center2 = [x2, y2];
    vec = (center2 - center1)/norm(center2 - center1);
    vec2 = [-vec(2), vec(1)];
    arrow_offset = 0.2 * r1; % Offset distance for arrows
    touch_point1 = center1 +  r1 * vec + arrow_offset * vec2;
    touch_point2 = center2 - r2 * vec + arrow_offset * vec2;
    quiver(touch_point1(1), touch_point1(2), touch_point2(1) - touch_point1(1), touch_point2(2) - touch_point1(2), 0, 'color',col, 'LineWidth', w, 'MaxHeadSize', 0.4);
end
end