%% Create 2-d polygons
% row -- number of polygons in a row
% col -- number of polygons in a column
% distance -- distance between adjacent two polygon
% M -- number of vertices (i.e., number of objectives)
function [Polygons] = CreatePolygons(row,col,distance,M)
    scenario = 1;
    % polygon centers
    N = row*col;
    Centers = zeros(N, 2);
    cnt = 1;
    for i = 1:row
        for j = 1:col
            Centers(cnt, :) = [distance*(j-1), distance*(i-1)];
            cnt = cnt + 1;
        end
    end
    
    % polygon vertices
    if scenario == 1
        Angle = 2*pi.*(1:M)/M;
    else
        if mod(M,2) == 0
            Angle = (2.*(1:M)-3).*pi./M;
        else
            Angle = (2.*(1:M)-2).*pi./M;
        end
    end
    
    VertexSet = zeros(length(Angle), 2, N);
    for i=1:N
        vertices = ([sin(Angle)', cos(Angle)']) + repmat(Centers(i, :), size(Angle, 2), 1);
        VertexSet(:, :, i) = vertices;
    end
    Polygons = VertexSet;
end

