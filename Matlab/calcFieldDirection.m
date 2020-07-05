displayELS;
c = zeros(2,2);
t = 0;
for i= 1:801
    for j = 1:801
        if ~isnan(E(i,j))
            t = t + 1;
            c(1,1) = c(1,1) + xx(i,j)*xx(i,j)*E(i,j);
            c(1,2) = c(1,2) + xx(i,j)*yy(i,j)*E(i,j);
            c(2,1) = c(2,1) + xx(i,j)*yy(i,j)*E(i,j);
            c(2,2) = c(2,2) + yy(i,j)*yy(i,j)*E(i,j);
        end
    end
end
c = c./t;
for ang = 1:180
    R = [cos(ang*pi/180) -sin(ang*pi/180);sin(ang*pi/180) cos(ang*pi/180)];
    new_c =R'*c*R;
    entry(ang) = abs(new_c(1,2));
end
[a,ang] = min(entry)
hold on
plot((-40:40)*cos((90-ang)*pi/180),(-40:40)*sin((90-ang)*pi/180),'k-')