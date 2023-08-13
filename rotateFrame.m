function frm = rotateFrame(crnr, hdng)


% Rotate the shape around its center point
cntr = mean(crnr, 1);
hdng = 360 - hdng;
frm = bsxfun(@minus, crnr, cntr) * [cosd(hdng), sind(hdng); -sind(hdng), cosd(hdng)];
frm(:,1) = frm(:,1) + cntr(1);
frm(:,2) = frm(:,2) + cntr(2);

